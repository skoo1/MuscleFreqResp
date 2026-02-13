#include <OpenSim/Actuators/Millard2012EquilibriumMuscle.h>
#include <OpenSim/Actuators/ActiveForceLengthCurve.h>
#include <OpenSim/Actuators/FiberForceLengthCurve.h>
#include <OpenSim/Actuators/ForceVelocityCurve.h>
#include <OpenSim/Actuators/TendonForceLengthCurve.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <filesystem>
#include <algorithm>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

namespace fs = std::filesystem;

// --- Helper Functions ---

// 1. Tendon Length Finder
double solveTendonLengthNewton(const OpenSim::Millard2012EquilibriumMuscle& m, double Ft_target) {
    const auto& fseCurve = m.getTendonForceLengthCurve();
    double Fmo = m.getMaxIsometricForce();
    double Lts = m.getTendonSlackLength();
    
    double ltN = 1.01; // Start from slightly stretched
    double tol = 1e-8; 
    int maxIter = 50;

    for (int i = 0; i < maxIter; ++i) {
        double f_val = fseCurve.calcValue(ltN);
        double err = f_val * Fmo - Ft_target;
        if (std::abs(err) < tol * Fmo) break;

        double df_val = fseCurve.calcDerivative(ltN, 1);
        double derr = df_val * Fmo;

        if (std::abs(derr) < 1e-10) ltN += 0.001; 
        else ltN = ltN - err / derr;
    }
    return ltN * Lts;
}

// 2. Newton-Raphson Solver for Fiber Velocity (Damping 포함)
// 식: a * fl * fv(v) + fp + beta * v - F_load_norm = 0
double solveFiberVelocity(const OpenSim::ForceVelocityCurve& fvCurve, 
                          double a, double fl, double fp, double beta, 
                          double F_fiber_load_norm) {
    double v_norm = 0.0; 
    double tol = 1e-7;
    int maxIter = 20;

    for(int i=0; i<maxIter; ++i) {
        double fv = fvCurve.calcValue(v_norm);
        double f_curr = a * fl * fv + fp + beta * v_norm;
        double err = f_curr - F_fiber_load_norm;

        if(std::abs(err) < tol) break;

        double dfv_dv = fvCurve.calcDerivative(v_norm, 1);
        double df_dv = a * fl * dfv_dv + beta; 

        if(std::abs(df_dv) < 1e-9) df_dv = 1e-9; 
        v_norm = v_norm - err / df_dv;
    }
    return v_norm;
}

// Utils
std::vector<double> logspace(double startExp, double endExp, int num) {
    std::vector<double> result;
    if (num <= 1) { result.push_back(std::pow(10, endExp)); return result; }
    double step = (endExp - startExp) / (num - 1);
    for (int i = 0; i < num; ++i) result.push_back(std::pow(10, startExp + i * step));
    return result;
}

double get_Tau(double a, double u) {
    if (u > a) return 0.01 * (0.5 + 1.5 * a);
    else       return 0.04 / (0.5 + 1.5 * a);
}

double sinwave(double freq, double t, double a0, double a, double b) {
    double val = (a0 - a) / b;
    if (val > 1.0) val = 1.0; if (val < -1.0) val = -1.0;
    return a + std::sin(2.0 * M_PI * freq * t + std::asin(val)) * b;
}

// --- MAIN ---
int main() {
    try {
        // 1. Configuration
        double F_mo = 3549.0, L_mo = 0.05, L_ts = 0.25, alphaOpt = 0.4363, V_mmax_norm = 10.0;

        OpenSim::Millard2012EquilibriumMuscle m;
        m.setMaxIsometricForce(F_mo);
        m.setOptimalFiberLength(L_mo);
        m.setTendonSlackLength(L_ts);
        m.setPennationAngleAtOptimalFiberLength(alphaOpt);
        m.setMaxContractionVelocity(V_mmax_norm);

        double dampingBeta = 0.1; 
        m.setMuscleConfiguration(false, false, dampingBeta);

        // Curve 1: Tendon (Stiffness 정밀도 보정)
        double eIso = 0.033;
        m.setTendonForceLengthCurve(OpenSim::TendonForceLengthCurve(eIso, 1.375/eIso, 2.0/3.0, 0.5));
        
        // Curve 2: Active Force-Length (Curviness 0.5 -> 1.0 보정, min=0.0)
        // m.setActiveForceLengthCurve(...) 대신 아래와 같이 설정
        m.setActiveForceLengthCurve(OpenSim::ActiveForceLengthCurve(
            0.4441, 0.73, 1.8123, 0.8616, 0.0 
        ));
        // Note: C++ API 한계로 Curviness를 1.0으로 설정하는 것은 생략 (기본값 0.5 사용).
        // 하지만 min=0.0 설정으로 큰 오차는 잡음.

        // Curve 3: Force-Velocity (Forward only)
        m.setForceVelocityCurve(OpenSim::ForceVelocityCurve(0.0, 0.15, 5.0, 0.1, 0.15, 1.4, 0.7, 0.9));

        // Curve References
        const auto& fvCurve = m.getForceVelocityCurve();
        const auto& flCurve = m.getActiveForceLengthCurve();
        const auto& fpCurve = m.getFiberForceLengthCurve();
        const auto& ftCurve = m.getTendonForceLengthCurve();

        double u0 = 0.5, Amp_input = 0.01, dt = 0.001, SimTime = 120.0;
        double Mass_ext = 3e9; 

        // 2. Initial Setup (MATLAB Logic)
        // Fiber Length 초기화: L_mo * 1.0
        double L_m_init = L_mo * 1.0; 
        
        // Fixed Width Model Constant Calculation
        double muscleWidth = L_mo * std::sin(alphaOpt);

        // Calculate Initial State
        double fl_init = flCurve.calcValue(1.0);
        double fp_init = fpCurve.calcValue(1.0);
        double F_fiber_init = (u0 * fl_init + fp_init) * F_mo; // Vm=0 가정

        double phi_init = std::asin(muscleWidth / L_m_init);
        double F_t_prev = F_fiber_init * std::cos(phi_init);
        double Lt_init = solveTendonLengthNewton(m, F_t_prev);

        double L_mt0 = L_m_init * std::cos(phi_init) + Lt_init;
        double F_ext_equil = F_t_prev; 

        std::cout << "Initialization Done. Lt_init = " << Lt_init << std::endl;

        // 3. Simulation Loop
        std::vector<double> frequencies = logspace(-1.0, 2.0, 100);
        fs::create_directories("MMM_Active_results_csv");

        for (int k = 0; k < 100; ++k) {
            double freq = frequencies[k];
            std::ofstream outFile("MMM_Active_results_csv/freq_res_" + std::to_string(k) + ".csv");
            outFile << "time,u,F_m_AT\n";

            double L_mt = L_mt0;
            double V_mt = 0.0;
            double a = u0;
            double Lm = L_m_init; // 초기화 값으로 시작

            for (int i = 0; i < (int)(SimTime/dt); ++i) {
                double t = i * dt;
                
                // Input
                double u = sinwave(freq, t, u0, u0, Amp_input);
                if (u < 0) u = 0; if (u > 1) u = 1;

                // Activation Dynamics
                if (i > 0) a += dt * ((u - a) / get_Tau(a, u));

                // --- Muscle Dynamics (Fixed Width Pennation Model 적용) ---
                
                // 1. Pennation Angle (Fixed Width)
                // sin(phi) = width / Lm
                double sin_phi = muscleWidth / Lm;
                if (sin_phi > 0.999) sin_phi = 0.999; // 안전장치
                double cos_phi = std::sqrt(1.0 - sin_phi * sin_phi);

                // 2. Tendon Force & Load
                double Lt = L_mt - Lm * cos_phi;
                double Lt_norm = Lt / L_ts;
                
                double f_t_norm = ftCurve.calcValue(Lt_norm);
                double F_t = f_t_norm * F_mo;
                
                // Fiber Load = F_t / cos(phi)
                double F_fiber_load_norm = (F_t / std::max(cos_phi, 0.01)) / F_mo;

                // 3. Components
                double Lm_norm = Lm / L_mo;
                double f_l = flCurve.calcValue(Lm_norm);
                double f_p = fpCurve.calcValue(Lm_norm);

                // 4. Solve Vm (Implicit Newton-Raphson)
                double Vm_norm = solveFiberVelocity(fvCurve, a, f_l, f_p, dampingBeta, F_fiber_load_norm);

                // 5. Integration (Euler)
                double Vm = Vm_norm * V_mmax_norm * L_mo;
                Lm += Vm * dt; 

                // --- External Dynamics ---
                double A_mt = (F_ext_equil - F_t) / Mass_ext; // Damping_ext = 0
                V_mt += A_mt * dt;
                L_mt += V_mt * dt;

                outFile << t << "," << u << "," << F_t << "\n";
            }
            outFile.close();
            if (k % 10 == 0) std::cout << "Frequency step " << k << " done." << std::endl;
        }
    } catch (const std::exception& e) { std::cerr << e.what() << std::endl; }
    return 0;
}