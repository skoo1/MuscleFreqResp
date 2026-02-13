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
#include <iomanip> 

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

namespace fs = std::filesystem;

// --- 구조체: 계산 결과 반환 ---
struct MuscleInfo {
    double F_m_AT; // Tendon 방향으로 투영된 Fiber 힘
    double F_t;    // Tendon 힘
};

// --- 핵심 함수: 근육 상태 계산 ---
MuscleInfo calcMuscleInfo(
    double L_mt,        // Musculotendon Length
    double L_m_AT,      // Fiber Length along Tendon
    double V_m_AT,      // Fiber Velocity along Tendon
    double a,           // Activation
    double muscleWidth, // Fixed Width (Constant)
    double F_mo, double L_mo, double L_ts, double V_mmax, double dampingBeta,
    const OpenSim::ActiveForceLengthCurve& flCurve,
    const OpenSim::FiberForceLengthCurve& fpCurve,
    const OpenSim::ForceVelocityCurve& fvCurve,
    const OpenSim::TendonForceLengthCurve& ftCurve
) {
    // 1. Geometry (Fixed Width Model)
    // Lm^2 = Lm_AT^2 + width^2
    double Lm = std::sqrt(L_m_AT * L_m_AT + muscleWidth * muscleWidth);
    double cos_phi = L_m_AT / Lm;
    
    // Kinematics: Vm = Vm_AT * cos(phi)
    double Vm = V_m_AT * cos_phi;

    // 2. Tendon Force
    double Lt = L_mt - L_m_AT;
    double Lt_norm = Lt / L_ts;
    double F_t = ftCurve.calcValue(Lt_norm) * F_mo;

    // 3. Fiber Force Components
    double Lm_norm = Lm / L_mo;
    double Vm_norm = Vm / (L_mo * V_mmax);

    double fl = flCurve.calcValue(Lm_norm);
    double fp = fpCurve.calcValue(Lm_norm);
    double fv = fvCurve.calcValue(Vm_norm);

    // Total Fiber Force
    // [수정 포인트] dampingBeta가 0이면 감쇠력은 사라짐
    double F_fiber_norm = (a * fl * fv + fp + dampingBeta * Vm_norm);
    double F_fiber = F_fiber_norm * F_mo;

    // 4. Project to Tendon
    double F_m_AT = F_fiber * cos_phi;

    return {F_m_AT, F_t};
}

// Utils
std::vector<double> logspace(double startExp, double endExp, int num) {
    std::vector<double> result;
    if (num <= 1) { result.push_back(std::pow(10, endExp)); return result; }
    double step = (endExp - startExp) / (num - 1);
    for (int i = 0; i < num; ++i) result.push_back(std::pow(10, startExp + i * step));
    return result;
}

int main() {
    try {
        // ------------------------------------------------------------------
        // 1. Configuration (MATLAB과 동일)
        // ------------------------------------------------------------------
        double L_mtn_input = 1.05;
        double U_input = 1.49012e-08; 
        double Amp_input = 0.005;
        double SimTime_input = 120.0;
        double SimDt_input = 0.001;
        double FreqLow_input = 0.1;
        double FreqHigh_input = 100.0;
        int NumFreqSamples = 100;

        // Muscle Properties
        double F_mo = 3549.0;
        double L_mo = 0.05;
        double L_ts = 0.25;
        double alphaOpt = 0.4363;
        double V_mmax_norm = 10.0;

        // ------------------------------------------------------------------
        // 2. OpenSim Muscle Object Setup
        // ------------------------------------------------------------------
        OpenSim::Millard2012EquilibriumMuscle m;
        m.setMaxIsometricForce(F_mo);
        m.setOptimalFiberLength(L_mo);
        m.setTendonSlackLength(L_ts);
        m.setPennationAngleAtOptimalFiberLength(alphaOpt);
        m.setMaxContractionVelocity(V_mmax_norm);

        // [핵심 수정] MATLAB의 useFiberDamping = 0 설정 반영
        // Damping을 0.0으로 설정하여 Phase Shift를 제거함
        double dampingBeta = 0.0; 
        
        // [Curve 1] Tendon Curve (MATLAB 파라미터)
        double eIso = 0.033;
        double kIso = 1.375/eIso;
        m.setTendonForceLengthCurve(OpenSim::TendonForceLengthCurve(eIso, kIso, 2.0/3.0, 0.5));
        
        // [Curve 2] Active FL Curve (Min=0.0 필수)
        m.setActiveForceLengthCurve(OpenSim::ActiveForceLengthCurve(0.4441, 0.73, 1.8123, 0.8616, 0.0));
        
        // [Curve 3] Passive FL Curve (MATLAB 파라미터 강제 적용)
        m.setFiberForceLengthCurve(OpenSim::FiberForceLengthCurve(0.0, 0.7, 0.2, 2.0/0.7, 0.75));

        // [Curve 4] FV Curve
        m.setForceVelocityCurve(OpenSim::ForceVelocityCurve(0.0, 0.15, 5.0, 0.1, 0.15, 1.4, 0.7, 0.9));

        // Get References
        const auto& fvCurve = m.getForceVelocityCurve();
        const auto& flCurve = m.getActiveForceLengthCurve();
        const auto& fpCurve = m.getFiberForceLengthCurve();
        const auto& ftCurve = m.getTendonForceLengthCurve();

        // ------------------------------------------------------------------
        // 3. Initialization
        // ------------------------------------------------------------------
        double L_mt0_ref = L_ts + L_mo * std::cos(alphaOpt);
        double L_mt_target_val = L_mt0_ref * L_mtn_input;
        double muscleWidth = L_mo * std::sin(alphaOpt);

        double low_L_m_AT = 0.0;
        double high_L_m_AT = L_mt_target_val - 1e-5;
        double L_m_AT_init = 0.0;
        double F_t_init = 0.0;

        double err = 1000.0;
        int init_iter = 0;
        while (std::abs(err) > 1e-8 && init_iter < 1000) {
            L_m_AT_init = 0.5 * (low_L_m_AT + high_L_m_AT);
            auto info = calcMuscleInfo(L_mt_target_val, L_m_AT_init, 0.0, U_input, 
                                       muscleWidth, F_mo, L_mo, L_ts, V_mmax_norm, dampingBeta,
                                       flCurve, fpCurve, fvCurve, ftCurve);
            err = info.F_t - info.F_m_AT;
            F_t_init = info.F_t;
            if (err < 0) high_L_m_AT = L_m_AT_init;
            else         low_L_m_AT = L_m_AT_init;
            init_iter++;
        }
        
        std::cout << "Init Done. F_t_init: " << std::setprecision(12) << F_t_init << std::endl;

        // ------------------------------------------------------------------
        // 4. Frequency Sweep Loop
        // ------------------------------------------------------------------
        std::vector<double> frequencies = logspace(std::log10(FreqLow_input), std::log10(FreqHigh_input), NumFreqSamples);
        fs::create_directories("MMM_Passive_results_csv");

        double amp_val = Amp_input * L_mt0_ref * L_mtn_input; 

        for (int k = 0; k < NumFreqSamples; ++k) {
            double freq = frequencies[k];
            std::ofstream outFile("MMM_Passive_results_csv/freq_res_" + std::to_string(k) + ".csv");
            outFile << std::setprecision(16);
            outFile << "time,L_mt,F_m_AT\n";

            double L_m_AT = L_m_AT_init;
            int time_len = (int)(SimTime_input / SimDt_input); 

            for (int i = 0; i < time_len; ++i) {
                double t = i * SimDt_input;
                
                // Motion Inputs
                double L_mt = L_mt_target_val + amp_val * std::sin(2 * M_PI * freq * t);
                
                double L_m_AT_old = L_m_AT;

                if (i > 0) {
                    int count = 0;
                    int max_iter = 50;
                    double error1 = 1.0;
                    double delta = 1e-7;

                    while (std::abs(error1) > 1e-8 && count < max_iter) {
                        count++;
                        
                        double V_m_AT = (L_m_AT - L_m_AT_old) / SimDt_input;
                        auto info1 = calcMuscleInfo(L_mt, L_m_AT, V_m_AT, U_input, 
                                                   muscleWidth, F_mo, L_mo, L_ts, V_mmax_norm, dampingBeta,
                                                   flCurve, fpCurve, fvCurve, ftCurve);
                        error1 = info1.F_t - info1.F_m_AT;
                        if (std::abs(error1) <= 1e-8) break;

                        double L_m_AT_perturb = L_m_AT + delta;
                        double V_m_AT_perturb = (L_m_AT_perturb - L_m_AT_old) / SimDt_input;
                        auto info2 = calcMuscleInfo(L_mt, L_m_AT_perturb, V_m_AT_perturb, U_input, 
                                                   muscleWidth, F_mo, L_mo, L_ts, V_mmax_norm, dampingBeta,
                                                   flCurve, fpCurve, fvCurve, ftCurve);
                        double error2 = info2.F_t - info2.F_m_AT;

                        double J = (error2 - error1) / delta;
                        if (std::abs(J) < 1e-14) break;
                        L_m_AT -= error1 / J;

                        if (L_m_AT < 1e-6) L_m_AT = 1e-6;
                        if (L_m_AT > L_mt) L_m_AT = L_mt - 1e-6;
                    }
                }

                double V_m_AT_final = 0.0;
                if (i > 0) V_m_AT_final = (L_m_AT - L_m_AT_old) / SimDt_input;
                
                auto infoEq = calcMuscleInfo(L_mt, L_m_AT, V_m_AT_final, U_input, 
                                            muscleWidth, F_mo, L_mo, L_ts, V_mmax_norm, dampingBeta,
                                            flCurve, fpCurve, fvCurve, ftCurve);

                outFile << t << "," << L_mt << "," << infoEq.F_t << "\n";
            }
            outFile.close();
            if (k % 10 == 0) std::cout << "Frequency step " << k << " done." << std::endl;
        }
    } catch (const std::exception& e) { std::cerr << e.what() << std::endl; }
    return 0;
}