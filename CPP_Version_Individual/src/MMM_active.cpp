// By Minseung Kim, Seungwoo Yoon and Seungbum Koo
// KAIST, Daejeon, South Korea
// February 23, 2026

#include "SimUtils.h"
#include "Millard2012EquilibriumMuscleHelper.h"
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

// --- MAIN RUN FUNCTION ---
int run_MMM_active(std::string MMM_type, double L_mn_input) {
    // ==================================================================
    // 0. Select Muscle Model Type
    // ==================================================================
    // std::string MMM_type = "Classic";
    // std::string MMM_type = "DEq";
    // std::string MMM_type = "Rigid";

    std::cout << "Starting MMM Active Simulation (" << MMM_type << ")..." << std::endl;

    try {
        // ------------------------------------------------------------------
        // 1. Configuration
        // ------------------------------------------------------------------
        // double L_mn_input = 1.0;
        double U_input = 0.5;
        double Amp_input = 0.01;
        double Mass_input = 3e9;
        double Damping_input = 0.0;
        double SimTime_input = 100.0;
        double SimDt_input = 0.001;
        double FreqLow_input  = 0.1;
        double FreqHigh_input = 100;
        double NumFreqSamples = 1000;

        double Mass_ext = Mass_input;
        double Damping_ext = Damping_input;
        double SimTime = SimTime_input;
        double dt = SimDt_input;
        double FreqLow = FreqLow_input;
        double FreqHigh = FreqHigh_input;

        // Soleus Muscle Parameters
        double F_mo = 3549.0;
        double L_mo = 0.05;
        double L_ts = 0.25;
        double alphaOpt = 0.4363;
        double V_mmax_norm = 10.0;

        // Initialize Millard2012 Muscle Model (OpenSim Object used for Curves)
        OpenSim::Millard2012EquilibriumMuscle m;
        m.setMaxIsometricForce(F_mo);
        m.setOptimalFiberLength(L_mo);
        m.setTendonSlackLength(L_ts);
        m.setPennationAngleAtOptimalFiberLength(alphaOpt);
        m.setMaxContractionVelocity(V_mmax_norm);

        // Model Branching Parameters
        double dampingBeta = 0.0;
        bool isRigid = false;

        if (MMM_type == "DEq") {
            dampingBeta = 0.1;
        } else if (MMM_type == "Rigid") {
            dampingBeta = 0.1;
            isRigid = true;
        }

        m.setMuscleConfiguration(isRigid, false, dampingBeta);

        // ------------------------------------------------------------------
        // [Curves Definition]
        // ------------------------------------------------------------------
        double t_strainAtOneNormForce    = 0.033;
        double t_stiffnessAtOneNormForce = 1.375 / 0.033;
        double t_normForceAtToeEnd       = 2.0 / 3.0;
        double t_curviness               = 0.5;
        m.setTendonForceLengthCurve(OpenSim::TendonForceLengthCurve(
            t_strainAtOneNormForce, t_stiffnessAtOneNormForce, t_normForceAtToeEnd, t_curviness
        ));

        double a_minActiveNormLength    = 0.4441;
        double a_transitionNormLength   = 0.73;
        double a_maxActiveNormLength    = 1.8123;
        double a_shallowAscendingSlope  = 0.8616;
        double a_minValue               = 0.0;
        m.setActiveForceLengthCurve(OpenSim::ActiveForceLengthCurve(
            a_minActiveNormLength, a_transitionNormLength, a_maxActiveNormLength, a_shallowAscendingSlope, a_minValue
        ));

        double fv_concSlopeAtVmax      = 0.0;
        double fv_concSlopeNearV0      = 0.15;
        double fv_isoSlope             = 5.0;
        double fv_eccSlopeAtVmax       = 0.1;
        double fv_eccSlopeNearV0       = 0.15;
        double fv_maxEccentricForce    = 1.4;
        double fv_concCurviness        = 0.7;
        double fv_eccCurviness         = 0.9;
        m.setForceVelocityCurve(OpenSim::ForceVelocityCurve(
            fv_concSlopeAtVmax, fv_concSlopeNearV0, fv_isoSlope, fv_eccSlopeAtVmax, fv_eccSlopeNearV0, fv_maxEccentricForce, fv_concCurviness, fv_eccCurviness
        ));

        const auto& fvCurve = m.getForceVelocityCurve();
        const auto& flCurve = m.getActiveForceLengthCurve();
        const auto& fpCurve = m.getFiberForceLengthCurve();
        const auto& ftCurve = m.getTendonForceLengthCurve();

        // ------------------------------------------------------------------
        // 2. Initialization (Finding Equilibrium State)
        // ------------------------------------------------------------------
        double L_m_init = L_mo * L_mn_input;
        double L_m_height = L_mo * std::sin(alphaOpt);
        double fl_init = flCurve.calcValue(L_mn_input);
        double fp_init = fpCurve.calcValue(L_mn_input);
        double F_m_init = (U_input * fl_init + fp_init) * F_mo;

        double sin_phi_init = L_m_height / L_m_init;
        if (sin_phi_init > 1.0) sin_phi_init = 1.0;
        double phi_init = std::asin(sin_phi_init);
        double F_t_prev = F_m_init * std::cos(phi_init);

        double L_t_init;
        if (isRigid) {
            L_t_init = L_ts;
        } else {
            L_t_init = calc_tendon_length_from_force(m, F_t_prev);
        }

        double L_mt0 = L_m_init * std::cos(phi_init) + L_t_init;
        double F_ext_equil = F_t_prev;

        std::cout << "Initialization Done. L_t_init = " << L_t_init << std::endl;

        // ------------------------------------------------------------------
        // 3. Frequency Sweep Loop
        // ------------------------------------------------------------------
        std::vector<double> frequencies = logspace(std::log10(FreqLow), std::log10(FreqHigh), NumFreqSamples);

        std::string dir_name = "MMM_" + MMM_type + "_Active_csv_" + std::to_string(L_mn_input);
        std::filesystem::create_directories(dir_name);

        std::ofstream freqFile(dir_name + "/frequency_list.csv");
        if (freqFile.is_open()) {
            for (size_t k = 0; k < frequencies.size(); ++k) {
                freqFile << k << "," << frequencies[k] << "\n";
            }
            freqFile.close();
        }

        for (int k = 0; k < NumFreqSamples; ++k) {
            double freq = frequencies[k];
            std::ofstream outFile(dir_name + "/freq_res_" + std::to_string(k) + ".csv");
            outFile << "time,u,F_m_AT\n";

            double L_mt = L_mt0;
            double V_mt = 0.0;
            double a = U_input;
            double L_m = L_m_init;
            int time_len = (int)(SimTime / dt);

            // Time Integration Loop
            for (int i = 0; i < time_len; ++i) {
                double t = i * dt;

                // 1. Input: Sinusoidal Excitation (u)
                double u = sinwave(freq, t, U_input, U_input, Amp_input);
                if (u < 0) u = 0; if (u > 1) u = 1;

                // 2. Activation Dynamics
                if (i > 0) a += dt * ((u - a) / get_Tau_MMM(a, u));

                // 3. Muscle Dynamics (Branching by Tendon Type)
                double L_m_prev = L_m;
                double V_mt_prev = V_mt;
                double F_t_final = 0.0;

                if (!isRigid) {
                    // ==========================================
                    // ELASTIC TENDON (Classic, DEq) - Newton Raphson
                    // ==========================================
                    double L_m_curr = L_m;
                    int max_iter = 50;
                    double tol = 1e-8;
                    double delta = 1e-7;
                    double error = 1.0;
                    int iter = 0;

                    while (std::abs(error) > tol && iter < max_iter) {
                        iter++;
                        auto calc_force_error = [&](double L_m_val) -> double {
                            double sin_phi = L_m_height / L_m_val;
                            if (sin_phi > 1.0) sin_phi = 1.0;
                            double cos_phi = std::sqrt(1.0 - sin_phi * sin_phi);

                            double L_t_val = L_mt - L_m_val * cos_phi;
                            double L_tn_val = L_t_val / L_ts;
                            double F_t_val = ftCurve.calcValue(L_tn_val) * F_mo;

                            double V_m_val = (L_m_val - L_m_prev) / dt;
                            double V_mn_val = V_m_val / (V_mmax_norm * L_mo);

                            double L_mn_val = L_m_val / L_mo;
                            double fl = flCurve.calcValue(L_mn_val);
                            double fp = fpCurve.calcValue(L_mn_val);
                            double fv = fvCurve.calcValue(L_mn_val);

                            double F_mn_val = a * fl * fv + fp + dampingBeta * V_mn_val;
                            double F_m_AT_val = (F_mn_val * F_mo) * cos_phi;
                            return F_t_val - F_m_AT_val;
                        };

                        double err1 = calc_force_error(L_m_curr);
                        double err2 = calc_force_error(L_m_curr + delta);
                        double J = (err2 - err1) / delta;
                        if (std::abs(J) < 1e-14) J = 1e-14;

                        L_m_curr -= err1 / J;

                        if (L_m_curr < 0.001) L_m_curr = 0.001;
                        if (L_m_curr > L_mt)  L_m_curr = L_mt - 0.001;
                        error = err1;
                    }

                    L_m = L_m_curr;
                    double sin_phi = L_m_height / L_m;
                    if(sin_phi > 1.0) sin_phi = 1.0;
                    double cos_phi = std::sqrt(1.0 - sin_phi * sin_phi);

                    double L_t = L_mt - L_m * cos_phi;
                    double L_tn = L_t / L_ts;
                    F_t_final = ftCurve.calcValue(L_tn) * F_mo;

                } else {
                    // ==========================================
                    // RIGID TENDON (Rigid) - Direct Kinematic Calc
                    // ==========================================
                    double L_m_AT = L_mt - L_ts;
                    if (L_m_AT < 0.001) L_m_AT = 0.001;

                    L_m = std::sqrt(L_m_AT * L_m_AT + L_m_height * L_m_height);

                    double V_m_val = (L_m - L_m_prev) / dt;
                    double V_mn = V_m_val / (V_mmax_norm * L_mo);

                    double sin_phi = L_m_height / L_m;
                    if (sin_phi > 1.0) sin_phi = 1.0;
                    double cos_phi = std::sqrt(1.0 - sin_phi * sin_phi);

                    double L_mn = L_m / L_mo;
                    double fl = flCurve.calcValue(L_mn);
                    double fp = fpCurve.calcValue(L_mn);
                    double fv = fvCurve.calcValue(V_mn);

                    double F_fiber_norm = a * fl * fv + fp + dampingBeta * V_mn;
                    F_t_final = (F_fiber_norm * F_mo) * cos_phi;

                    if (F_t_final < 0.0) F_t_final = 0.0;
                }

                // 4. External Dynamics
                double A_mt = (F_ext_equil - F_t_final - V_mt_prev * Damping_ext) / Mass_ext;
                V_mt += A_mt * dt;
                L_mt += V_mt * dt;

                // Output Data
                outFile << t << "," << u << "," << F_t_final << "\n";
            }
            outFile.close();
            if (k % 10 == 0) std::cout << "Frequency step " << k << " done." << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }
    return 0;
}