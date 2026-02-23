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
#include <algorithm>
#include <iomanip>

int run_MMM_passive(std::string MMM_type) {
    // ==================================================================
    // 0. Select Muscle Model Type
    // ==================================================================
    // std::string MMM_type = "Classic"; 
    // std::string MMM_type = "DEq";
    // std::string MMM_type = "Rigid";

    std::cout << "Starting MMM Passive Simulation (" << MMM_type << ")..." << std::endl;
    
    try {
        // ------------------------------------------------------------------
        // 1. Configuration (Parameters matching the MATLAB simulation)
        // ------------------------------------------------------------------
        double L_mtn_input = 1.05;      // Target Normalized Musculotendon Length multiplier
        double U_input = 0.0;           // Input Excitation (Very small value effectively 0 for passive)
        double Amp_input = 0.005;       // Amplitude of the oscillation (normalized)
        double SimTime_input = 120.0;   // Total simulation time (seconds)
        double SimDt_input = 0.001;     // Simulation time step (seconds)
        double FreqLow_input = 0.1;     // Start frequency for sweep (Hz)
        double FreqHigh_input = 100.0;  // End frequency for sweep (Hz)
        int NumFreqSamples = 100;       // Number of frequency steps

        // Frequency Sweep Settings
        double SimTime = SimTime_input;
        double dt = SimDt_input;
        double FreqLow = FreqLow_input;          
        double FreqHigh = FreqHigh_input;      

        // Muscle Properties (Soleus muscle parameters)
        double F_mo = 3549.0;           
        double L_mo = 0.05;             
        double L_ts = 0.25;             
        double alphaOpt = 0.4363;       
        double V_mmax_norm = 10.0;      
        
        // Model Branching Parameters
        double dampingBeta = 0.0;
        bool isRigid = false;

        if (MMM_type == "DEq") {
            dampingBeta = 0.1;
        } else if (MMM_type == "Rigid") {
            dampingBeta = 0.1;
            isRigid = true;
        }

        // ------------------------------------------------------------------
        // 2. OpenSim Muscle Object Setup
        // ------------------------------------------------------------------
        OpenSim::Millard2012EquilibriumMuscle m;
        m.setMaxIsometricForce(F_mo);
        m.setOptimalFiberLength(L_mo);
        m.setTendonSlackLength(L_ts);
        m.setPennationAngleAtOptimalFiberLength(alphaOpt);
        m.setMaxContractionVelocity(V_mmax_norm);
        m.setMuscleConfiguration(isRigid, false, dampingBeta);

        // -------------------------------------------------------------------------
        // [Curve 1] Tendon Force-Length Curve
        // -------------------------------------------------------------------------
        double t_strainAtOneNormForce    = 0.033;             // e0
        double t_stiffnessAtOneNormForce = 1.375 / 0.033;     // k_iso
        double t_normForceAtToeEnd       = 2.0 / 3.0;         // F_toe
        double t_curviness               = 0.5;               
        m.setTendonForceLengthCurve(OpenSim::TendonForceLengthCurve(
            t_strainAtOneNormForce, t_stiffnessAtOneNormForce, t_normForceAtToeEnd, t_curviness
        ));

        // -------------------------------------------------------------------------
        // [Curve 2] Active Force-Length Curve
        // -------------------------------------------------------------------------
        double a_minActiveNormLength    = 0.4441;
        double a_transitionNormLength   = 0.73;
        double a_maxActiveNormLength    = 1.8123;
        double a_shallowAscendingSlope  = 0.8616;
        double a_minValue               = 0.0;
        m.setActiveForceLengthCurve(OpenSim::ActiveForceLengthCurve(
            a_minActiveNormLength, a_transitionNormLength, a_maxActiveNormLength, a_shallowAscendingSlope, a_minValue
        ));

        // -------------------------------------------------------------------------
        // [Curve 3] Passive Force-Length Curve
        // -------------------------------------------------------------------------
        double p_strainAtZeroForce       = 0.0;
        double p_strainAtOneNormForce    = 0.7;          
        double p_stiffnessAtLowForce     = 0.2;          
        double p_stiffnessAtOneNormForce = 2.0 / 0.7;    
        double p_curviness               = 0.75;
        m.setFiberForceLengthCurve(OpenSim::FiberForceLengthCurve(
            p_strainAtZeroForce, p_strainAtOneNormForce, p_stiffnessAtLowForce, p_stiffnessAtOneNormForce, p_curviness
        ));

        // -------------------------------------------------------------------------
        // [Curve 4] Force-Velocity Curve
        // -------------------------------------------------------------------------
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
        // 3. Initialization (Branching for Rigid Tendon)
        // ------------------------------------------------------------------
        double L_m_height = L_mo * std::sin(alphaOpt);
        double L_mt0_ref = L_ts + L_mo * std::cos(alphaOpt);
        double L_mt_target = L_mt0_ref * L_mtn_input;

        double L_m_AT_init = 0.0;

        if (isRigid) {
            L_m_AT_init = L_mt_target - L_ts;
        } else {
            double low_L_m_AT = 0.0;
            double high_L_m_AT = L_mt_target - 1e-5;

            for(int iter=0; iter<100; ++iter) {
                L_m_AT_init = 0.5 * (low_L_m_AT + high_L_m_AT);

                double L_m = std::sqrt(L_m_AT_init * L_m_AT_init + L_m_height * L_m_height);
                double cos_phi = L_m_AT_init / L_m;
                double L_t = L_mt_target - L_m_AT_init;

                double f_t_norm = ftCurve.calcValue(L_t / L_ts);
                double f_p_norm = fpCurve.calcValue(L_m / L_mo);

                double err = (f_t_norm * F_mo) - (f_p_norm * F_mo * cos_phi);

                if (err < 0) high_L_m_AT = L_m_AT_init;
                else         low_L_m_AT = L_m_AT_init;

                if(std::abs(err) < 1e-5) break;
            }
        }
        std::cout << "Init Done. L_m_AT_init = " << L_m_AT_init << std::endl;

        // ------------------------------------------------------------------
        // 4. Frequency Sweep Loop
        // ------------------------------------------------------------------
        std::vector<double> frequencies = logspace(std::log10(FreqLow_input), std::log10(FreqHigh_input), NumFreqSamples);

        std::string res_folder = "MMM_" + MMM_type + "_Passive_csv";
        std::filesystem::create_directories(res_folder);

        std::ofstream freqFile(res_folder + "/frequency_list.csv");
        if (freqFile.is_open()) {
            for (size_t k = 0; k < frequencies.size(); ++k) {
                freqFile << k << "," << frequencies[k] << "\n";
            }
            freqFile.close();
        }

        double amp_val = Amp_input * L_mt0_ref * L_mtn_input;

        for (int k = 0; k < NumFreqSamples; ++k) {
            double freq = frequencies[k];

            std::ofstream outFile(res_folder + "/freq_res_" + std::to_string(k) + ".csv");
            outFile << std::setprecision(16);
            outFile << "time,L_mt,F_m_AT\n";

            double L_m_AT = L_m_AT_init; 
            double a = U_input;
            int time_len = (int)(SimTime / dt);

            // Time Loop for Simulation
            for (int i = 0; i < time_len; ++i) {
                double t = i * SimDt_input;

                // 1. Calculate Input: L_mt oscillation (Sine wave)
                double L_mt = L_mt_target + amp_val * std::sin(2.0 * M_PI * freq * t);

                double L_m_prev_sq = L_m_AT * L_m_AT + L_m_height * L_m_height;
                double L_m_prev_val = std::sqrt(L_m_prev_sq);
                double F_t_final = 0.0;

                // 2. Muscle Dynamics (Branching by Tendon Type)
                if (!isRigid) {
                    // ==========================================
                    // ELASTIC TENDON (Classic, DEq) - Newton Raphson
                    // ==========================================
                    if (i > 0) {
                        double L_m_AT_curr = L_m_AT; 
                        int max_iter = 50;
                        double tol = 1e-8;
                        double delta = 1e-7;
                        double error = 1.0;
                        int iter = 0;

                        while (std::abs(error) > tol && iter < max_iter) {
                            iter++;
                            auto calc_force_error = [&](double l_at_val) -> double {
                                double L_m_curr_val = std::sqrt(l_at_val * l_at_val + L_m_height * L_m_height);
                                double cos_phi = l_at_val / L_m_curr_val;
                                
                                double L_t = L_mt - l_at_val;
                                double f_t_norm = ftCurve.calcValue(L_t / L_ts);
                                double F_t = f_t_norm * F_mo;

                                double V_m = (L_m_curr_val - L_m_prev_val) / SimDt_input;
                                double V_m_norm = V_m / (V_mmax_norm * L_mo);

                                double L_m_norm = L_m_curr_val / L_mo;
                                double f_l = flCurve.calcValue(L_m_norm);
                                double f_p = fpCurve.calcValue(L_m_norm);
                                double f_v = fvCurve.calcValue(V_m_norm);

                                double F_fib_norm = a * f_l * f_v + f_p + dampingBeta * V_m_norm;
                                double F_fib_proj = (F_fib_norm * F_mo) * cos_phi;

                                return F_t - F_fib_proj;
                            };

                            double err1 = calc_force_error(L_m_AT_curr);
                            double err2 = calc_force_error(L_m_AT_curr + delta);
                            
                            double J = (err2 - err1) / delta;
                            if (std::abs(J) < 1e-14) J = 1e-14;

                            L_m_AT_curr -= err1 / J;

                            if (L_m_AT_curr < 1e-6) L_m_AT_curr = 1e-6;
                            if (L_m_AT_curr > L_mt) L_m_AT_curr = L_mt - 1e-6;

                            error = err1;
                        }
                        L_m_AT = L_m_AT_curr; 
                    }
                    
                    double L_t_final = L_mt - L_m_AT;
                    F_t_final = ftCurve.calcValue(L_t_final / L_ts) * F_mo;

                } else {
                    // ==========================================
                    // RIGID TENDON (Rigid) - Direct Kinematic Calc
                    // ==========================================
                    L_m_AT = L_mt - L_ts;
                    if (L_m_AT < 0.001) L_m_AT = 0.001; 

                    double L_m_curr_val = std::sqrt(L_m_AT * L_m_AT + L_m_height * L_m_height);
                    double cos_phi = L_m_AT / L_m_curr_val;

                    double V_m = 0.0;
                    if (i > 0) {
                        V_m = (L_m_curr_val - L_m_prev_val) / SimDt_input;
                    }
                    double V_m_norm = V_m / (V_mmax_norm * L_mo);
                    
                    double L_m_norm = L_m_curr_val / L_mo;
                    double f_l = flCurve.calcValue(L_m_norm);
                    double f_p = fpCurve.calcValue(L_m_norm);
                    double f_v = fvCurve.calcValue(V_m_norm);

                    double F_fib_norm = a * f_l * f_v + f_p + dampingBeta * V_m_norm;
                    F_t_final = (F_fib_norm * F_mo) * cos_phi;

                    if (F_t_final < 0.0) F_t_final = 0.0;
                }

                // 3. Output Data Writing
                outFile << t << "," << L_mt << "," << F_t_final << "\n";
            }
            outFile.close();
            if (k % 10 == 0) std::cout << "Frequency step " << k << " done." << std::endl;
        }
    } catch (const std::exception& e) { 
        std::cerr << e.what() << std::endl; 
    }
    return 0;
}