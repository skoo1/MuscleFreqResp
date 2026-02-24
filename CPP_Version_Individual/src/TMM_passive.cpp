// By Minseung Kim, Seungwoo Yoon and Seungbum Koo
// KAIST, Daejeon, South Korea
// February 23, 2026

#include "SimUtils.h"
#include "Thelen2003MuscleWrapper.h"
#include "Thelen2003MuscleHelper.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <filesystem>
#include <iomanip>

// ---------------------------------------------------------------------------
// Main Simulation Function for Thelen 2003 Muscle Model (Passive State)
// ---------------------------------------------------------------------------
// This function performs a passive stretch simulation where the muscle length
// is externally controlled (prescribed motion), and the resulting force is calculated.
int run_TMM_passive() {
    std::cout << "Starting TMM Passive Simulation..." << std::endl;
    try {
        // --- 1. Simulation Configuration ---
        // Parameters match the MATLAB simulation for validation.

        double L_mtn_input = 1.05;    // Target Normalized Musculotendon Length multiplier
        double U_input     = 0.0;     // Excitation Input (u): Set to 0.0 for passive test
        double Amp_input   = 0.005;   // Amplitude of length oscillation (normalized)
        double SimTime_input = 100.0; // Total simulation time (s)
        double SimDt_input = 0.001;   // Time step (s)
        double FreqLow_input  = 0.1;  // Start Frequency (Hz)
        double FreqHigh_input = 100;  // End Frequency (Hz)
        int NumFreqSamples = 1000;    // Number of frequency steps

        // Frequency Sweep Settings
        double SimTime = SimTime_input;
        double dt = SimDt_input;
        double FreqLow = FreqLow_input;         // Start Frequency (Hz)
        double FreqHigh = FreqHigh_input;      // End Frequency (Hz)

       // Soleus Muscle Parameters
        double F_mo = 3549.0;           // Max Isometric Force (N)
        double L_mo = 0.05;             // Optimal Fiber Length (m)
        double L_ts = 0.25;             // Tendon Slack Length (m)
        double alphaOpt = 0.4363;       // Pennation Angle at Optimal Length (rad)
        double V_mmax_norm = 10.0;      // Max Contraction Velocity (L_opt/s)

        // --- 2. Muscle Object Setup ---
        // Create an instance of the Thelen2003MuscleWrapper class.
        Thelen2003MuscleWrapper m;
        m.set_max_isometric_force(F_mo);              // F_max (N)
        m.set_optimal_fiber_length(L_mo);             // L_opt (m)
        m.set_tendon_slack_length(L_ts);              // L_slack (m)
        m.set_pennation_angle_at_optimal(alphaOpt);   // Alpha_opt (rad)
        m.set_max_contraction_velocity(V_mmax_norm);  // V_max (L_opt/s)

        // Set Thelen 2003 Specific Model Parameters
        m.set_KshapeActive(0.45);        // Active shape factor
        m.set_KshapePassive(4.0);        // Passive shape factor (exponential)
        m.set_FmaxMuscleStrain(0.6);     // Muscle strain at max passive force
        m.set_FmaxTendonStrain(0.033);   // Tendon strain at max isometric force (e0)
        m.set_Af(0.3);                   // Force-velocity shape factor
        m.set_Flen(1.4);                 // Max normalized eccentric force
        m.set_activation_time_constant(0.01);
        m.set_deactivation_time_constant(0.04);

        // Initialize internal properties
        m.finalizeFromProperties();

        // --- Derived Constants ---
        // Muscle Width (Constant Volume assumption)
        double L_m_height = L_mo * std::sin(alphaOpt);

        // --- 3. Initialization (Finding Operating Point) ---
        // Calculate the reference and target musculotendon lengths.

        // Reference Length at Rest: L_mt0 = L_opt * cos(alpha_opt) + L_slack
        double L_mt0_ref = L_mo * std::cos(alphaOpt) + L_ts;

        // Target Mean Length for Simulation
        double L_mt_target = L_mtn_input * L_mt0_ref;

        // Bisection Search for Initial Fiber Length (L_mn)
        // We need to find L_mn such that the static forces balance at L_mt_target.
        // Equation: F_tendon - F_passive * cos(alpha) = 0 (Active force is 0)
        double low_L_mn = 0.0;
        double high_L_mn = 2.0;
        double L_mn_init = 0.0;
        double error1 = 1000.0;

        while (std::abs(error1) > 1e-11) {
            L_mn_init = 0.5 * (low_L_mn + high_L_mn);
            double L_m = L_mn_init * L_mo;

            // Calculate Pennation Angle
            double sin_phi = L_m_height / L_m;
            if (sin_phi > 1.0) sin_phi = 1.0;
            double alpha = std::asin(sin_phi);

            // Calculate Tendon Length
            double L_t = L_mt_target - L_m * std::cos(alpha);
            double L_tn = L_t / L_ts; // Normalized Tendon Length

            // Calculate Forces
            double F_tn = m.calcfse_(L_tn);      // Tendon Force
            double F_pn = m.calcfpe_(L_mn_init); // Passive Fiber Force

            // Error Calculation (Equilibrium check)
            error1 = F_tn - F_pn * std::cos(alpha);

            // Update Bisection Bounds
            if (error1 < 0) high_L_mn = L_mn_init;
            else            low_L_mn  = L_mn_init;
        }

        std::cout << "Init Done. L_mn0=" << L_mn_init << std::endl;

        // --- 4. Frequency Sweep Loop ---
        // Generate logarithmic frequency steps
        std::vector<double> frequencies = logspace(std::log10(FreqLow), std::log10(FreqHigh), NumFreqSamples);
        std::filesystem::create_directories("TMM_Passive_result_csv");

        // Save frequency list
        std::ofstream freqFile("TMM_Passive_result_csv/frequency_list.csv");
        if (freqFile.is_open()) {
            for (size_t k = 0; k < frequencies.size(); ++k) {
                freqFile << k << "," << frequencies[k] << "\n";
            }
            freqFile.close();
        }

        // Amplitude of length oscillation in meters
        double amp_val = Amp_input * L_mt0_ref * L_mtn_input;

        // Loop through each frequency
        for (int k = 0; k < NumFreqSamples; ++k) {
            double freq = frequencies[k];
            std::ofstream outFile("TMM_Passive_result_csv/freq_res_" + std::to_string(k) + ".csv");
            outFile << std::setprecision(16);
            outFile << "time,L_mt,F_m_AT\n"; // Header: Time, Input Length, Output Force

            // Reset state variables
            double L_mn = L_mn_init;
            double a = U_input; // Fixed at 0.0
            int time_len = (int)(SimTime / dt);

            // Time Integration Loop
            for (int i = 0; i < time_len; ++i) {
                double t = i * dt;

                // 1. Input: Prescribed Motion (L_mt oscillation)
                // L_mt(t) = Mean + Amplitude * sin(2*pi*f*t)
                double L_mt = L_mt_target + amp_val * std::sin(2.0 * M_PI * freq * t);

                // 2. Implicit Solver (Newton-Raphson)
                // Find L_mn that satisfies equilibrium for the current prescribed L_mt.
                // Note: The loop starts solving from the second step (i>0) to use velocity.

                if (i > 0) {
                    double L_mn_prev = L_mn;
                    double L_mnc = L_mn; // Initial guess
                    int max_iter = 50;
                    double error = 1.0;
                    double delta = 1e-7; // Step for numerical derivative
                    int iter = 0;
                    double tol = 1e-8;

                    while (std::abs(error) > tol && iter < max_iter) {
                        iter++;

                        // Error Calculation Lambda
                        auto calc_error = [&](double L_mn_val) -> double {
                            // Current Fiber Length
                            double L_m_val = L_mn_val * L_mo;

                            // Pennation Angle
                            double sin_phi = L_m_height / L_m_val;
                            if (sin_phi > 1.0) sin_phi = 1.0;
                            double cos_phi = std::cos(std::asin(sin_phi));

                            // Tendon Force Calculation
                            double L_t_val = L_mt - L_m_val * cos_phi;
                            double L_tn_val = L_t_val / L_ts;
                            double F_tn_val = m.calcfse_(L_tn_val);

                            // Fiber Velocity Calculation
                            // V_mn = (dL_mn / dt) / V_max
                            // Note: Uses backward difference (Current Guess - Previous State)
                            double V_mn_val = ((L_mn_val - L_mn_prev) * L_mo / dt) / (V_mmax_norm * L_mo);

                            // Fiber Force Components
                            double fl = m.calcfal_(L_mn_val);  // Active Force-Length Multiplier
                            double fp = m.calcfpe_(L_mn_val);  // Normalized Passive Force
                            double fv = calc_force_velocity(m, V_mn_val, a); // Active Force-Velocity Multiplier

                            // Total Fiber Force
                            double F_mn_val = a * fl * fv + fp;

                            // Equilibrium Error
                            return F_tn_val - F_mn_val * cos_phi;
                        };

                        double err1 = calc_error(L_mnc);
                        double err2 = calc_error(L_mnc + delta);

                        // Calculate Jacobian (Numerical Derivative)
                        double J = (err2 - err1) / delta;
                        if (std::abs(J) < 1e-14) J = 1e-14;

                        // Newton Update
                        L_mnc -= err1 / J;

                        // Apply Bounds
                        if (L_mnc < 0.01) L_mnc = 0.01;
                        if (L_mnc > 2.0) L_mnc = 2.0;

                        error = err1;
                    }
                    L_mn = L_mnc; // Update State
                }

                // 3. Output Calculation (Post-Convergence)
                // Re-calculate final forces based on the converged length
                double L_m_final = L_mn * L_mo;
                double sin_phi = L_m_height / L_m_final;
                if(sin_phi > 1.0) sin_phi = 1.0;
                double cos_phi = std::cos(std::asin(sin_phi));

                double L_t_final = L_mt - L_m_final * cos_phi;
                double L_tn_final = L_t_final / L_ts;
                double F_t_final = m.calcfse_(L_tn_final) * F_mo;

                // Log Data: Time, Input (Length), Output (Force)
                outFile << t << "," << L_mt << "," << F_t_final << "\n";
            }
            outFile.close();
            if (k % 10 == 0) std::cout << "Step " << k << " done." << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }
    return 0;
}