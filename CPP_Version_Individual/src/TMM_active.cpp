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
// Main Simulation Function for Thelen 2003 Muscle Model (Active State)
// ---------------------------------------------------------------------------
// This function performs a forward dynamics simulation of a muscle-tendon unit
// attached to a mass, specifically testing the active force generation.
int run_TMM_active(double L_mn_input) {
    std::cout << "Starting TMM Active Simulation..." << std::endl;
    try {
        // --- 1. Simulation Configuration ---
        // Basic parameters defining the simulation environment and muscle inputs.

        // double L_mn_input = 1.0;   // Target Normalized Muscle Length multiplier (at rest)
        double U_input = 0.5;         // Excitation Input (u): Constant activation level for active test
        double Amp_input = 0.01;      // Amplitude of input oscillation (if applicable)
        double Mass_input = 3e9;      // External Mass (kg): Very large mass simulates isometric conditions
        double Damping_input = 0.0;   // External Damping
        double SimTime_input = 100.0; // Total simulation time (s)
        double SimDt_input = 0.001;   // Time step (s)
        double FreqLow_input  = 0.1;  // Start Frequency (Hz)
        double FreqHigh_input = 100;  // End Frequency (Hz)
        double NumFreqSamples = 1000; // Number of frequency steps

        // External dynamics
        double Mass_ext = Mass_input;
        double Damping_ext = Damping_input;

        // Frequency Sweep Settings for Bode Analysis
        double SimTime = SimTime_input;
        double dt = SimDt_input;
        double FreqLow = FreqLow_input;         // Start Frequency (Hz)
        double FreqHigh = FreqHigh_input;      // End Frequency (Hz)

        // --- 2. Muscle Object Setup ---
        // Create an instance of the Thelen2003MuscleWrapper class.
        // This wrapper exposes protected OpenSim functions for direct calculation.
        Thelen2003MuscleWrapper m;

       // Soleus Muscle Parameters
        double F_mo = 3549.0;           // Max Isometric Force (N)
        double L_mo = 0.05;             // Optimal Fiber Length (m)
        double L_ts = 0.25;             // Tendon Slack Length (m)
        double alphaOpt = 0.4363;       // Pennation Angle at Optimal Length (rad)
        double V_mmax_norm = 10.0;      // Max Contraction Velocity (L_opt/s)

        // Set Standard Muscle Properties (Based on Soleus muscle data)
        m.set_max_isometric_force(F_mo);              // F_max (N)
        m.set_optimal_fiber_length(L_mo);             // L_opt (m)
        m.set_tendon_slack_length(L_ts);              // L_slack (m)
        m.set_pennation_angle_at_optimal(alphaOpt);   // Alpha_opt (rad)
        m.set_max_contraction_velocity(V_mmax_norm);  // V_max (L_opt/s)

        // Set Thelen 2003 Specific Model Parameters
        // These parameters define the shape of force-length and force-velocity curves.
        m.set_KshapeActive(0.45);        // Shape factor for active force-length curve
        m.set_KshapePassive(4.0);        // Shape factor for passive force-length curve (exponential)
        m.set_FmaxMuscleStrain(0.6);     // Muscle strain at max passive force
        m.set_FmaxTendonStrain(0.033);   // Tendon strain at max isometric force (e0)
        m.set_Af(0.3);                   // Force-velocity shape factor (concentric)
        m.set_Flen(1.4);                 // Max normalized eccentric force
        m.set_activation_time_constant(0.01);   // Activation time constant (tau_act)
        m.set_deactivation_time_constant(0.04); // Deactivation time constant (tau_deact)

        // Initialize internal properties (OpenSim boilerplate)
        m.finalizeFromProperties();

        // --- 3. Derived Constants ---
        // Muscle Width (Constant Volume assumption)
        // h = L_opt * sin(alpha_opt)
        double L_m_height = L_mo * std::sin(alphaOpt);

        // --- 4. Initialization (Finding Equilibrium State) ---
        // We need to find the initial state (lengths and forces) that satisfies equilibrium.

        // Initial normalized fiber length guess
        double L_mn0 = L_mn_input;
        double L_m0 = L_mn0 * L_mo;

        // Calculate initial forces assuming zero velocity
        double f_l0 = m.calcfal_(L_mn0);    // Active force multiplier
        double F_pn0 = m.calcfpe_(L_mn0);   // Passive force multiplier
        double f_v0 = 1.0;                  // Velocity multiplier (at V=0)

        // Initial pennation angle
        double alpha0 = std::asin(L_m_height / L_m0);

        // Calculate required tendon force for equilibrium: F_tendon = F_fiber * cos(alpha)
        double F_tn_init = (U_input * f_l0 * f_v0 + F_pn0) * std::cos(alpha0);

        // Inverse calculation: Find tendon strain given the force
        double eps_t_init = calc_tendon_strain_from_force(m, F_tn_init);

        // Calculate initial tendon length: L_t = L_slack * (1 + strain)
        double L_t0 = L_ts * (1.0 + eps_t_init);

        // Initial total musculotendon length: L_mt = L_m * cos(alpha) + L_t
        double L_mt0 = L_m0 * std::cos(alpha0) + L_t0;

        // Equilibrium external force (Force exerted by the muscle on the mass)
        double F_ext_equil = F_mo * F_tn_init;

        std::cout << "Init: L_mt0=" << L_mt0 << ", F_equil=" << F_ext_equil << std::endl;

        // --- 5. Frequency Sweep Loop ---
        // Generate logarithmic frequency steps
        std::vector<double> frequencies = logspace(std::log10(FreqLow), std::log10(FreqHigh), NumFreqSamples);

        std::string dir_name = "TMM_Active_result_csv_" + std::to_string(L_mn_input);
        std::filesystem::create_directories(dir_name);

        // Save frequency list for post-processing
        std::ofstream freqFile(dir_name + "/frequency_list.csv");
        if (freqFile.is_open()) {
            for (size_t k = 0; k < frequencies.size(); ++k) {
                freqFile << k << "," << frequencies[k] << "\n";
            }
            freqFile.close();
        }

        // Loop through each frequency
        for (int k = 0; k < NumFreqSamples; ++k) {
            double freq = frequencies[k];
            std::ofstream outFile(dir_name + "/freq_res_" + std::to_string(k) + ".csv");
            outFile << std::setprecision(16);
            outFile << "time, u, F_m_AT\n";

            // Reset state variables for each frequency run
            double L_mt = L_mt0;    // Start from equilibrium length
            double V_mt = 0.0;      // Start from static condition
            double L_mn = L_mn0;    // Start from normalized fiber length 1.0
            double a = U_input;     // Activation
            int time_len = (int)(SimTime / dt);

            // Time Integration Loop
            for (int i = 0; i < time_len; ++i) {
                double t = i * dt;

                // 1. Input: Sinusoidal Excitation (u)
                // u = u0 + amp * sin(2*pi*f*t)
                double u = sinwave(freq, t, U_input, U_input, Amp_input);
                if (u < 0) u = 0; if (u > 1) u = 1; // Clamp between 0 and 1

                // 2. Activation Dynamics
                // First-order differential equation: da/dt = (u - a) / tau
                if (i > 0) {
                    double tau = calc_tau(m, a, u); // Calculate time constant
                    double da_dt = (u - a) / tau;
                    a += dt * da_dt;
                }

                // 3. Implicit Solver (Newton-Raphson)
                // Find normalized fiber length (L_mnc) that satisfies force equilibrium.
                // Equation: F_tendon - F_fiber_projected = 0
                double L_mn_prev = L_mn;
                double V_mt_prev = V_mt;
                double L_mnc = L_mn; // Initial guess = previous step
                int max_iter = 50;
                double error = 1.0;
                double delta = 1e-7; // Step for numerical derivative
                int iter = 0;

                while (std::abs(error) > 1e-8 && iter < max_iter) {
                    iter++;

                    // Lambda function to calculate error at a given length (L_mn_val)
                    auto calc_error = [&](double L_mn_val) -> double {
                        // Current Fiber Length
                        double L_m_val = L_mn_val * L_mo;

                        // Pennation Angle Geometry
                        double sin_phi = L_m_height / L_m_val;
                        if (sin_phi > 1.0) sin_phi = 1.0;
                        double cos_phi = std::cos(std::asin(sin_phi));

                        // Calculate Tendon Length: L_t = L_mt - L_m * cos(phi)
                        double L_t_val = L_mt - L_m_val * cos_phi;

                        // Normalized Tendon Length (L_tn = 1.0 + strain)
                        // Note: Thelen's calcfse takes normalized length, not strain directly.
                        double L_tn_val = L_t_val / L_ts;

                        // Tendon Force
                        double F_tn_val = m.calcfse_(L_tn_val);

                        // Fiber Velocity
                        // V_mn = (dL_mn / dt) / V_max
                        double V_mn_val = ((L_mn_val - L_mn_prev) * L_mo / dt) / (V_mmax_norm * L_mo);

                        // Fiber Force Components
                        double fl = m.calcfal_(L_mn_val); // active force-length multiplier
                        double fp = m.calcfpe_(L_mn_val); // passive force normalized

                        // Calculate active Force-Velocity multiplier (Inverse logic compared to OpenSim)
                        double fv = calc_force_velocity(m, V_mn_val, a);

                        // Total Fiber Force (Normalized)
                        double F_mn_val = a * fl * fv + fp;

                        // Error: Tendon Force - Projected Fiber Force
                        return F_tn_val - F_mn_val * cos_phi;
                    };

                    // Newton-Raphson Step
                    double err1 = calc_error(L_mnc);
                    double err2 = calc_error(L_mnc + delta);

                    // Jacobian (Derivative)
                    double J = (err2 - err1) / delta;
                    if (std::abs(J) < 1e-14) J = 1e-14; // Prevent division by zero

                    // Update Guess
                    L_mnc -= err1 / J;

                    // Apply Bounds
                    if (L_mnc < 0.01) L_mnc = 0.01;
                    if (L_mnc > 2.0) L_mnc = 2.0;

                    error = err1;
                }
                L_mn = L_mnc; // Update state

                // 4. Output Force Calculation (Post-Convergence)
                double L_m_final = L_mn * L_mo;
                double sin_phi = L_m_height / L_m_final;
                if (sin_phi > 1.0) sin_phi = 1.0;
                double cos_phi = std::cos(std::asin(sin_phi));
                double L_t_final = L_mt - L_m_final * cos_phi;

                // Final Tendon Force (in Newtons)
                double L_tn_final = L_t_final / L_ts;
                double F_t_final = m.calcfse_(L_tn_final) * F_mo;

                // 5. External Dynamics Integration (Mass)
                // F = ma -> a = F_net / m
                double A_mt = (F_ext_equil - F_t_final - V_mt_prev * Damping_ext) / Mass_ext;
                V_mt += A_mt * dt;
                L_mt += V_mt * dt;

                // Log Data
                outFile << t << "," << u << "," << F_t_final << "\n";
            }
            outFile.close();
            if (k % 10 == 0) std::cout << "Step " << k << " done." << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }
    return 0;
}