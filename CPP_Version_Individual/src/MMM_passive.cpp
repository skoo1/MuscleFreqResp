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

int run_MMM_passive() {
    std::cout << "Starting MMM Passive Simulation..." << std::endl;
    try {
        // ------------------------------------------------------------------
        // 1. Configuration (Parameters matching the MATLAB simulation)
        // ------------------------------------------------------------------
        double L_mtn_input = 1.05;      // Target Normalized Musculotendon Length multiplier
        double U_input = 0.0;   // Input Excitation (Very small value effectively 0 for passive)
        double Amp_input = 0.005;       // Amplitude of the oscillation (normalized)
        double SimTime_input = 120.0;   // Total simulation time (seconds)
        double SimDt_input = 0.001;     // Simulation time step (seconds)
        double FreqLow_input = 0.1;     // Start frequency for sweep (Hz)
        double FreqHigh_input = 100.0;  // End frequency for sweep (Hz)
        int NumFreqSamples = 100;       // Number of frequency steps

        // Frequency Sweep Settings
        double SimTime = SimTime_input;
        double dt = SimDt_input;
        double FreqLow = FreqLow_input;         // Start Frequency (Hz)
        double FreqHigh = FreqHigh_input;      // End Frequency (Hz)

        // Muscle Properties (Soleus muscle parameters)
        double F_mo = 3549.0;           // Max Isometric Force (N)
        double L_mo = 0.05;             // Optimal Fiber Length (m)
        double L_ts = 0.25;             // Tendon Slack Length (m)
        double alphaOpt = 0.4363;       // Pennation Angle at Optimal Length (rad)
        double V_mmax_norm = 10.0;      // Max Contraction Velocity (L_opt/s)
        double dampingBeta = 0.0;
        // ------------------------------------------------------------------
        // 2. OpenSim Muscle Object Setup
        // ------------------------------------------------------------------
        // Initialize the Millard2012EquilibriumMuscle object
        OpenSim::Millard2012EquilibriumMuscle m;
        m.setMaxIsometricForce(F_mo);
        m.setOptimalFiberLength(L_mo);
        m.setTendonSlackLength(L_ts);
        m.setPennationAngleAtOptimalFiberLength(alphaOpt);
        m.setMaxContractionVelocity(V_mmax_norm);

        // -------------------------------------------------------------------------
        // [Curve 1] Tendon Force-Length Curve
        // -------------------------------------------------------------------------
        // Defines the nonlinear elastic properties of the tendon.
        // - strainAtOneNormForce: Tendon strain at normalized force 1.0 (F_mo). (Typically 0.033 ~ 0.04)
        // - stiffnessAtOneNormForce: Slope (stiffness) of the curve at normalized force 1.0.
        // - normForceAtToeEnd: Normalized force limit where the nonlinear 'toe' region ends and the linear region begins (0.0 ~ 1.0).
        // - curviness: Curvature of the toe region (0.0: linear, 1.0: highly curved).
        double t_strainAtOneNormForce    = 0.033;             // e0
        double t_stiffnessAtOneNormForce = 1.375 / 0.033;     // k_iso (Linear stiffness)
        double t_normForceAtToeEnd       = 2.0 / 3.0;         // F_toe (Toe limit)
        double t_curviness               = 0.5;               // Curve shape

        m.setTendonForceLengthCurve(OpenSim::TendonForceLengthCurve(
            t_strainAtOneNormForce,
            t_stiffnessAtOneNormForce,
            t_normForceAtToeEnd,
            t_curviness
        ));

        // -------------------------------------------------------------------------
        // [Curve 2] Active Force-Length Curve
        // -------------------------------------------------------------------------
        // Defines the active force generation capability of the fiber based on its length.
        // - minActiveNormLength: Minimum normalized fiber length where active force starts (Start of ascending limb).
        // - transitionNormLength: Transition point from steep ascending limb to shallow ascending limb.
        // - maxActiveNormLength: Maximum normalized fiber length where active force becomes zero (End of descending limb).
        // - shallowAscendingSlope: Slope of the shallow ascending limb.
        // - minValue: Minimum value of the curve (Set to 0.0 or a small value for numerical stability).
        double a_minActiveNormLength    = 0.4441;
        double a_transitionNormLength   = 0.73;
        double a_maxActiveNormLength    = 1.8123;
        double a_shallowAscendingSlope  = 0.8616;
        double a_minValue               = 0.0;

        // [Curve 3] Passive Force-Length Curve
        // -------------------------------------------------------------------------
        // Defines the passive elastic force (parallel elastic element) of the fiber.
        // - strainAtZeroForce: Fiber strain at which passive force starts to develop (Usually 0.0).
        // - strainAtOneNormForce: Fiber strain at normalized force 1.0 (F_mo).
        // - stiffnessAtLowForce: Stiffness (slope) at the beginning of force generation (low force).
        // - stiffnessAtOneNormForce: Stiffness (slope) at normalized force 1.0.
        // - curviness: Curvature of the passive curve.
        double p_strainAtZeroForce       = 0.0;
        double p_strainAtOneNormForce    = 0.7;          // e0_m (Fiber strain at Fmax)
        double p_stiffnessAtLowForce     = 0.2;          // k_low
        double p_stiffnessAtOneNormForce = 2.0 / 0.7;    // k_iso (Stiffness at Fmax)
        double p_curviness               = 0.75;

        m.setFiberForceLengthCurve(OpenSim::FiberForceLengthCurve(
            p_strainAtZeroForce,
            p_strainAtOneNormForce,
            p_stiffnessAtLowForce,
            p_stiffnessAtOneNormForce,
            p_curviness
        ));

        // -------------------------------------------------------------------------
        // [Curve 4] Force-Velocity Curve
        // -------------------------------------------------------------------------
        // Defines the force-velocity relationship (Hill's curve). (Total 8 parameters)
        // [Concentric: Shortening]
        // - concSlopeAtVmax: Slope at maximum shortening velocity (Vmax). (Usually 0).
        // - concSlopeNearV0: Slope of the concentric region near zero velocity.
        // - isoSlope: Slope at isometric state (V=0). (Used > 1.0 for damping-like stability).
        // [Eccentric: Lengthening]
        // - eccSlopeAtVmax: Slope at maximum lengthening velocity.
        // - eccSlopeNearV0: Slope of the eccentric region near zero velocity.
        // - maxEccentricForce: Maximum normalized force during eccentric contraction (F_len, typically 1.4 ~ 1.8).
        // [Curvature]
        // - concCurviness: Curvature of the concentric region.
        // - eccCurviness: Curvature of the eccentric region.
        double fv_concSlopeAtVmax      = 0.0;
        double fv_concSlopeNearV0      = 0.15;
        double fv_isoSlope             = 5.0;     // Damping like behavior at V=0
        double fv_eccSlopeAtVmax       = 0.1;
        double fv_eccSlopeNearV0       = 0.15;
        double fv_maxEccentricForce    = 1.4;     // F_len
        double fv_concCurviness        = 0.7;
        double fv_eccCurviness         = 0.9;

        m.setForceVelocityCurve(OpenSim::ForceVelocityCurve(
            fv_concSlopeAtVmax,
            fv_concSlopeNearV0,
            fv_isoSlope,
            fv_eccSlopeAtVmax,
            fv_eccSlopeNearV0,
            fv_maxEccentricForce,
            fv_concCurviness,
            fv_eccCurviness
        ));

        // Get references to the curves for direct calculation in the loop
        const auto& fvCurve = m.getForceVelocityCurve();
        const auto& flCurve = m.getActiveForceLengthCurve();
        const auto& fpCurve = m.getFiberForceLengthCurve();
        const auto& ftCurve = m.getTendonForceLengthCurve();

        // ------------------------------------------------------------------
        // 3. Initialization
        // ------------------------------------------------------------------
        // Fixed width
        double L_m_height = L_mo * std::sin(alphaOpt);

        // Calculate reference length at rest (L_mt0 = L_ts + L_mo * cos(alpha))
        double L_mt0_ref = L_ts + L_mo * std::cos(alphaOpt);

        // Target Musculotendon Length
        double L_mt_target = L_mt0_ref * L_mtn_input;

        // Bisection method to find the initial Fiber Length (projected along tendon)
        // that satisfies equilibrium (F_tendon = F_fiber_projected) at the target L_mt.
        double low_L_m_AT = 0.0;
        double high_L_m_AT = L_mt_target - 1e-5;
        double L_m_AT_init = 0.0;

        // Simple Bisection
        for(int iter=0; iter<100; ++iter) {
            L_m_AT_init = 0.5 * (low_L_m_AT + high_L_m_AT);

            // Calc Forces (Static)
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
        std::cout << "Init Done. L_m_AT_init = " << L_m_AT_init << std::endl;
        // ------------------------------------------------------------------
        // 4. Frequency Sweep Loop
        // ------------------------------------------------------------------
        // Create logarithmic frequency steps
        std::vector<double> frequencies = logspace(std::log10(FreqLow_input), std::log10(FreqHigh_input), NumFreqSamples);

        // Create directory for output files
        std::filesystem::create_directories("MMM_Passive_results_csv");

        // Generate frequency_list.csv for MATLAB analysis
        // This file maps the index 'k' to the actual frequency value.
        std::ofstream freqFile("MMM_Passive_results_csv/frequency_list.csv");
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

            // Open output CSV file for this frequency
            std::ofstream outFile("MMM_Passive_results_csv/freq_res_" + std::to_string(k) + ".csv");
            outFile << std::setprecision(16);
            outFile << "time,L_mt,F_m_AT\n"; // Header: Time, Input (Length), Output (Force)

            double L_m_AT = L_m_AT_init; // Initialize state with equilibrium value
            double a = U_input;
            double amp_val = Amp_input * L_mt0_ref * L_mtn_input;
            int time_len = (int)(SimTime / dt);

            // Time Loop for Simulation
            for (int i = 0; i < time_len; ++i) {
                double t = i * SimDt_input;

                // 1. Calculate Input: L_mt oscillation (Sine wave)
                double L_mt = L_mt_target + amp_val * std::sin(2.0 * M_PI * freq * t);

                double L_m_prev_sq = L_m_AT * L_m_AT + L_m_height * L_m_height;
                double L_m_prev_val = std::sqrt(L_m_prev_sq);

                if (i > 0) {
                    double L_m_AT_curr = L_m_AT; // Initial guess
                    int max_iter = 50;
                    double tol = 1e-8;
                    double delta = 1e-7;
                    double error = 1.0;
                    int iter = 0;

                    while (std::abs(error) > tol && iter < max_iter) {
                        iter++;

                        auto calc_force_error = [&](double l_at_val) -> double {
                            // A. Geometry
                            // L_m^2 = L_m_AT^2 + w^2
                            double L_m_curr_val = std::sqrt(l_at_val * l_at_val + L_m_height * L_m_height);
                            double cos_phi = l_at_val / L_m_curr_val;

                            // B. Tendon Force
                            double L_t = L_mt - l_at_val;
                            double f_t_norm = ftCurve.calcValue(L_t / L_ts);
                            double F_t = f_t_norm * F_mo;

                            // C. Fiber Kinematics (Implicit: V = dL_fiber/dt)
                            double V_m = (L_m_curr_val - L_m_prev_val) / SimDt_input;
                            double V_m_norm = V_m / (V_mmax_norm * L_mo);

                            // D. Fiber Forces
                            double L_m_norm = L_m_curr_val / L_mo;
                            double f_l = flCurve.calcValue(L_m_norm);
                            double f_p = fpCurve.calcValue(L_m_norm);
                            double f_v = fvCurve.calcValue(V_m_norm);

                            // F_fiber = a * fl * fv + fp + damping
                            double F_fib_norm = a * f_l * f_v + f_p + dampingBeta * V_m_norm;
                            double F_fib_proj = (F_fib_norm * F_mo) * cos_phi;

                            return F_t - F_fib_proj;
                        };

                        // Newton Step
                        double err1 = calc_force_error(L_m_AT_curr);
                        double err2 = calc_force_error(L_m_AT_curr + delta);

                        double J = (err2 - err1) / delta;
                        if (std::abs(J) < 1e-14) J = 1e-14;

                        L_m_AT_curr -= err1 / J;

                        // Bounds
                        if (L_m_AT_curr < 1e-6) L_m_AT_curr = 1e-6;
                        if (L_m_AT_curr > L_mt) L_m_AT_curr = L_mt - 1e-6;

                        error = err1;
                    }
                    L_m_AT = L_m_AT_curr; // Update State
                }

                // 3. Output Calculation
                double L_t_final = L_mt - L_m_AT;
                double F_t_final = ftCurve.calcValue(L_t_final / L_ts) * F_mo;

                outFile << t << "," << L_mt << "," << F_t_final << "\n";
            }
            outFile.close();
            if (k % 10 == 0) std::cout << "Frequency step " << k << " done." << std::endl;
        }
    } catch (const std::exception& e) { std::cerr << e.what() << std::endl; }
    return 0;
}