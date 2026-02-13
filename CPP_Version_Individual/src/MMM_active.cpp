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
int run_MMM_active() {
    std::cout << "Starting MMM Active Simulation..." << std::endl;
    try {
        // ------------------------------------------------------------------
        // 1. Configuration
        // ------------------------------------------------------------------

        double L_mn_input = 1.0;      // Target Normalized Muscle Length multiplier (at rest)
        double U_input = 0.5;         // Excitation Input (u): Constant activation level for active test
        double Amp_input = 0.01;      // Amplitude of input oscillation (if applicable)
        double Mass_input = 3e9;      // External Mass (kg): Very large mass simulates isometric conditions
        double Damping_input = 0.0;    // External Damping
        double SimTime_input = 120.0; // Total simulation time (s)
        double SimDt_input = 0.001;   // Time step (s)
        double FreqLow_input  = 0.1;  // Start Frequency (Hz)
        double FreqHigh_input = 100;  // End Frequency (Hz)
        double NumFreqSamples = 100;  // Number of frequency steps

        // External dynamics
        double Mass_ext = Mass_input;
        double Damping_ext = Damping_input;

        // Simulation Control Parameters
        double SimTime = SimTime_input;    // Total simulation time
        double dt = SimDt_input;           // Time step
        double FreqLow = FreqLow_input;    // Start Frequency (Hz)
        double FreqHigh = FreqHigh_input;  // End Frequency (Hz)

        // Soleus Muscle Parameters
        double F_mo = 3549.0;           // Max Isometric Force (N)
        double L_mo = 0.05;             // Optimal Fiber Length (m)
        double L_ts = 0.25;             // Tendon Slack Length (m)
        double alphaOpt = 0.4363;       // Pennation Angle at Optimal Length (rad)
        double V_mmax_norm = 10.0;      // Max Contraction Velocity (L_opt/s)

        // Initialize Millard2012 Muscle Model
        OpenSim::Millard2012EquilibriumMuscle m;
        m.setMaxIsometricForce(F_mo);
        m.setOptimalFiberLength(L_mo);
        m.setTendonSlackLength(L_ts);
        m.setPennationAngleAtOptimalFiberLength(alphaOpt);
        m.setMaxContractionVelocity(V_mmax_norm);

        // Set Damping to 0.0
        bool ignoreTendonCompliance = false;
        bool ignoreActivationDynamics = false;
        double dampingBeta = 0.0;
        m.setMuscleConfiguration(false, false, dampingBeta);

        // ------------------------------------------------------------------
        // [Curve 1] Tendon Force-Length Curve
        // ------------------------------------------------------------------
        // Defines the non-linear elastic properties of the tendon (Toe region + Linear region).
        // - strainAtOneNormForce: Tendon strain at normalized force 1.0 (F_mo). (e0)
        // - stiffnessAtOneNormForce: Stiffness (slope) at normalized force 1.0.
        // - normForceAtToeEnd: Normalized force limit where the nonlinear 'toe' region ends.
        // - curviness: Curvature of the toe region (0.0: linear, 1.0: highly curved).
        double t_strainAtOneNormForce    = 0.033;
        double t_stiffnessAtOneNormForce = 1.375 / 0.033;
        double t_normForceAtToeEnd       = 2.0 / 3.0;
        double t_curviness               = 0.5;

        m.setTendonForceLengthCurve(OpenSim::TendonForceLengthCurve(
            t_strainAtOneNormForce,
            t_stiffnessAtOneNormForce,
            t_normForceAtToeEnd,
            t_curviness
        ));

        // ------------------------------------------------------------------
        // [Curve 2] Active Force-Length Curve
        // ------------------------------------------------------------------
        // Defines the active force generation capability based on fiber length.
        // - minActiveNormLength: Normalized length where active force starts (Ascending limb start).
        // - transitionNormLength: Length where the steep ascending limb transitions to shallow.
        // - maxActiveNormLength: Normalized length where active force ends (Descending limb end).
        // - shallowAscendingSlope: Slope of the shallow ascending limb.
        // - minValue: Minimum value of the curve (0.0 is critical for correct physics).
        double a_minActiveNormLength    = 0.4441;
        double a_transitionNormLength   = 0.73;
        double a_maxActiveNormLength    = 1.8123;
        double a_shallowAscendingSlope  = 0.8616;
        double a_minValue               = 0.0;

        m.setActiveForceLengthCurve(OpenSim::ActiveForceLengthCurve(
            a_minActiveNormLength,
            a_transitionNormLength,
            a_maxActiveNormLength,
            a_shallowAscendingSlope,
            a_minValue
        ));

        // ------------------------------------------------------------------
        // [Curve 3] Force-Velocity Curve
        // ------------------------------------------------------------------
        // Defines the force modulation based on contraction velocity (Hill's Curve).
        // - concSlopeAtVmax: Slope at maximum shortening velocity (usually 0).
        // - concSlopeNearV0: Slope near zero velocity during shortening.
        // - isoSlope: Slope at isometric state (V=0). High value (e.g., 5.0) adds numerical stability.
        // - eccSlopeAtVmax: Slope at maximum lengthening velocity.
        // - eccSlopeNearV0: Slope near zero velocity during lengthening.
        // - maxEccentricForce: Max force multiplier during lengthening (F_len).
        // - concCurviness: Curvature of the shortening region.
        // - eccCurviness: Curvature of the lengthening region.
        double fv_concSlopeAtVmax      = 0.0;
        double fv_concSlopeNearV0      = 0.15;
        double fv_isoSlope             = 5.0;
        double fv_eccSlopeAtVmax       = 0.1;
        double fv_eccSlopeNearV0       = 0.15;
        double fv_maxEccentricForce    = 1.4;
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

        // Get references to curve objects for direct calculation
        const auto& fvCurve = m.getForceVelocityCurve();
        const auto& flCurve = m.getActiveForceLengthCurve();
        const auto& fpCurve = m.getFiberForceLengthCurve();
        const auto& ftCurve = m.getTendonForceLengthCurve();

        // ------------------------------------------------------------------
        // 2. Initialization (Finding Equilibrium State)
        // ------------------------------------------------------------------

        // Assume initial fiber length is at optimal length
        double L_m_init = L_mo * 1.0;

        // Calculate constant muscle width (Fixed Width Pennation Model)
        double L_m_height = L_mo * std::sin(alphaOpt);

        // Calculate initial fiber force components
        double fl_init = flCurve.calcValue(1.0);
        double fp_init = fpCurve.calcValue(1.0);
        // Initial fiber force (Assuming velocity is zero)
        double F_m_init = (U_input * fl_init + fp_init) * F_mo;

        // Calculate initial pennation angle
        double phi_init = std::asin(L_m_height / L_m_init);

        // Calculate required tendon force for equilibrium: F_tendon = F_fiber * cos(phi)
        double F_t_prev = F_m_init * std::cos(phi_init);

        // Solve for initial tendon length given the force
        double L_t_init = calc_tendon_length_from_force(m, F_t_prev);

        // Calculate initial total musculotendon length
        double L_mt0 = L_m_init * std::cos(phi_init) + L_t_init;
        double F_ext_equil = F_t_prev; // Equilibrium force for external mass

        std::cout << "Initialization Done. L_t_init = " << L_t_init << std::endl;

        // ------------------------------------------------------------------
        // 3. Frequency Sweep Loop
        // ------------------------------------------------------------------
        // Generate frequencies from 0.1 Hz to 100 Hz
        std::vector<double> frequencies = logspace(std::log10(FreqLow), std::log10(FreqHigh), NumFreqSamples);
        std::filesystem::create_directories("MMM_Active_results_csv");

        // [ADD] Generate frequency_list.csv for MATLAB
        std::ofstream freqFile("MMM_Active_results_csv/frequency_list.csv");
        if (freqFile.is_open()) {
            for (size_t k = 0; k < frequencies.size(); ++k) {
                freqFile << k << "," << frequencies[k] << "\n";
            }
            freqFile.close();
        }

        for (int k = 0; k < 100; ++k) {
            double freq = frequencies[k];
            std::ofstream outFile("MMM_Active_results_csv/freq_res_" + std::to_string(k) + ".csv");
            outFile << "time,u,F_m_AT\n";

            double L_mt = L_mt0;
            double V_mt = 0.0;
            double a = U_input;
            double L_m = L_m_init; // Reset fiber length for each frequency
            int time_len = (int)(SimTime / dt);

            // Time Integration Loop
            for (int i = 0; i < time_len; ++i) {
                double t = i * dt;

                // 1. Input: Sinusoidal Excitation (u)
                double u = sinwave(freq, t, U_input, U_input, Amp_input);
                if (u < 0) u = 0; if (u > 1) u = 1;

                // 2. Activation Dynamics (First-order ODE)
                // da/dt = (u - a) / tau
                if (i > 0) a += dt * ((u - a) / get_Tau_MMM(a, u));

                // 3. Muscle Dynamics (Fixed Width Pennation Model)

                double L_m_prev = L_m;
                double V_mt_prev = V_mt;

                // Solver variable setting
                double L_m_curr = L_m;  // Initial guess
                int max_iter = 50;
                double tol = 1e-8;
                double delta = 1e-7;
                double error = 1.0;
                int iter = 0;

                // Newton-Raphson Loop
                while (std::abs(error) > tol && iter < max_iter) {
                    iter++;

                    auto calc_force_error = [&](double l_val) -> double {
                        // A. Geometry Update
                        double sin_phi = L_m_height / l_val;
                        if (sin_phi > 1.0) sin_phi = 1.0; 
                        double cos_phi = std::sqrt(1.0 - sin_phi * sin_phi);

                        // B. Tendon Force
                        double L_t = L_mt - l_val * cos_phi;
                        double L_t_norm = L_t / L_ts;
                        double F_tn = ftCurve.calcValue(L_t_norm);
                        double F_t = F_tn * F_mo;

                        // C. Fiber Kinematics (Implicit: V = dL/dt)
                        double V_m_val = (l_val - L_m_prev) / dt;
                        double V_mn = V_m_val / (V_mmax_norm * L_mo);

                        // D. Fiber Forces
                        double L_m_norm_val = l_val / L_mo;
                        double fl = flCurve.calcValue(L_m_norm_val);
                        double fp = fpCurve.calcValue(L_m_norm_val);
                        double fv = fvCurve.calcValue(V_mn); // 계산된 속도 반영

                        // E. Equilibrium Error
                        // F_fiber = a * fl * fv + fp + damping
                        double F_fiber_norm = a * fl * fv + fp + dampingBeta * V_mn;
                        double F_fiber_proj = (F_fiber_norm * F_mo) * cos_phi;

                        return F_t - F_fiber_proj;
                    };

                    // Newton Step
                    double err1 = calc_force_error(L_m_curr);
                    double err2 = calc_force_error(L_m_curr + delta); // Perturbation

                    double J = (err2 - err1) / delta; // Jacobian
                    if (std::abs(J) < 1e-14) J = 1e-14;

                    L_m_curr -= err1 / J; // Update L_m

                    // Bounds Check
                    if (L_m_curr < 0.001) L_m_curr = 0.001; 
                    if (L_m_curr > 0.5)   L_m_curr = 0.5;   

                    error = err1;
                }
                
                L_m = L_m_curr;

                // Geometry
                double sin_phi = L_m_height / L_m;
                if(sin_phi > 1.0) sin_phi = 1.0;
                double cos_phi = std::sqrt(1.0 - sin_phi * sin_phi);
                
                // Tendon Force
                double L_t = L_mt - L_m * cos_phi;
                double L_t_norm = L_t / L_ts;
                double F_t_final = ftCurve.calcValue(L_t_norm) * F_mo;

                // External Dynamics (Single Mass)
                // F = ma -> a = F/m -> v = v + a*dt -> x = x + v*dt
                double A_mt = (F_ext_equil - F_t_final - V_mt_prev * Damping_ext) / Mass_ext;
                V_mt += A_mt * dt;
                L_mt += V_mt * dt;

                // Output Data
                outFile << t << "," << u << "," << F_t_final << "\n";
            }
            outFile.close();
            if (k % 10 == 0) std::cout << "Frequency step " << k << " done." << std::endl;
        }
    } catch (const std::exception& e) { std::cerr << e.what() << std::endl; }
    return 0;
}