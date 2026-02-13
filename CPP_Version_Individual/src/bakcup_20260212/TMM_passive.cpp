#include "Thelen2003MuscleWrapper.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <filesystem>
#include <iomanip>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

namespace fs = std::filesystem;

// Utils
std::vector<double> logspace(double startExp, double endExp, int num) {
    std::vector<double> result;
    if (num <= 1) { result.push_back(std::pow(10, endExp)); return result; }
    double step = (endExp - startExp) / (num - 1);
    for (int i = 0; i < num; ++i) result.push_back(std::pow(10, startExp + i * step));
    return result;
}

// ---------------------------------------------------------------------------
// Main Simulation
// ---------------------------------------------------------------------------
int main() {
    try {
        // --- 1. Configuration (Matching MATLAB Code) ---
        double L_mtn_input = 1.05; // Target L_mt multiplier
        double U_input     = 0.0;  // Passive Test
        double Amp_input   = 0.005;
        double SimTime     = 120.0;
        double dt          = 0.001;
        double FreqLow     = 0.1;
        double FreqHigh    = 100.0;
        int NumFreqSamples = 100;

        // --- 2. Muscle Object Setup ---
        Thelen2003MuscleWrapper m;
        m.set_max_isometric_force(3549.0);
        m.set_optimal_fiber_length(0.05);
        m.set_tendon_slack_length(0.25);
        m.set_pennation_angle_at_optimal(0.4363);
        m.set_max_contraction_velocity(10.0);

        // Thelen Parameters
        m.set_KshapeActive(0.45);
        m.set_KshapePassive(4.0);
        m.set_FmaxMuscleStrain(0.6);
        m.set_FmaxTendonStrain(0.033);
        m.set_Af(0.3);
        m.set_Flen(1.4);
        m.set_activation_time_constant(0.01);
        m.set_deactivation_time_constant(0.04);
        m.finalizeFromProperties(); 

        // Derived Constants
        double L_mo = m.get_optimal_fiber_length();
        double L_ts = m.get_tendon_slack_length();
        double alphaOpt = m.get_pennation_angle_at_optimal();
        double F_mo = m.get_max_isometric_force();
        double V_mmax = m.get_max_contraction_velocity() * L_mo; 
        
        // Fixed Height for Pennation
        double L_m_height = L_mo * std::sin(alphaOpt); 

        // --- 3. Initialization (Finding Operating Point) ---
        // L_mt0 = L_mo * cos(alphaOpt) + L_ts
        double L_mt0_ref = L_mo * std::cos(alphaOpt) + L_ts;
        double L_mt_target_val = L_mtn_input * L_mt0_ref;

        // Bisection Search for L_mn (Static Equilibrium at L_mt_target_val)
        double low_L_mn = 0.0;
        double high_L_mn = 2.0;
        double L_mn_init = 0.0;
        double error1 = 1000.0;

        while (std::abs(error1) > 1e-11) {
            L_mn_init = 0.5 * (low_L_mn + high_L_mn);
            double L_m = L_mn_init * L_mo;
            
            // Pennation Angle
            double sin_phi = L_m_height / L_m;
            if (sin_phi > 1.0) sin_phi = 1.0;
            double alpha = std::asin(sin_phi);
            
            // Tendon Length
            double L_t = L_mt_target_val - L_m * std::cos(alpha);
            double tlN = L_t / L_ts; // Normalized Tendon Length (1+eps)
            
            // Forces
            double F_tn = m.calcfse_(tlN);
            double F_pn = m.calcfpe_(L_mn_init);
            
            // Error = Tendon - Passive * cos(alpha) (Since U=0, Active=0)
            error1 = F_tn - F_pn * std::cos(alpha);
            
            if (error1 < 0) high_L_mn = L_mn_init;
            else            low_L_mn  = L_mn_init;
        }

        std::cout << "Init Done. L_mn0=" << L_mn_init << std::endl;

        // --- 4. Frequency Sweep ---
        std::vector<double> frequencies = logspace(std::log10(FreqLow), std::log10(FreqHigh), NumFreqSamples);
        fs::create_directories("TMM_Passive_results_csv");

        // Generate frequency_list.csv
        std::ofstream freqFile("TMM_Passive_results_csv/frequency_list.csv");
        if (freqFile.is_open()) {
            for (size_t k = 0; k < frequencies.size(); ++k) {
                freqFile << k << "," << frequencies[k] << "\n";
            }
            freqFile.close();
        }

        double amp_val = Amp_input * L_mt0_ref * L_mtn_input;

        for (int k = 0; k < NumFreqSamples; ++k) {
            double freq = frequencies[k];
            std::ofstream outFile("TMM_Passive_results_csv/freq_res_" + std::to_string(k) + ".csv");
            outFile << std::setprecision(16);
            outFile << "time,L_mt,F_m_AT\n"; // Header changed: L_mt is input

            double L_mn = L_mn_init; 
            double a = U_input; // Fixed at 0.0
            int time_len = (int)(SimTime / dt);

            // Time Loop
            for (int i = 0; i < time_len; ++i) {
                double t = i * dt;

                // 1. Input: L_mt oscillation
                // MATLAB uses cos/sin presets, here we compute on the fly.
                // L_mt_preset = L_mt0 * L_mtn_target + amp * sin(...)
                // Note: MATLAB loop starts at i=2 (t=dt). i=1 is initial.
                // We compute for all t >= 0.

                double L_mt = L_mt_target_val + amp_val * std::sin(2.0 * M_PI * freq * t);

                // 2. Implicit Solver (Newton-Raphson)
                // Skip solver for i=0 (t=0) as it is already equilibrated, unless we want to be safe.
                // MATLAB: i=2 to end.

                if (i > 0) {
                    double L_mnc = L_mn; 
                    int max_iter = 50;
                    double error = 1.0;
                    double delta = 1e-7;
                    int iter = 0;
                    double tol = 1e-8;

                    while (std::abs(error) > tol && iter < max_iter) {
                        iter++;

                        auto calc_error = [&](double l_val) -> double {
                            double L_m = l_val * L_mo;
                            double sin_phi = L_m_height / L_m;
                            if (sin_phi > 1.0) sin_phi = 1.0;
                            double cos_phi = std::cos(std::asin(sin_phi));

                            // Tendon Force
                            double L_t = L_mt - L_m * cos_phi;
                            double tlN = L_t / L_ts; // (1 + eps)
                            double F_tn_val = m.calcfse_(tlN);

                            // Fiber Velocity
                            // MATLAB: V_mn = ((L_mnc - L_mn) * L_mo / dt) / V_mmax;
                            // Note: 'L_mn' here is the state from previous step
                            double V_mn = ((l_val - L_mn) * L_mo / dt) / V_mmax;

                            // Fiber Force
                            double fal = m.calcfal_(l_val);
                            double fpe = m.calcfpe_(l_val);
                            double fv = m.calc_force_velocity_(V_mn, a);

                            double F_fib_norm = a * fal * fv + fpe;
                            return F_tn_val - F_fib_norm * cos_phi;
                        };

                        double err1 = calc_error(L_mnc);

                        // Perturbation
                        double err2 = calc_error(L_mnc + delta);

                        double J = (err2 - err1) / delta;
                        if (std::abs(J) < 1e-14) J = 1e-14;

                        L_mnc -= err1 / J;

                        // Bounds
                        if (L_mnc < 0.01) L_mnc = 0.01;
                        if (L_mnc > 2.0) L_mnc = 2.0;

                        error = err1;
                    }
                    L_mn = L_mnc;
                }

                // 3. Output Calculation
                double L_m_final = L_mn * L_mo;
                double sin_phi = L_m_height / L_m_final;
                if(sin_phi > 1.0) sin_phi = 1.0;
                double cos_phi = std::cos(std::asin(sin_phi));

                double L_t_final = L_mt - L_m_final * cos_phi;
                double tlN_final = L_t_final / L_ts;
                double F_t_final = m.calcfse_(tlN_final) * F_mo;

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