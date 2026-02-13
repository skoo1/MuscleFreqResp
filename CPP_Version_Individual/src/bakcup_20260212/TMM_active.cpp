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

double sinwave(double freq, double t, double a0, double a, double b) {
    double val = (a0 - a) / b;
    if (val > 1.0) val = 1.0; if (val < -1.0) val = -1.0;
    double phase = std::asin(val);
    return a + std::sin(2.0 * M_PI * freq * t + phase) * b;
}

// ---------------------------------------------------------------------------
// Main Simulation
// ---------------------------------------------------------------------------
int main() {
    try {
        // --- Configuration ---
        double L_mn_input = 1.0;
        double U_input = 0.5;
        double Amp_input = 0.01;
        double Mass_input = 3e9;
        double SimTime = 120.0;
        double dt = 0.001;
        double FreqLow = 0.1;
        double FreqHigh = 100.0;
        int NumFreqSamples = 100;

        // --- Muscle Object Setup ---
        Thelen2003MuscleWrapper m;

        m.set_max_isometric_force(3549.0);
        m.set_optimal_fiber_length(0.05);
        m.set_tendon_slack_length(0.25);
        m.set_pennation_angle_at_optimal(0.4363);
        m.set_max_contraction_velocity(10.0);

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
        double L_m_height = L_mo * std::sin(alphaOpt);

        // --- Initialization ---
        double L_mn0 = L_mn_input;
        double L_m0 = L_mn0 * L_mo;

        double f_l0 = m.calcfal_(L_mn0);
        double F_pn0 = m.calcfpe_(L_mn0);
        double f_v0 = 1.0;

        double alpha0 = std::asin(L_m_height / L_m0);

        double F_tn_init = (U_input * f_l0 * f_v0 + F_pn0) * std::cos(alpha0);
        double eps_t_init = m.calc_tendon_strain_from_force_(F_tn_init);
        double L_t0 = L_ts * (1.0 + eps_t_init);

        double L_mt0 = L_m0 * std::cos(alpha0) + L_t0;
        double F_ext_equil = F_mo * F_tn_init;

        std::cout << "Init: L_mt0=" << L_mt0 << ", F_equil=" << F_ext_equil << std::endl;

        // --- Frequency Sweep ---
        std::vector<double> frequencies = logspace(std::log10(FreqLow), std::log10(FreqHigh), NumFreqSamples);
        fs::create_directories("TMM_Active_results_csv");

        std::ofstream freqFile("TMM_Active_results_csv/frequency_list.csv");
        if (freqFile.is_open()) {
            for (size_t k = 0; k < frequencies.size(); ++k) {
                freqFile << k << "," << frequencies[k] << "\n";
            }
            freqFile.close();
        }

        for (int k = 0; k < NumFreqSamples; ++k) {
            double freq = frequencies[k];
            std::ofstream outFile("TMM_Active_results_csv/freq_res_" + std::to_string(k) + ".csv");
            outFile << std::setprecision(16);
            outFile << "time,u,F_m_AT\n";

            double L_mt = L_mt0;
            double V_mt = 0.0;
            double L_mn = L_mn0;
            double a = U_input;
            int time_len = (int)(SimTime / dt);

            for (int i = 0; i < time_len; ++i) {
                double t = i * dt;

                // 1. Input
                double u = sinwave(freq, t, U_input, U_input, Amp_input);
                if (u < 0) u = 0; if (u > 1) u = 1;

                // 2. Activation
                if (i > 0) {
                    double tau = m.calc_tau_(a, u);
                    double da_dt = (u - a) / tau;
                    a += dt * da_dt;
                }

                // 3. Implicit Solver
                double L_mnc = L_mn;
                int max_iter = 50;
                double error = 1.0;
                double delta = 1e-7;
                int iter = 0;

                while (std::abs(error) > 1e-8 && iter < max_iter) {
                    iter++;

                    auto calc_error = [&](double l_val) -> double {
                        double L_m = l_val * L_mo;
                        double sin_phi = L_m_height / L_m;
                        if (sin_phi > 1.0) sin_phi = 1.0;
                        double cos_phi = std::cos(std::asin(sin_phi));

                        double L_t = L_mt - L_m * cos_phi;
                        
                        // OpenSim's calcfse takes 'tlN' (Normalized Tendon Length)
                        // tlN = L_t / L_ts = 1.0 + eps_t
                        double tlN = L_t / L_ts;
                        
                        double F_tn_val = m.calcfse_(tlN);

                        double V_mn = ((l_val - L_mn) * L_mo / dt) / V_mmax;

                        double fal = m.calcfal_(l_val);
                        double fpe = m.calcfpe_(l_val);
                        double fv = m.calc_force_velocity_(V_mn, a);

                        double F_fib_norm = a * fal * fv + fpe;
                        return F_tn_val - F_fib_norm * cos_phi;
                    };

                    double err1 = calc_error(L_mnc);
                    double err2 = calc_error(L_mnc + delta);

                    double J = (err2 - err1) / delta;
                    if (std::abs(J) < 1e-14) J = 1e-14;

                    L_mnc -= err1 / J;
                    if (L_mnc < 0.01) L_mnc = 0.01;
                    if (L_mnc > 2.0) L_mnc = 2.0;

                    error = err1;
                }
                L_mn = L_mnc;

                // 4. Output Force
                double L_m_final = L_mn * L_mo;
                double sin_phi = L_m_height / L_m_final;
                if (sin_phi > 1.0) sin_phi = 1.0;
                double cos_phi = std::cos(std::asin(sin_phi));
                double L_t_final = L_mt - L_m_final * cos_phi;
                
                // Pass tlN (Normalized Length) to matching calcfse_
                double tlN_final = L_t_final / L_ts; 
                double F_t_final = m.calcfse_(tlN_final) * F_mo;

                // 5. External Dynamics
                double A_mt = (F_ext_equil - F_t_final) / Mass_input;
                V_mt += A_mt * dt;
                L_mt += V_mt * dt;

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