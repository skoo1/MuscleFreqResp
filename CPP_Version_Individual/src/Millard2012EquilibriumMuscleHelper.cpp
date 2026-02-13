#include "Millard2012EquilibriumMuscleHelper.h"

// Calculates the tendon length that corresponds to a specific target tendon force.
// Used during initialization to find the starting state of the muscle-tendon unit.
double calc_tendon_length_from_force(const OpenSim::Millard2012EquilibriumMuscle& m, double Ft_target) {
    const auto& fseCurve = m.getTendonForceLengthCurve();
    double Fmo = m.getMaxIsometricForce();
    double Lts = m.getTendonSlackLength();

    // Initial guess for normalized tendon length (slightly stretched)
    double ltN = 1.01;
    double tol = 1e-8;  // Tolerance for convergence
    int maxIter = 50;   // Maximum number of iterations

    for (int i = 0; i < maxIter; ++i) {
        // Calculate tendon force at current length guess
        double f_val = fseCurve.calcValue(ltN);

        // Calculate error: Calculated Force - Target Force
        double err = f_val * Fmo - Ft_target;

        // Check for convergence
        if (std::abs(err) < tol * Fmo) break;

        // Calculate derivative (stiffness) for Newton step
        double df_val = fseCurve.calcDerivative(ltN, 1);
        double derr = df_val * Fmo;

        // Update length guess
        // Avoid division by zero if stiffness is very small
        if (std::abs(derr) < 1e-10) ltN += 0.001;
        else ltN = ltN - err / derr;
    }
    // Return actual tendon length (Normalized Length * Slack Length)
    return ltN * Lts;
}

// Calculates activation time constant (tau) based on excitation (u) and activation (a).
// Used for first-order activation dynamics: da/dt = (u - a) / tau
double get_Tau_MMM(double a, double u) {
    if (u > a) return 0.01 * (0.5 + 1.5 * a); // Activation (fast)
    else       return 0.04 / (0.5 + 1.5 * a); // Deactivation (slow)
}