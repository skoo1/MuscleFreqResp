// By Minseung Kim, Seungwoo Yoon and Seungbum Koo
// KAIST, Daejeon, South Korea
// February 23, 2026

#include "Thelen2003MuscleHelper.h"

// Inverse Tendon: Force -> Strain (Necessary for Initialization)
// Note: Since we use the exact complex math for calcfse, we need an inverse that matches it.
// However, for simplicity and stability in initialization, we use the standard inverse form here.
// This function is NOT in Thelen2003Muscle.cpp (OpenSim solves equilibrium generally).
// We implement the inverse logic corresponding to the calcfse logic.
double calc_tendon_strain_from_force(const Thelen2003MuscleWrapper& m, double fse) 
{
    // Re-calculate constants exactly as in calcfse
    double e0 = m.get_FmaxTendonStrain();
    double kToe = 3.0;
    double Ftoe = 33.0 / 100.0;
    double t1 = exp(0.3e1);
    double eToe = (0.99e2 * e0 * t1) / (0.166e3 * t1 - 0.67e2);
    t1 = exp(0.3e1);
    double klin = (0.67e2 / 0.100e3) * 1.0 / (e0 - (0.99e2 * e0 * t1) / (0.166e3 * t1 - 0.67e2));

    if (fse <= 0.0) return 0.0;

    if (fse <= Ftoe) {
        // Inverse of: fse = (Ftoe/(exp(kToe)-1.0))*(exp(kToe*x/eToe)-1.0)
        return (eToe / kToe) * log(fse * (exp(kToe) - 1.0) / Ftoe + 1.0);
    } else {
        // Inverse of: fse = klin*(x-eToe)+Ftoe
        return (fse - Ftoe) / klin + eToe;
    }
}

// Force-Velocity Multiplier
// Note: Thelen2003Muscle does NOT have a function that returns Force given Velocity.
// It has 'calcdlceN' (Velocity given Force).
// We strictly invert the logic found in 'Thelen2003Muscle::calcdlceN'
double calc_force_velocity(const Thelen2003MuscleWrapper& m, double v_mn, double activation) 
{
    double aF = m.get_Af();
    double fLen = m.get_Flen();
    double v_scale = 1.0 / (0.25 + 0.75 * activation);
    double v_norm_scaled = v_mn * v_scale;
    double fv = 0;

    if (v_mn <= 0) { // Concentric
        fv = (1. + v_norm_scaled) / (1. - v_norm_scaled / aF);
    } else { // Eccentric
        double someth = v_norm_scaled * (2. + 2. / aF) / (fLen - 1.);
        fv = (1. + fLen * someth) / (1. + someth);
    }
    return fv;
}

// Returns Activation/Deactivation Time Constant
// from excitation level (u) and activation level (a)
double calc_tau(const Thelen2003MuscleWrapper& m, double a, double u) 
{
    if (u > a) return m.get_activation_time_constant() * (0.5 + 1.5 * a);
    else       return m.get_deactivation_time_constant() / (0.5 + 1.5 * a);
}