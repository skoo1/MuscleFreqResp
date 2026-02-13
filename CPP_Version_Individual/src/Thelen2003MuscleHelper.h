// Thelen2003MuscleHelper.h
#pragma once
#include "Thelen2003MuscleWrapper.h"

// Inverse Tendon: Force -> Strain (Necessary for Initialization)
// Note: Since we use the exact complex math for calcfse, we need an inverse that matches it.
// However, for simplicity and stability in initialization, we use the standard inverse form here.
// This function is NOT in Thelen2003Muscle.cpp (OpenSim solves equilibrium generally).
// We implement the inverse logic corresponding to the calcfse logic above.
double calc_tendon_strain_from_force(const Thelen2003MuscleWrapper& m, double fse);

// Force-Velocity Multiplier
// Note: Thelen2003Muscle does NOT have a function that returns Force given Velocity.
// It has 'calcdlceN' (Velocity given Force).
// We strictly invert the logic found in 'Thelen2003Muscle::calcdlceN'
double calc_force_velocity(const Thelen2003MuscleWrapper& m, double v_mn, double activation);

// Activation Time Constant
double calc_tau(const Thelen2003MuscleWrapper& m, double a, double u);