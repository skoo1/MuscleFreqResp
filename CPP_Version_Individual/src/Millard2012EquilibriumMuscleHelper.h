#include <OpenSim/Actuators/Millard2012EquilibriumMuscle.h>

// Calculates the tendon length that corresponds to a specific target tendon force.
// Used during initialization to find the starting state of the muscle-tendon unit.
double calc_tendon_length_from_force(const OpenSim::Millard2012EquilibriumMuscle& m, double Ft_target);

// Calculates activation time constant (tau) based on excitation (u) and activation (a).
// Used for first-order activation dynamics: da/dt = (u - a) / tau
double get_Tau_MMM(double a, double u);