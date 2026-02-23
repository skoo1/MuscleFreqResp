// By Minseung Kim, Seungwoo Yoon and Seungbum Koo
// KAIST, Daejeon, South Korea
// February 23, 2026

// Thelen2003MuscleWrapper.h
#pragma once
#include "OpenSim/Actuators/Thelen2003Muscle.h"
#include <cmath>
#include <vector>

// ============================================================================
// [Wrapper Class Explanation]
// The Thelen2003MuscleWrapper class inherits from OpenSim::Thelen2003Muscle.
// 
// Purpose:
// The original force calculation functions (calcfse, calcfal, calcfpe) in 
// OpenSim::Thelen2003Muscle are declared as 'private' or 'protected', 
// making them inaccessible to external custom solvers.
// 
// Implementation Note:
// To ensure results are mathematically identical to OpenSim, the code for 
// these functions was COPIED VERBATIM (directly copy-pasted) from the 
// original 'Thelen2003Muscle.cpp' source file. 
//
// The only modifications are:
// 1. Changing access level to 'public'.
// 2. Appending an underscore ('_') to the function names to avoid conflict.
// ============================================================================


class Thelen2003MuscleWrapper : public OpenSim::Thelen2003Muscle {
public:
    Thelen2003MuscleWrapper() : OpenSim::Thelen2003Muscle() {}

    //Tendon related helper functions
    double calcfse_(double tlN) const;

    //Active force length functions
    double calcfal_( double lceN) const;

    //Parallel element functions    
    double calcfpe_(double lceN) const;
};
