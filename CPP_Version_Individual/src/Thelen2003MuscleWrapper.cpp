#include "Thelen2003MuscleWrapper.h"

using namespace std;
using namespace OpenSim;
using namespace SimTK;

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


//=============================================================================
// TENDON RELATED FUNCTIONS
//=============================================================================
double Thelen2003MuscleWrapper::calcfse_(const double tlN) const 
{
    double x = tlN-1;
    double e0 = get_FmaxTendonStrain();
    
    /*The paper reports etoe = 0.609e0, however, this is a severely rounded off
        The exact answer, to SimTK::Eps is   
        etoe =  99*e0*e^3 / ( 166*e^3 - 67)
        klin =  67 /( 100*(e0 - (99*e0*e^3)/(166*e^3-67)) )
        See thelenINIT_20120127.mw for details
    */    
    double kToe = 3.0;
    double Ftoe = 33.0/100.0;
    double t1   = exp(0.3e1);
    double eToe = (0.99e2*e0*t1) / (0.166e3*t1 - 0.67e2);
    t1 = exp(0.3e1);
    double klin = (0.67e2/0.100e3) 
                * 1.0/(e0 - (0.99e2*e0*t1) / (0.166e3*t1 - 0.67e2));

    //Compute tendon force
    double fse = 0;
    if (x > eToe){
        fse = klin*(x-eToe)+Ftoe;
    }else if (x>0.0){ 
        fse =(Ftoe/(exp(kToe)-1.0))*(exp(kToe*x/eToe)-1.0);
    }else{
        fse=0.;}

    return fse;
}

//==============================================================================
// ACTIVE FORCE LENGTH FUNCTIONS
//==============================================================================
double Thelen2003MuscleWrapper::calcfal_(const double lceN) const{       
    double kShapeActive = get_KshapeActive();   
    double x=(lceN-1.)*(lceN-1.);
    double fal = exp(-x/kShapeActive);
    return fal;
}

//=============================================================================
// FIBER PARALLEL ELEMENT HELPER FUNCTIONS
//=============================================================================
double Thelen2003MuscleWrapper::calcfpe_(const double lceN) const {
    double fpe = 0;
    double e0 = get_FmaxMuscleStrain();
    double kpe = get_KshapePassive();

    //Compute the passive force developed by the muscle
    if(lceN > 1.0){
        double t5 = exp(kpe * (lceN - 0.10e1) / e0);
        double t7 = exp(kpe);
        fpe = (t5 - 0.10e1) / (t7 - 0.10e1);
    }
    return fpe;
}
