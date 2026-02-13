#include "OpenSim/Actuators/Thelen2003Muscle.h"

class Thelen2003MuscleWrapper : public OpenSim::Thelen2003Muscle {
public:
    Thelen2003MuscleWrapper() : OpenSim::Thelen2003Muscle() {}

    double calcFm_(double ma, double fal, double fv, 
                        double fpe, double fiso) const;

    double calcActiveFm_(double ma, double fal, 
                        double fv, double fiso) const;

    //Stiffness related functions
    double calcDFmDlce_(double lce, double a,  double fv, 
                        double fiso, double ofl) const;

    double calcDFmATDlce_(double lce, double phi, double cosphi, 
                        double Fm, double d_Fm_d_lce, double penHeight) const;

    double calcDFseDlce_(double tl, double lce, double phi, double cosphi, 
                        double fiso, double tsl, double vol) const;

    double calcDFseDtl_(double tl, double fiso, double tsl) const;
    
    
    //Tendon related helper functions
    double calcfse_(double tlN) const;
    double calcDfseDtlN_(double tlN) const;
    double calcfsefisoPE_(double tlN) const;

    //Active force length functions
    double calcfal_( double lceN) const;
    double calcDfalDlceN_( double lceN) const;

    //Parallel element functions    
    double calcfpe_(double lceN) const;
    double calcDfpeDlceN_(double lceN) const;
    double calcfpefisoPE_(double lceN) const;    

    //Force velocity functions      
    double calcdlceN_(double act,double fal, double actFalFv) const;
    double calcfv_(double aFse, double aFpe, double aFal,
                        double aCosPhi, double aAct) const;
    double calcfvInv_(double aAct,  double aFal, double dlceN, 
                        double tolerance, int maxIterations) const;
    double calcDdlceDaFalFv_(double aAct, double fal, 
                        double aFalFv) const;

    //Returns true if the fiber state is currently clamped to prevent the 
    //fiber from attaining a length that is too short.
    bool isFiberStateClamped_(const SimTK::State& s, double dlceN) const;

	// Force-Velocity Multiplier (Standard Thelen Inversion)
    double calc_force_velocity_(double v_mn, double activation) const;

    // Activation Time Constant
    double Thelen2003MuscleWrapper::calc_tau_(double a, double u) const;

    // Inverse Tendon: Force -> Strain (Necessary for Initialization)
    double Thelen2003MuscleWrapper::calc_tendon_strain_from_force_(double fse) const;
};
