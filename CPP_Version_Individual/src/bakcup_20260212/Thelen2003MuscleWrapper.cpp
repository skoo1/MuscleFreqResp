#include "Thelen2003MuscleWrapper.h"

using namespace std;
using namespace OpenSim;
using namespace SimTK;

double Thelen2003MuscleWrapper::calcFm_(double ma, double fal, double fv, double fpe,
                                 double fiso) const
{
    double Fm = (ma*fal*fv + fpe)*fiso;
    return Fm;
}

double Thelen2003MuscleWrapper::calcActiveFm_(double ma, double fal, double fv,
                                 double fiso) const
{
    double aFm = (ma*fal*fv)*fiso;
    return aFm;
}

double Thelen2003MuscleWrapper::calcDFmDlce_(double lce, double ma, double fv, 
                                      double fiso, double ofl) const
{

            double lceN = lce/ofl;
            double dfal_d_lceN = calcDfalDlceN_(lceN);
            double dfpe_d_lceN = calcDfpeDlceN_(lceN);

            double dFm_d_lce = ((ma*fv)*dfal_d_lceN + dfpe_d_lceN)*fiso
                              *(1/ofl);                   
            return dFm_d_lce;
}

double Thelen2003MuscleWrapper::calcDFmATDlce_(double lce, double phi, double cosphi, 
    double Fm, double d_Fm_d_lce, double penHeight) const
{
    //SINGULARITY: when vol*vol/(lce*lce) = 1,same as phi=pi/2
    double tmp1 = penHeight*penHeight;
    double tmp2 = lce*lce;
    double tmp3 = tmp2*lce;
    double dcosphi_d_lce = (tmp1 /(tmp3*pow((1-(tmp1/tmp2)),0.5 ) ));


    double d_FmAT_d_lce = d_Fm_d_lce*cosphi + Fm*dcosphi_d_lce;
    return d_FmAT_d_lce;
}

double Thelen2003MuscleWrapper::calcDFseDlce_(double tl, double lce, double phi, double cosphi,
                                      double fiso, double tsl, double vol) const
{
    double tlN = tl/tsl;
        //SINGULARITY: When lce = 0
    double tmp1 = vol/lce;        

    //SINGULARITY: when vol/lce = 1 - equivalent to when phi = pi/2
    double dphi_d_lce = -vol / ( lce*lce * pow( (1-(tmp1)*(tmp1)),0.5)); 
    double dtl_d_lce  = -cos(phi) + lce*sin(phi)*dphi_d_lce;

    double dfse_d_tlN  = calcDfseDtlN_(tlN); 
            tmp1 = (fiso/tsl);
    double dFt_d_lce = dfse_d_tlN*dtl_d_lce*tmp1;  
    return dFt_d_lce;
}

double Thelen2003MuscleWrapper::calcDFseDtl_(double tl, double fiso, double tsl) const
{
    double dfse_d_tlN  = calcDfseDtlN_(tl/tsl);
    double tmp1 = (fiso/tsl);
    double dFt_d_tl= dfse_d_tlN*tmp1;
    return dFt_d_tl;
}


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


double Thelen2003MuscleWrapper::calcDfseDtlN_(const double tlN) const {
    double x = tlN-1;
    double e0 = get_FmaxTendonStrain();
    
    /*The paper reports etoe = 0.609e0, however, this is a severely rounded off
    result of the exact answer:    
        etoe =  99*e0*e^3 / ( 166*e^3 - 67)
        See thelenINIT_20120127.mw for details
    */
    double kToe = 3.0;
    double Ftoe = 33.0/100.0;

    double t1   = exp(0.3e1);
    double eToe = (0.99e2*e0*t1) / (0.166e3*t1 - 0.67e2);

    /*The paper reports etoe = 0.609e0, however, this is a severely rounded off
    result of the exact answer:    
        klin =  67 /( 100*(e0 - (99*e0*e^3)/(166*e^3-67)) )
        See thelenINIT_20120127.mw for details
    */
    t1 = exp(0.3e1);
    double klin = (0.67e2/0.100e3) 
                * 1.0/(e0 - (0.99e2*e0*t1) / (0.166e3*t1 - 0.67e2));

    //Compute tendon force
    double dfse_d_dtlN = 0;
    if (x > eToe){
        dfse_d_dtlN = klin;
    }else if (x>0.0){ 
        dfse_d_dtlN =(Ftoe/(exp(kToe)-1.0)) * (kToe/eToe) * (exp(kToe*x/eToe));
    }else{
        dfse_d_dtlN=0.;}

    return dfse_d_dtlN;
}

double Thelen2003MuscleWrapper::calcfsefisoPE_(double tendonStrain) const
{

    double tendon_strain =  tendonStrain;
    double fmaxTendonStrain = get_FmaxTendonStrain();       

    //Future optimization opportunity: precompute kToe, fToe, eToe and klin
    //when the muscle is initialized. Store these values rather than 
    //computing them every time.

    double kToe = 3.0;
    double Ftoe = 33.0/100.0;

    double t1   = exp(0.3e1);
    double eToe = (0.99e2*fmaxTendonStrain*t1) / (0.166e3*t1 - 0.67e2);

    t1 = exp(0.3e1);
    double klin = (0.67e2/0.100e3) 
                * 1.0/(fmaxTendonStrain - (0.99e2*fmaxTendonStrain*t1) 
                / (0.166e3*t1 - 0.67e2));

    //Compute the energy stored in the tendon. 
    //Integrals computed symbolically in muscle_kepew_20111021.mw just to check
    double tendonPE = 0.0;
    double lenR        = getTendonSlackLength();
    double lenTdn    = (tendon_strain+1)*lenR;
    double lenToe    = (eToe+1.0)*lenR;    
    double fiso        = getMaxIsometricForce();

    if (tendon_strain>eToe){
       //compute the energy stored in the toe portion of the tendon strain curve
        double len = lenToe;
        double toePE_len = (fiso*Ftoe/(exp(kToe)-1.0))
                            *((lenR*eToe/kToe)
                            *exp(kToe*(len-lenR)/(lenR*eToe)) - len);
        len =  lenR;
        double toePE_0    =  (fiso*Ftoe/(exp(kToe)-1.0))
                            *((lenR*eToe/kToe)
                            *exp(kToe*(len-lenR)/(lenR*eToe)) - len);
        // double toePEtest = toePE_len-toePE_0;

        //compute the energy stored in the linear section of the 
        //tendon strain curve from ..... 0 to len
        len = lenTdn;
        double linPE_len = (1.0/2.0)*(fiso*klin*(len*len)/lenR) 
                           + fiso*len*(klin*(-1.0-eToe)+Ftoe);
        //ditto from 0 .... eToe
        len = lenToe;
        double linPE_eToe= (1.0/2.0)*(fiso*klin*(len*len)/lenR) 
                            + fiso*len*(klin*(-1.0-eToe)+Ftoe);       
        
        //compute the total potential energy stored in the tendon
         tendonPE =(toePE_len-toePE_0) + (linPE_len-linPE_eToe);
    }else if (tendon_strain>0.0){ 
        //PE from 0 .... len
        double len = lenTdn;
        double toePE_len = (fiso*Ftoe/(exp(kToe)-1.0)) * ((lenR*eToe/kToe) 
            * exp(kToe*(len-lenR)/(lenR*eToe)) - len);
        //PE from 0 .... eToe
        len = lenR;
        double toePE_0    =  (fiso*Ftoe/(exp(kToe)-1.0)) * ((lenR*eToe/kToe) 
            * exp(kToe*(len-lenR)/(lenR*eToe)) - len);

        //Compute the total PE stored in the tendon
        tendonPE = toePE_len-toePE_0;
    }else{
        tendonPE = 0.0;
    }
    
    
    return tendonPE;
}

//==============================================================================
//
// ACTIVE FORCE LENGTH FUNCTIONS
//
//==============================================================================
double Thelen2003MuscleWrapper::calcfal_(const double lceN) const{       
    double kShapeActive = get_KshapeActive();   
    double x=(lceN-1.)*(lceN-1.);
    double fal = exp(-x/kShapeActive);
    return fal;
}
double Thelen2003MuscleWrapper::calcDfalDlceN_(const double lceN) const {
    double kShapeActive = get_KshapeActive();   
    double t1 = lceN - 0.10e1;
    double t2 = 0.1e1 / kShapeActive;
    double t4 = t1 * t1;
    double t6 = exp(-t4 * t2);
    double dfal_d_lceN = -0.2e1 * t1 * t2 * t6;
    return dfal_d_lceN;
}

//=============================================================================
//
// FIBER PARALLEL ELEMENT HELPER FUNCTIONS
//
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

double Thelen2003MuscleWrapper::calcDfpeDlceN_(const double lceN) const {
    double dfpe_d_lceN = 0;
    double e0 = get_FmaxMuscleStrain();
    double kpe = get_KshapePassive();

    if(lceN > 1.0){
        double t1 = 0.1e1 / e0;
        double t6 = exp(kpe * (lceN - 0.10e1) * t1);
        double t7 = exp(kpe);
        dfpe_d_lceN = kpe * t1 * t6 / (t7 - 0.10e1);
    }
    return dfpe_d_lceN;
}

double Thelen2003MuscleWrapper::calcfpefisoPE_(double lceN) const
{
    double fmaxMuscleStrain = get_FmaxMuscleStrain();
    double kShapePassive = get_KshapePassive();

    double musclePE = 0.0;
    //Compute the potential energy stored in the muscle
    if(lceN > 1.0){
        //Shorter variable names to make the equations readable.
        double lenR = getOptimalFiberLength();
        double fiso = getMaxIsometricForce();
        double len = lceN*lenR;        
        double kpe = kShapePassive;
        double e0 = fmaxMuscleStrain;


        //PE storage at current stretch
        double fpePE_len = (fiso/(exp(kpe)-1))
            *( (lenR*e0/kpe)*exp( (kpe/e0)*( (len/lenR)-1)) - len); 
        
        //PE stored between 0 and 1 for the exponential function that is 
        //used to represent fpe for normalized fiber lengths > 1.
        len = lenR;
        double fpePE_0 = (fiso/(exp(kpe)-1))
            *( (lenR*e0/kpe)*exp( (kpe/e0)*( (len/lenR)-1)) - len); 

        musclePE = fpePE_len - fpePE_0;

    }else{
        musclePE = 0.0;
    }
    return musclePE;
}


//=============================================================================
//
// FIBER FORCE - VELOCITY CURVE FUNCTIONS
//
//=============================================================================

double Thelen2003MuscleWrapper::calcdlceN_(double act,double fal,double actFalFv) const
{
    //The variable names have all been switched to closely match 
    //with the notation in Thelen 2003.
    double dlceN = 0.0;      //contractile element velocity    
    double af   = get_Af();

    double a    = act;
    double afl  = a*fal; //afl = a*fl
    double Fm   = actFalFv;     //Fm = a*fl*fv    
    double flen = get_Flen();
    // double Fmlen_afl = flen*afl;

    double dlcedFm = 0.0; //partial derivative of contractile element
                          // velocity w.r.t. Fm

    double b = 0;
    double db= 0;

    double Fm_asyC = 0;           //Concentric contraction asymptote
    double Fm_asyE = afl*flen;    
                                //Eccentric contraction asymptote
    double asyE_thresh = get_fv_linear_extrap_threshold();

    //If fv is in the appropriate region, use 
    //Thelen 2003 Eqns 6 & 7 to compute dlceN
    if (Fm > Fm_asyC && Fm < Fm_asyE*asyE_thresh){

        if( Fm <= afl ){        //Muscle is concentrically contracting
            b = afl + Fm/af;
            db= 1/af;
        }else{                    //Muscle is eccentrically contracting
            b = ((2+2/af)*(afl*flen-Fm))/(flen-1); 
            db= ((2+2/af)*(-1))/(flen-1); 
        }

        dlceN = (0.25 + 0.75*a)*(Fm-afl)/b; 
        //Scaling by VMAX is left out, and is post multiplied outside 
        //of the function


    }else{  //Linear extrapolation
            double Fm0 = 0.0; //Last Fm value from the Thelen curve

            //Compute d and db/dFm from Eqn 7. of Thelen2003
            //for the last
            if(Fm <= Fm_asyC){ //Concentrically contracting
                Fm0 = Fm_asyC;
                b = afl + Fm0/af;
                db= 1/af;               
            }else{             //Eccentrically contracting
                Fm0 = asyE_thresh*Fm_asyE;
                b = ((2+2/af)*(afl*flen-Fm0))/(flen-1); 
                db= ((2+2/af)*(-1))/(flen-1); 
            }

            //Compute the last dlceN value that falls in the region where
            //Thelen 2003 Eqn. 6 is valid
            double dlce0 = (0.25 + 0.75*a)*(Fm0-afl)/b;

            //Compute the dlceN/dfm of Eqn. 6 of Thelen 2003 at the last
            //valid point
            dlcedFm = (0.25 + 0.75*a)*(1)/b 
                    - ((0.25 + 0.75*a)*(Fm0-afl)/(b*b))*db;

            //Linearly extrapolate Eqn. 6 from Thelen 2003 to compute
            //the new value for dlceN/dFm
            dlceN = dlce0 + dlcedFm*(Fm-Fm0);            
        }
            
        return dlceN;
}

double Thelen2003MuscleWrapper::calcfv_(double aFse,     double aFpe, double aFal, 
                                        double aCosPhi,  double aAct) const
{
    //This only works for an equilibrium model, but its a lot less 
    //computationally expensive (and error prone) than trying to invert the 
    //weird function that defines the fv curve in the Thelen model
    double fv = ((aFse/aCosPhi) - aFpe)/(aAct*aFal);
    return fv;
}

double Thelen2003MuscleWrapper::calcDdlceDaFalFv_(double aAct, 
                                        double aFal, double aFalFv) const
{
    //The variable names have all been switched to closely match with 
    //the notation in Thelen 2003.
    // double dlceN = 0.0;      //contractile element velocity    
    double af   = get_Af();

    double a    = aAct;
    double afl  = aAct*aFal;  //afl = a*fl
    double Fm   = aFalFv;    //Fm = a*fl*fv    
    double flen = get_Flen();
    // double Fmlen_afl = flen*aAct*aFal;

    double dlcedFm = 0.0; //partial derivative of contractile element 
                          //velocity w.r.t. Fm

    double b = 0;
    double db= 0;

    double Fm_asyC = 0;           //Concentric contraction asymptote
    double Fm_asyE = aAct*aFal*flen;    
                                //Eccentric contraction asymptote
    double asyE_thresh = get_fv_linear_extrap_threshold();

    //If fv is in the appropriate region, use 
    //Thelen 2003 Eqns 6 & 7 to compute dlceN
    if (Fm > Fm_asyC && Fm < Fm_asyE*asyE_thresh){

        if( Fm <= afl ){        //Muscle is concentrically contracting
            b = afl + Fm/af;
            db= 1/af;
        }else{                    //Muscle is eccentrically contracting
            b = ((2+2/af)*(afl*flen-Fm))/(flen-1); 
            db= ((2+2/af)*(-1))/(flen-1); 
        }

        //This variable may have future use outside this function
        dlcedFm = (0.25 + 0.75*a)*(1)/b - ((0.25 + 0.75*a)*(Fm-afl)/(b*b))*db;            

    }else{  //Linear extrapolation
            double Fm0 = 0.0; //Last Fm value from the Thelen curve

            //Compute d and db/dFm from Eqn 7. of Thelen2003
            //for the last
            if(Fm <= Fm_asyC){ //Concentrically contracting
                Fm0 = Fm_asyC;
                b = afl + Fm0/af;
                db= 1/af;               
            }else{             //Eccentrically contracting
                Fm0 = asyE_thresh*Fm_asyE;
                b = ((2+2/af)*(afl*flen-Fm0))/(flen-1); 
                db= ((2+2/af)*(-1))/(flen-1); 
            }

            
            //Compute the dlceN/dfm of Eqn. 6 of Thelen 2003 at the last
            //valid point
            dlcedFm = (0.25 + 0.75*a)*(1)/b 
                - ((0.25 + 0.75*a)*(Fm0-afl)/(b*b))*db;
          
        }
            
        return dlcedFm;
}

// Compute the force-velocity multiplier by inverting Thelen 2003' f-v
// equations for fiber-velocity given the active fiber force (see calcdlceN()).
// This is here because it is non-trivial to correctly invert the piece-wise
// continuous force velocity curve specified by the modified Thelen2003Muscle
// force velocity curve. This converges quickly and is well tested.
double Thelen2003MuscleWrapper::calcfvInv_(double aAct,double aFal,double dlceN,
                                    double tolerance, int maxIterations) const
{
    double result = SimTK::NaN;
    double ferr=1;
    double iter= 0;

    double dlceN1 = 0;
    double dlceN1_d_Fm = 0;
    double fv = 1;
    double aFalFv = fv*aAct*aFal;
    double delta_aFalFv = 0;

    while(abs(ferr) >= tolerance && iter < maxIterations)
    {
        dlceN1 = calcdlceN_(aAct,aFal, aFalFv);
        ferr   = dlceN1-dlceN;
        dlceN1_d_Fm = calcDdlceDaFalFv_(aAct,aFal,aFalFv);

        if(abs(dlceN1_d_Fm) > SimTK::SignificantReal){
           delta_aFalFv = -ferr/(dlceN1_d_Fm);
           aFalFv = aFalFv + delta_aFalFv;
        }
        iter = iter+1;
    }

    if(abs(ferr) < tolerance){
        result = max(0.0, aFalFv/(aAct*aFal));
        return result;
    }

    OPENSIM_THROW_FRMOBJ(OpenSim::Exception, 
        "Solver for force-velocity multiplier failed to converge.");
}

bool Thelen2003MuscleWrapper::isFiberStateClamped_(const SimTK::State& s, double dlceN) const
{
    bool clamped = false;

    //Is the fiber length  clamped and it is shortening, then the fiber length
    //not valid
    if( (getStateVariableValue(s, STATE_FIBER_LENGTH_PATH)
            <= getMinimumFiberLength())
        && dlceN <= 0){
        clamped = true;
    }

    return clamped;     
}

// Force-Velocity Multiplier
// Note: Thelen2003Muscle does NOT have a function that returns Force given Velocity.
// It has 'calcdlceN' (Velocity given Force).
// We strictly invert the logic found in 'Thelen2003Muscle::calcdlceN'
double Thelen2003MuscleWrapper::calc_force_velocity_(double v_mn, double activation) const 
{
    double aF = get_Af();
    double fLen = get_Flen();
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

// Activation Time Constant
double Thelen2003MuscleWrapper::calc_tau_(double a, double u) const {
    if (u > a) return get_activation_time_constant() * (0.5 + 1.5 * a);
    else       return get_deactivation_time_constant() / (0.5 + 1.5 * a);
}

// Inverse Tendon: Force -> Strain (Necessary for Initialization)
// Note: Since we use the exact complex math for calcfse, we need an inverse that matches it.
// However, for simplicity and stability in initialization, we use the standard inverse form here.
// This function is NOT in Thelen2003Muscle.cpp (OpenSim solves equilibrium generally).
// We implement the inverse logic corresponding to the calcfse logic above.
double Thelen2003MuscleWrapper::calc_tendon_strain_from_force_(double fse) const {
    // Re-calculate constants exactly as in calcfse
    double e0 = get_FmaxTendonStrain();
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