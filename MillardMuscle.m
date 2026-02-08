% By Minseung Kim and Seungbum Koo
% KAIST, Daejeon, South Korea
% February 8, 2026

% This class requires Millard2012 Matlab library
% https://github.com/mjhmilla/Millard2012EquilibriumMuscleMatlabPort

classdef MillardMuscle < MuscleModel
    properties
        MMM             % Struct for time constants
        MuscleArch
        NormMuscleCurves
        ModelConfig
        MTInfo
    end
    
    methods
        function obj = MillardMuscle(name, f_mo, l_mo, l_ts, alphaOpt, mass, damping)
            obj@MuscleModel(name, f_mo, l_mo, l_ts, alphaOpt, mass, damping);
            obj.init_parameters();
        end
        
        function [L_mt, F_equil, Alpha] = initialize_static_given_Lm(obj, a, L_m)
            Alpha = obj.calc_pennation_angle(L_m);
            L_m_AT = L_m*cos(Alpha);
            L_mt_init = L_m_AT + obj.L_ts;
            pathState = [0; L_mt_init];
            muscleState = [0; L_m_AT];

            mtInfo = calcMillard2012DampedEquilibriumMuscleInfo( ...
                a, pathState, muscleState, ...
                obj.MuscleArch, obj.NormMuscleCurves, obj.ModelConfig);

            % mtInfo.initialization.err should be large because of the 
            % difference between F_m_AT and F_t

            F_m_AT  = mtInfo.muscleDynamicsInfo.fiberForceAlongTendon;
            F_mn_AT = F_m_AT / obj.MuscleArch.fiso;

            % Calculate the required normalized tendon length (ltN = lt / L_ts) using the inverse function.
            % Pass the pre-generated inverse tendon curve into the calcBezierYFcnXDerivative function.
            L_tn = calcBezierYFcnXDerivative(F_mn_AT, obj.NormMuscleCurves.tendonInverseCurve, 0);
            L_t  = L_tn * obj.MuscleArch.tendonSlackLength;

            pathState = [0; L_m_AT + L_t];
            muscleState = [0; L_m_AT];

            mtInfo = calcMillard2012DampedEquilibriumMuscleInfo( ...
                a, pathState, muscleState, ...
                obj.MuscleArch, obj.NormMuscleCurves, obj.ModelConfig);

            % mtInfo.initialization.err should be almost zero

            F_equil = mtInfo.muscleDynamicsInfo.fiberForceAlongTendon;
            L_mt    = L_m_AT + L_t;

            % update muscle state
            obj.MTInfo = mtInfo;
            obj.a    = a;
            obj.L_m  = mtInfo.muscleLengthInfo.fiberLength;
            obj.V_m  = mtInfo.fiberVelocityInfo.fiberVelocity;
            obj.L_mt = mtInfo.muscleLengthInfo.fiberLengthAlongTendon + mtInfo.muscleLengthInfo.tendonLength;
            obj.V_mt = mtInfo.fiberVelocityInfo.fiberVelocityAlongTendon + mtInfo.fiberVelocityInfo.tendonVelocity;
            obj.L_mn = mtInfo.muscleLengthInfo.normFiberLength;
            obj.F_t  = mtInfo.muscleDynamicsInfo.tendonForce;
            obj.F_m  = mtInfo.muscleDynamicsInfo.fiberForce;
        end

        function [L_m, F_t, Alpha] = initialize_static_given_Lmt(obj, a, L_mt)
            % Initially V_mt = 0
            V_mt = 0;
            pathState   = [V_mt; L_mt];

            low_L_m_AT          = 0.0;
            high_L_m_AT         = L_mt - 1e-5;
            error1              = 1000;

            while abs(real(error1)) > 1e-8
                L_m_AT = 1/2 * (low_L_m_AT + high_L_m_AT);
                muscleState = [0, L_m_AT]; % input V_m_AT and L_m_AT

                % the following calculates the error F_m_AT - F_t.
                mtInfo = calcMillard2012DampedEquilibriumMuscleInfo( ...
                    a, pathState, muscleState, ...
                    obj.MuscleArch, obj.NormMuscleCurves, obj.ModelConfig);

                F_m_AT = mtInfo.muscleDynamicsInfo.fiberForceAlongTendon;
                F_t    = mtInfo.muscleDynamicsInfo.tendonForce;

                error1 = F_t - F_m_AT;

                if (error1 < 0)
                    high_L_m_AT = L_m_AT;
                else
                    low_L_m_AT = L_m_AT;
                end
            end

            L_m_height = obj.L_mo * sin(obj.AlphaOpt);
            L_m = sqrt(L_m_AT^2 + L_m_height^2);
            Alpha = obj.calc_pennation_angle(L_m);

            obj.MTInfo = mtInfo;
            obj.a    = a;
            obj.L_m  = mtInfo.muscleLengthInfo.fiberLength;
            obj.V_m  = mtInfo.fiberVelocityInfo.fiberVelocity;
            obj.L_mt = mtInfo.muscleLengthInfo.fiberLengthAlongTendon + mtInfo.muscleLengthInfo.tendonLength;
            obj.V_mt = mtInfo.fiberVelocityInfo.fiberVelocityAlongTendon + mtInfo.fiberVelocityInfo.tendonVelocity;
            obj.L_mn = mtInfo.muscleLengthInfo.normFiberLength;
            obj.F_t  = mtInfo.muscleDynamicsInfo.tendonForce;
            obj.F_m  = mtInfo.muscleDynamicsInfo.fiberForce;
        end

        function F_t = updateDynamics(obj, dt, u, F_ext_equil, mass_ext, damping_ext)
            % Update Activation
            da_dt = obj.getActivationRate(u, obj.a);
            a = obj.a + dt * da_dt;
 
            % Calculate Dynamics using Millard API
            L_m    = obj.L_m;
            Alpha = obj.calc_pennation_angle(L_m);
            L_m_AT = L_m * cos(Alpha);
            pathState = [obj.V_mt; obj.L_mt];
            muscleState = L_m_AT; % Scalar input means compute derivative
            
            mtInfo = calcMillard2012DampedEquilibriumMuscleInfo(...
                a, pathState, muscleState, ...
                obj.MuscleArch, obj.NormMuscleCurves, obj.ModelConfig);
            
            F_t = mtInfo.muscleDynamicsInfo.tendonForce;

            % Update muscle length
            V_m_AT = mtInfo.state.derivative;
            L_m_AT = L_m_AT + V_m_AT * dt;
            L_m_height = obj.L_mo * sin(obj.AlphaOpt);
            L_m = sqrt(L_m_AT^2 + L_m_height^2);

            % Equilibrium with external forces
            A_mt    = (F_ext_equil - F_t - damping_ext * obj.V_mt) / mass_ext;
            V_mt    = obj.V_mt + A_mt * dt;
            L_mt    = obj.L_mt + V_mt * dt;

            obj.MTInfo = mtInfo;
            obj.a    = a;
            obj.L_m  = L_m;
            obj.V_m  = mtInfo.fiberVelocityInfo.fiberVelocity;
            obj.L_mt = L_mt;
            obj.V_mt = V_mt;
            obj.L_mn = mtInfo.muscleLengthInfo.normFiberLength;
            obj.F_t  = mtInfo.muscleDynamicsInfo.tendonForce;
            obj.F_m  = mtInfo.muscleDynamicsInfo.fiberForce;
        end

        function F_t = updateDynamicsQuasiStatic(obj, dt, u, L_mt, V_mt)
            obj.L_mt = L_mt;
            obj.V_mt = V_mt; % quasi-static condition
            
            % Activation dynamics
            da_dt = obj.getActivationRate(u, obj.a);
            a = obj.a + dt * da_dt;

            % Kinematically driven
            pathState   = [V_mt; L_mt];
            
            % Initial value
            L_m_AT      = obj.L_m * cos(obj.calc_pennation_angle(obj.L_m));
            L_m_AT_old  = L_m_AT;

            count       = 0;
            max_iter    = 50;
            err         = 1.0;
            delta       = 1e-7;
    
            while (abs(real(err)) > 1e-8 && count < max_iter)
                count   = count + 1;

                V_m_AT = (L_m_AT - L_m_AT_old)/dt;
                muscleState = [V_m_AT; L_m_AT]; % these two should be found
                
                % the following calculates force balance error F_m_AT - F_t
                mtInfo = calcMillard2012DampedEquilibriumMuscleInfo( ...
                    a, pathState, muscleState, ...
                    obj.MuscleArch, obj.NormMuscleCurves, obj.ModelConfig);

                F_m_AT = mtInfo.muscleDynamicsInfo.fiberForceAlongTendon;
                F_t    = mtInfo.muscleDynamicsInfo.tendonForce;
                err  = F_t - F_m_AT;  % error of force equilbrium
    
                L_m_AT_p = L_m_AT + delta;
                V_m_AT_p = (L_m_AT_p - L_m_AT_old)/dt;

                muscleState_p = [V_m_AT_p; L_m_AT_p]; % these two should be found

                mtInfo_p = calcMillard2012DampedEquilibriumMuscleInfo( ...
                    a, pathState, muscleState_p, ...
                    obj.MuscleArch, obj.NormMuscleCurves, obj.ModelConfig);
    
                F_m_AT_p = mtInfo_p.muscleDynamicsInfo.fiberForceAlongTendon;
                F_t_p    = mtInfo_p.muscleDynamicsInfo.tendonForce;
                err_p  = F_t_p - F_m_AT_p;
    
                J       = (err_p - err) / delta;
                if abs(J) < 1e-14
                    error('Cannot update: stiffness is close to zero.');
                end
    
                L_m_AT  = L_m_AT - (err / J);
                if (L_m_AT < 1e-6) || (L_m_AT > L_mt)
                    error('Newton solver failed: Solution %.4e is out of bounds [%.4e, %.4e].', ...
                        L_m_AT, low_L_m_AT, high_L_m_AT);
                end
            end

            muscleState = [V_m_AT; L_m_AT]; % these two should be found
            mtInfo_Eq = calcMillard2012DampedEquilibriumMuscleInfo( ...
                a, pathState, muscleState, ...
                obj.MuscleArch, obj.NormMuscleCurves, obj.ModelConfig);
            F_t    = mtInfo_Eq.muscleDynamicsInfo.tendonForce;

            L_m_height = obj.L_mo * sin(obj.AlphaOpt);
            L_m        = sqrt(L_m_AT^2 + L_m_height^2);
            Alpha      = obj.calc_pennation_angle(L_m);

            obj.a    = a;
            obj.L_m  = L_m;
            obj.V_m  = (L_m - obj.L_m)/dt;
            % obj.L_mt = L_mt; % it is set at the beginning and not modified
            % obj.V_mt = V_mt; % it is set at the beginning and not modified
            obj.L_mn = L_m/obj.L_mo;
            obj.F_t  = F_t;
            obj.F_m  = F_t / cos(Alpha);
        end

        function tau = get_Tau(obj, u, a)
            p = obj.MMM;

            if (u > a)
                tau = p.Tau_a * (0.5 + 1.5 * a);
            else
                tau = p.Tau_d / (0.5 + 1.5 * a);
            end
        end         
    end
    
    methods (Access = private)
        function init_parameters(obj)
            % Initialize using the provided init_MMM logic
            [obj.MuscleArch, obj.NormMuscleCurves, obj.ModelConfig, obj.MMM] = ...
                obj.init_MMM();
        end

        function [muscleArch, normMuscleCurves, modelConfig, MMM] ...
                = init_MMM(obj)

            muscleName = obj.MuscleName;
            muscleAbbr = obj.MuscleName;
            F_mo       = obj.F_mo;
            L_mo       = obj.L_mo;
            L_ts       = obj.L_ts;
            alphaOpt   = obj.AlphaOpt;
            V_mmax_norm = 10;

            % Default Parameters & Constants
            maximumPennationAngle = 89 * (pi/180); 

            % If this is greater than 0 this value will be used to make 
            % the tendon-force-length curve. Otherwise the default of 0.049 is taken.
            tendonStrainAtOneNormForce = 0.033; 
            flag_plotNormMuscleCurves = 0;
            flag_updateCurves = 1;

            % Normalized Muscle Curves Creation (Millard Model)
            normMuscleCurves = createDefaultNormalizedMuscleCurves(muscleName,...
                                                    tendonStrainAtOneNormForce,...
                                                    flag_updateCurves,...
                                                    flag_plotNormMuscleCurves);

            normMuscleCurves.tendonInverseCurve = ...
                createInverseBezierCurve(normMuscleCurves.tendonForceLengthCurve);

            % Muscle Architecture Setup (muscleArch)
            muscleArch = struct();
            muscleArch.name = muscleName;
            muscleArch.abbr = muscleAbbr;
            muscleArch.fiso = F_mo;
            muscleArch.optimalFiberLength = L_mo;
            muscleArch.maximumNormalizedFiberVelocity = V_mmax_norm;
            muscleArch.pennationAngle = alphaOpt;
            muscleArch.tendonSlackLength = L_ts;

            % Kinematics Limit Calculation
            minimumActiveFiberNormalizedLength = normMuscleCurves.activeForceLengthCurve.xEnd(1);

            minFiberKinematics = calcFixedWidthPennatedFiberMinimumLength(...
                        minimumActiveFiberNormalizedLength,...
                        maximumPennationAngle,...
                        muscleArch.optimalFiberLength,...
                        muscleArch.pennationAngle);
            
            muscleArch.minimumFiberLength = ...
                       minFiberKinematics.minimumFiberLength;
            muscleArch.minimumFiberLengthAlongTendon =...
                       minFiberKinematics.minimumFiberLengthAlongTendon;
            muscleArch.pennationAngleAtMinumumFiberLength = ...
                       minFiberKinematics.pennationAngleAtMinimumFiberLength;

            % Model Configuration (modelConfig)
            modelConfig = struct();
            modelConfig.useElasticTendon    = 1;
            modelConfig.useFiberDamping     = 0;  
            modelConfig.damping             = 0.1;
            modelConfig.minActivation       = 1.49012e-08;
            modelConfig.iterMax             = 10000;
            modelConfig.tol                 = 1e-10;

            % Additional Constants (MMM)
            MMM = struct();
            MMM.Tau_a = 0.01;
            MMM.Tau_d = 0.04;
        end
    end
end