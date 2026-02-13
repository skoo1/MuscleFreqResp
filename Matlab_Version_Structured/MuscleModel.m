% By Minseung Kim and Seungbum Koo
% KAIST, Daejeon, South Korea
% February 8, 2026

classdef (Abstract) MuscleModel < handle
    % MUSCLEMODEL Abstract base class for muscle models
    %   Defines the interface for muscle properties and dynamics updates.

    properties
        MuscleName
        MuscleType  % Millard Classic/DEq/Rigid
        F_mo        % Max Isometric Force (N)
        L_mo        % Optimal Fiber Length (m)
        L_ts        % Tendon Slack Length (m)
        AlphaOpt    % Pennation Angle at L_mo (rad)
        Mass        % Muscle Mass (kg)
        Damping     % Viscous Damping (Ns/m)
        V_m_max     % Maximum Muscle Velocity (m/s)

        % Internal States
        a           % Activation (0-1)
        L_m         % Fiber Length (m)
        V_m         % Fiber Velocity (m/s)
        L_mt        % Muscle-Tendon Length (m)
        V_mt        % Muscle-Tendon velocity (m/s)
        L_mn        % Fiber Length normalized
        F_t         % Tendon Force (N)
        F_m         % Fiber Force (N)
    end

    methods
        function obj = MuscleModel(name, F_mo, L_mo, L_ts, alphaOpt, mass_not_used, damping_not_used, type_str)
            obj.MuscleName = name;
            obj.MuscleType = type_str;
            obj.F_mo = F_mo;
            obj.L_mo = L_mo;
            obj.L_ts = L_ts;
            obj.AlphaOpt = alphaOpt;
            obj.Mass = mass_not_used; % Not used in this study
            obj.Damping = damping_not_used; % Not used in this study
            obj.V_m_max = 10 * L_mo; % according to literature

            % Default Initialization
            obj.a = 0;
            obj.L_mt = L_mo * cos(alphaOpt) + L_ts;
            obj.L_m  = L_mo * cos(alphaOpt); % Default projected length
        end

        % --- Abstract Methods (Must be implemented by subclasses) ---

        % Initialize internal states (find equilibrium)
        [L_mt, F_equil, Alpha] = initialize_static_given_Lm(obj, a0, L_m0);
        [L_m0, F_t0, Alpha0]   = initialize_static_given_Lmt(obj, a0, L_mt0);
        
        % Calculate Forces given inputs (Dynamic Simulation)
        % Returns tendon force (F_t) to be used in equations of motion
        F_tendon = updateDynamics(obj, dt, u, F_ext_equil, mass_ext, damping_ext);
        F_tendon = updateDynamicsQuasiStatic(obj, dt, u, L_mt, V_mt);

        % Get muscle specific activation/deactivation time constant
        tau = get_Tau(obj, u, a);
    end

    methods (Access = protected)
        % Common helper
        function da_dt = getActivationRate(obj, u, a)
            % da/dt = (u - a) / tau
            tau = obj.get_Tau(u, a);
            da_dt = (u - a) / tau;
        end

        function alpha = calc_pennation_angle(obj, L_m)
            L_m_height = obj.L_mo * sin(obj.AlphaOpt);
            if L_m < L_m_height
                alpha = 89 * (pi/180);
            else
                alpha = asin(L_m_height / L_m);
            end
        end

    end
end