% By Minseung Kim and Seungbum Koo
% KAIST, Daejeon, South Korea
% February 8, 2026

classdef ThelenMuscle < MuscleModel
    properties
        TMM % Struct containing Thelen model specific parameters
    end

    methods
        function obj = ThelenMuscle(name, F_mo, L_mo, L_ts, alphaOpt, mass, damping, ~)
            % Call superclass constructor
            obj@MuscleModel(name, F_mo, L_mo, L_ts, alphaOpt, mass, damping);
            obj.init_parameters();
        end

        function [L_mt, F_equil, Alpha] = initialize_static_given_Lm(obj, a, L_m)
            Alpha = obj.calc_pennation_angle(L_m);

            % Static Muscle Forces (V_m = 0 -> f_v = 1)
            L_mn = L_m/obj.L_mo;
            f_l  = obj.active_muscle_force_length_multiplier(L_mn);
            f_v  = 1.0; % f_v is 1.0 at static
            F_pn = obj.passive_muscle_force_normalized(L_m);
            F_an = a * f_l * f_v;

            % Projected Muscle Force (Normalized)
            F_mn_AT = (F_an + F_pn) * cos(Alpha);

            % Equilibrium: F_tendon MUST equal F_muscle_proj
            F_tn = F_mn_AT;
            F_equil = F_tn * obj.F_mo;

            % Inverse Tendon Model: Find strain (eps_t) from Force (F_tn)
            eps_t = obj.get_tendon_strain_from_force(F_tn);

            % Calculate L_mt
            L_t = obj.L_ts * (1 + eps_t);
            L_mt = L_m * cos(Alpha) + L_t;

            % update muscle state
            obj.a    = a;
            obj.L_m  = L_m;
            obj.V_m  = 0; % static condition
            obj.L_mt = L_mt;
            obj.V_mt = 0; % static condition
            obj.L_mn = L_mn;
            obj.F_t  = F_tn * obj.F_mo;
            obj.F_m  = obj.F_t / cos(Alpha);
        end

        function [L_m, F_t, Alpha] = initialize_static_given_Lmt(obj, a, L_mt)
            low_L_mn      = 0.0;
            high_L_mn     = 2.0;
            err        = 1000;

            while (abs(real(err)) > 1e-11)
                L_mn   = 1/2 * (low_L_mn + high_L_mn);
                L_m    = L_mn * obj.L_mo;
                Alpha  = obj.calc_pennation_angle(L_m);
                L_t    = L_mt - L_m * cos(Alpha);
                eps_t  = L_t / obj.L_ts - 1;

                F_tn   = obj.get_tendon_force_normalized_from_strain(eps_t);
                F_pn   = obj.passive_muscle_force_normalized(L_mn);
                % Assume a = 0, thus muscle active force F_an = zero 
                F_an   = 0;

                err = F_tn - (F_an + F_pn) * cos(Alpha);

                if err < 0
                    high_L_mn = L_mn;
                else
                    low_L_mn  = L_mn;
                end
            end

            F_t = F_tn * obj.F_mo;
            obj.a    = a;
            obj.L_m  = L_m;
            obj.V_m  = 0; % initially, static condition
            obj.L_mt = L_mt;
            obj.V_mt = 0; % initially, static condition
            obj.L_mn = L_mn;
            obj.F_t  = F_t;
            obj.F_m  = obj.F_t / cos(Alpha);
        end

        function F_t = updateDynamics(obj, dt, u, F_ext_equil, mass_ext, damping_ext)
            % --- External interaction: Environment involves F_ext_equil, mass_ext, and damping_ext ---
            % Calculate muscle-tendon force (F_t) or projected fiber force (F_m_AT) 
            % based on the control input (u).
            % Internal states, such as fiber length (L_m) and tendon length (L_t), 
            % are updated within the model.

            % Update Activation
            % da_dt  = (u - obj.a) / obj.get_Tau(u, obj.a);
            da_dt = obj.getActivationRate(u, obj.a);
            a = obj.a + dt * da_dt;

            % Solve for Fiber Length (Newton-Raphson)
            % Solver finds L_mn that satisfies F_t = (F_p + F_a) * cos(alpha)
            % Uses Implicit Euler for velocity: V_mn = (L_mn - L_mn_prev)/dt
            L_mnc     = obj.L_mn; % Initial guess
            L_mn_prev = obj.L_mn;
            L_mt      = obj.L_mt;

            max_iter  = 50;
            iter      = 0;
            err     = 1.0;
            tol       = 1e-8;
            delta     = 1e-7;

            while (abs(err) > tol && iter < max_iter)
                iter = iter + 1;

                % Eval error at current guess
                err = obj.calcForceError(L_mnc, L_mn_prev, L_mt, a, dt);

                % Eval error at perturbed guess
                err_p = obj.calcForceError(L_mnc + delta, L_mn_prev, L_mt, a, dt);

                % Jacobian
                J = (err_p - err) / delta;
                if abs(J) < 1e-14, J = 1e-14; end

                L_mnc = L_mnc - (err / J);

                % Clamp
                if L_mnc < 0.01, L_mnc = 0.01; elseif L_mnc > 2.0, L_mnc = 2.0; end
            end

            % Update Internal States: Finalize Muscle-Tendon force at equilibrium
            L_m  = L_mnc * obj.L_mo;
            alpha = obj.calc_pennation_angle(L_m);
            L_t   = L_mt - L_m * cos(alpha);
            eps_t = L_t / obj.L_ts - 1;
            F_tn  = obj.get_tendon_force_normalized_from_strain(eps_t);
            F_t   = F_tn * obj.F_mo;
            F_m   = F_t / cos(alpha);

            % Update External Dynamics: Interaction with environment (Mass-Damper)
            % Compute acceleration based on net force (External - Muscle Tension)
            A_mt = (F_ext_equil - F_t - damping_ext * obj.V_mt) / mass_ext;
            V_mt = obj.V_mt + A_mt * dt;
            L_mt = L_mt + V_mt * dt;

            obj.a    = a;
            obj.V_m  = (L_m - obj.L_m) / dt;
            obj.L_m  = L_m;
            obj.V_mt = V_mt;
            obj.L_mt = L_mt;
            obj.L_mn = L_mnc;
            obj.F_t  = F_t;
            obj.F_m  = F_m;
        end

        function F_t = updateDynamicsQuasiStatic(obj, dt, u, L_mt, V_mt)
            obj.L_mt = L_mt;
            obj.V_mt = V_mt; % quasi-static condition

            % activation dynamics
            da_dt = obj.getActivationRate(u, obj.a);
            a = obj.a + dt * da_dt;

            % initial guess and value
            L_mnc     = obj.L_mn;
            L_mn_prev = obj.L_mn;
            % L_mt      = L_mt; % preset

            max_iter  = 50;
            iter      = 0;
            err     = 1.0;
            tol       = 1e-8;
            delta     = 1e-7;

            while (abs(err) > tol && iter < max_iter)
                iter = iter + 1;

                % Eval error at current guess
                err = obj.calcForceError(L_mnc, L_mn_prev, L_mt, a, dt);

                % Eval error at perturbed guess
                err_p = obj.calcForceError(L_mnc + delta, L_mn_prev, L_mt, a, dt);

                % Jacobian
                J = (err_p - err) / delta;
                if abs(J) < 1e-14, J = 1e-14; end

                L_mnc = L_mnc - (err / J);

                % Clamp
                if L_mnc < 0.01, L_mnc = 0.01; elseif L_mnc > 2.0, L_mnc = 2.0; end
            end

            if iter >= max_iter && abs(err) > tol
                error('Newton solver failed to converge within %d iterations (Error: %.2e).', max_iter, err);
            end

            % Update Internal States: Finalize Muscle-Tendon force at equilibrium
            L_m  = L_mnc * obj.L_mo;
            alpha = obj.calc_pennation_angle(L_m);
            L_t   = L_mt - L_m * cos(alpha);
            eps_t = L_t / obj.L_ts - 1;
            F_tn  = obj.get_tendon_force_normalized_from_strain(eps_t);
            F_t   = F_tn * obj.F_mo;
            F_m   = F_t / cos(alpha);

            obj.a    = a;
            obj.V_m  = (L_m - obj.L_m) / dt;
            obj.L_m  = L_m;
            % obj.V_mt = V_mt; % it is set at the beginning and not modified
            % obj.L_mt = L_mt; % it is set at the beginning and not modified
            obj.L_mn = L_mnc;
            obj.F_t  = F_t;
            obj.F_m  = F_m;
        end

        function tau = get_Tau(obj, u, a)
            p = obj.TMM;

            if (u > a)
                tau = p.Tau_a * (0.5 + 1.5 * a);
            else
                tau = p.Tau_d / (0.5 + 1.5 * a);
            end
        end
    end

    methods (Access = private)
        function init_parameters(obj)
            % Default Thelen Parameters
            obj.TMM.A_f     = 0.3;
            obj.TMM.F_mnlen = 1.4;
            obj.TMM.Gamma   = 0.45;
            obj.TMM.K_pe    = 4.0;
            obj.TMM.EPSmo   = 0.6;
            obj.TMM.EPSto   = 0.033;
            obj.TMM.EPSttoe = 0.609 * obj.TMM.EPSto;
            obj.TMM.F_ttoe  = 0.33;
            obj.TMM.K_toe   = 3;
            obj.TMM.K_lin   = 1.712 / obj.TMM.EPSto;
            obj.TMM.Tau_a   = 0.01;
            obj.TMM.Tau_d   = 0.04;
        end

        function err = calcForceError(obj, L_mn, L_mn_prev, L_mt, a, dt)
            L_m = L_mn * obj.L_mo;

            % Dynamic Pennation
            alpha = obj.calc_pennation_angle(L_m);
            L_t   = L_mt - L_m * cos(alpha);
            eps_t = L_t / obj.L_ts - 1;
            F_tn  = obj.get_tendon_force_normalized_from_strain(eps_t);
            V_mn = ((L_mn - L_mn_prev) * obj.L_mo / dt) / (obj.V_m_max);
            f_l  = obj.active_muscle_force_length_multiplier(L_mn);
            f_v  = obj.active_muscle_force_velocity_multiplier(V_mn, a);
            F_pn = obj.passive_muscle_force_normalized(L_mn);
            F_an = a * f_l * f_v;

            err = F_tn - (F_pn + F_an) * cos(alpha);
        end

        function F_tn = get_tendon_force_normalized_from_strain(obj, eps_t)
            p = obj.TMM;
            if eps_t <= p.EPSttoe
                F_tn = p.F_ttoe / (exp(p.K_toe) - 1) * (exp(p.K_toe * eps_t / p.EPSttoe) - 1);
            else
                F_tn = p.K_lin * (eps_t - p.EPSttoe) + p.F_ttoe;
            end
            if F_tn < 0, F_tn = 0; end
        end

        function eps_t = get_tendon_strain_from_force(obj, F_tn)
            p = obj.TMM;
            % tendon strain with exponential behavior
            eps_t1 = p.EPSttoe / p.K_toe * log(F_tn / p.F_ttoe * (exp(p.K_toe) - 1) + 1);
            % tendon strain with linear behavior
            eps_t2 = (F_tn - p.F_ttoe) / p.K_lin + p.EPSttoe;

            if (F_tn <= p.F_ttoe)
                eps_t = eps_t1;
            else
                eps_t = eps_t2;
            end
        end

        function f_l = active_muscle_force_length_multiplier(obj, L_mn)
            p = obj.TMM;
            f_l = exp(-((L_mn - 1).^2) / p.Gamma);
        end

        function f_v = active_muscle_force_velocity_multiplier(obj, V_mn, a)
            p = obj.TMM;
            v1 = V_mn / (0.25 + 0.75 * a);

            if (V_mn <= 0)
                f_v = (1 + v1) / (1 - v1 / p.A_f);
            else
                % When V is positive
                someth = v1 * (2 + 2 / p.A_f) / (p.F_mnlen - 1);
                f_v = (1 + p.F_mnlen * someth) / (1 + someth);
            end
        end

        function F_pn = passive_muscle_force_normalized(obj, L_mn)
            p = obj.TMM;
            if (L_mn > 1.0)
                F_pn = (exp((p.K_pe * (L_mn - 1)) / p.EPSmo) - 1) / (exp(p.K_pe) - 1);
            else
                F_pn = 0;
            end
        end

    end
end