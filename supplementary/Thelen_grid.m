clear; clc;

V_mmax      = 10;
L_mo        = 0.05;
L_ts        = 0.25;      
F_mo        = 3549.0;   
alpha       = 0;         
activations = [0.5, 1.0];

n_len = 101;        % length grid count
n_vel = 101;        % velocity grid count

norm_len_min = 0.5; 
norm_len_max = 1.5;
norm_vel_min = -1.0; 
norm_vel_max = 1.0;

norm_fiber_length   = linspace(norm_len_min, norm_len_max, n_len);
norm_fiber_velocity = linspace(norm_vel_min, norm_vel_max, n_vel);

assert(numel(norm_fiber_length) == n_len);
assert(numel(norm_fiber_velocity) == n_vel);

N = numel(activations) * n_len * n_vel;

activation_col          = zeros(N,1);
norm_fiber_length_col   = zeros(N,1);
norm_fiber_velocity_col = zeros(N,1);
fiber_velocity_m_per_s  = zeros(N,1);

activeAT_N   = zeros(N,1);
passiveAT_N  = zeros(N,1);
totalAT_N    = zeros(N,1);

activeAT_norm  = zeros(N,1);
passiveAT_norm = zeros(N,1);
totalAT_norm   = zeros(N,1);

Vmax_mps = V_mmax * L_mo;     
cosAlpha = cos(alpha);        

idx = 1;

for ia = 1:numel(activations)
    a = activations(ia);

    for iL = 1:n_len
        Lmn = norm_fiber_length(iL);

        fl  = active_muscle_force_length_multiplier(Lmn);

        fpe = passive_muscle_force_length_multiplier(Lmn);

        for iV = 1:n_vel
            Vn = norm_fiber_velocity(iV);

            factive     = active_force(Vn, fl, a);
            fpassive    = fpe;

            aAT_norm    = factive * cosAlpha;
            pAT_norm    = fpassive * cosAlpha;
            tAT_norm    = aAT_norm + pAT_norm;

            activation_col(idx)          = a;
            norm_fiber_length_col(idx)   = Lmn;
            norm_fiber_velocity_col(idx) = Vn;
            fiber_velocity_m_per_s(idx)  = Vn * 10;

            activeAT_norm(idx)  = aAT_norm;
            passiveAT_norm(idx) = pAT_norm;
            totalAT_norm(idx)   = tAT_norm;

            activeAT_N(idx)  = aAT_norm * F_mo;
            passiveAT_N(idx) = pAT_norm * F_mo;
            totalAT_N(idx)   = tAT_norm * F_mo;

            idx = idx + 1;
        end
    end
end

T = table( ...
    activation_col, ...
    norm_fiber_length_col, ...
    norm_fiber_velocity_col, ...
    fiber_velocity_m_per_s, ...
    activeAT_N, passiveAT_N, totalAT_N, ...
    activeAT_norm, passiveAT_norm, totalAT_norm);

T.Properties.VariableNames = { ...
    'activation', ...
    'norm_fiber_length', ...
    'norm_fiber_velocity', ...
    'fiber_velocity_m_per_s', ...
    'activeAT_N', 'passiveAT_N', 'totalAT_N', ...
    'activeAT_norm', 'passiveAT_norm', 'totalAT_norm'};

out_name = 'Thelen_grid.csv';

writetable(T, out_name);

fprintf('Saved: %s (rows=%d)\n', out_name, height(T));

%% Functions

function F_pe = passive_muscle_force_length_multiplier(L_mn)      
    K_pe  = 5.0;    
    EPSmo = 0.6; 
    if (L_mn > 1.0)
        F_pe = (exp((K_pe * (L_mn - 1)) / EPSmo) - 1) / (exp(K_pe) - 1);
    else
        F_pe = 0;
    end
end

function F_mn = active_force(V_temp, f_l, a)                                                
    A_f     = 0.3;                                                                          
    F_mnlen = 1.4;      

    velocity_norm = V_temp / (0.25 + 0.75 * a);
   
    if (V_temp <= 0)
        % When V is negative
        velocity_factor = (1 + velocity_norm) / (1 - velocity_norm / A_f);                 
    else
        % When V is positive
        someth = velocity_norm * (2 + 2 / A_f) / (F_mnlen - 1);
        velocity_factor = (1 + F_mnlen * someth) / (1 + someth);
    end
    F_mn = a * f_l * velocity_factor;
end

function f_l0 = active_muscle_force_length_multiplier(L_mn)                                 
    Gamma = 0.45;                                                                           
    f_l0 = exp(-((L_mn - 1).^2) / Gamma);
end