% By Minseung Kim and Seungbum Koo
% January 26, 2026

clear; clc;

% Sim configuration
% L_mtn_intput = 1.05; % for passive test
L_mn_input     = 1.0;
V_m_input      = 0.0;
U_input        = 0.5;
Amp_input      = 0.01;
Mass_input     = 300000000.0;
Damping_imput  = 0.0;
SimTime_input  = 120;
SimDt_input    = 0.001;
FreqLow_input  = 0.1;
FreqHigh_input = 100;
NumFreqSamples = 100;

% Muscle properties
muscleName = 'Soleus';
F_mo       = 3549.0;
L_mo       = 0.05;
L_ts       = 0.25;
alphaOpt   = 0.4363;
V_mmax     = 10 * L_mo;

% External dynamics parameters
mass = Mass_input;
damp = Damping_imput;

% Operating point
L_mn_target = L_mn_input;
V_m_target  = V_m_input;
u0 = U_input;

% Load Thelen-Muscle model parameters
TMM = init_TMM();

% Muscle parameter initialization
L_mn0   = L_mn_target;
L_m0    = L_mn0*L_mo;
L_mt0   = L_m0 * cos(alphaOpt) + L_ts;
V_mt0   = V_m_target;
V_m0    = V_m_target;
L_m_hight  = L_mo * sin(alphaOpt);

% Calculate the equilibrium external force at initial state
f_l0    = active_muscle_force_length_multiplier(L_mn0, TMM);
F_pn0   = passive_muscle_force_normalized(L_mn0, TMM);
f_v0    = 1;
F_ext_equil = F_mo * (u0 * f_l0 *f_v0 + F_pn0) * cos(alphaOpt); 

% Dynamics simulation parameters
totalTime  = SimTime_input;
dt         = SimDt_input;
time_array = 0 : dt : (totalTime - dt);
time_len   = length(time_array);
steps      = NumFreqSamples;

% Frequency sampling parameters
freqlb = FreqLow_input;
freqhb = FreqHigh_input;
frequencies  = logspace(log10(freqlb), log10(freqhb), steps);
freq_len     = length(frequencies);          

% Storage initialization
visual_result = cell(1, 4);
pha          = zeros(1, length(frequencies));
sin_f        = zeros(1, length(frequencies));
mag          = zeros(1, length(frequencies));
exc          = zeros(1, length(frequencies));

parfor k = 1 : length(frequencies)
    tic;
    freq = frequencies(k);

    F_m_AT    = zeros(1, time_len);

    amp = Amp_input;     % amplitude of the sine shape excitation signal
    u    = sinwave(freq, time_array, u0, u0, amp);
    u(u < 0) = 0; 
    u(u > 1) = 1;

    a    = ones(1, length(u));
    a(1) = u(1);

    V_mt      = V_mt0;
    L_mt      = L_mt0;
    L_mn      = L_mn0;

    for i = 2 : time_len
        % da_dt  = (u(i-1) - a(i-1)) / get_Tau(a(i-1), u(i-1), TMM); % one-step (0.001 ms) lag
        da_dt  = (u(i) - a(i-1)) / get_Tau(a(i-1), u(i), TMM);
        a(i) = a(i-1) + dt * da_dt;

        % parfor temporary variable initialization
        L_mnc = 0; F_tn = 0; F_pn = 0; F_an = 0; Alpha = 0;

        % -----------------------------------------------------------------
        % Newton-Raphson Solver Start
        % -----------------------------------------------------------------
        L_mnc      = L_mn;  % Initial guess: value from previous step
        max_iter   = 50;
        iter_count = 0;
        error1     = 1.0;
        delta      = 1e-7;  % Perturbation step size
        tol        = 1e-8;  % Tolerance

        % Calculate the force measured at the sensor or F_m_AT(i)
        % as a response to the input signal u(i)
        while ( abs(real(error1)) > tol && iter_count < max_iter)
            iter_count = iter_count + 1;

            % --- 1. Calculate Error at Current State (L_mnc) ---
            L_m    = L_mnc * L_mo;
            Alpha  = asin(L_m_hight/L_m);
            L_t    = L_mt - (L_mnc * L_mo) * cos(Alpha);
            V_mn   = ((L_mnc - L_mn) * L_mo / dt) / V_mmax;
            eps_t  = L_t / L_ts - 1;
            F_tn   = tendon_force_normalized(eps_t, TMM);
            f_l    = active_muscle_force_length_multiplier(L_mnc, TMM);
            f_v    = active_muscle_force_velocity_multiplier(V_mn, a(i), TMM);
            F_an   = a(i) * f_l * f_v;
            F_pn   = passive_muscle_force_normalized(L_mnc, TMM);

            error1 = F_tn - (F_pn + F_an) * cos(Alpha);

            % --- 2. Calculate Jacobian using Perturbation (L_mnc + delta) ---
            L_mnc_p = L_mnc + delta;
            L_m_p   = L_mnc_p * L_mo;
            Alpha_p = asin(L_m_hight/L_m_p);
            L_t_p   = L_mt - (L_mnc_p * L_mo) * cos(Alpha_p);
            V_mn_p  = ((L_mnc_p - L_mn) * L_mo / dt) / V_mmax;    % Velocity changes with L_mnc
            eps_t_p = L_t_p / L_ts - 1;
            F_tn_p  = tendon_force_normalized(eps_t_p, TMM);
            f_l_p   = active_muscle_force_length_multiplier(L_mnc_p, TMM);
            f_v_p   = active_muscle_force_normalized(V_mn_p, TMM);
            F_an_p  = a(i) * f_l_p * f_v_p;
            F_pn_p  = passive_muscle_force_normalized(L_mnc_p, TMM);

            error2  = F_tn_p - (F_pn_p + F_an_p) * cos(Alpha_p);

            % --- 3. Update State ---
            J = (error2 - error1) / delta; % Numerical Gradient

            if abs(J) < 1e-14
                J = 1e-14; % Prevent divide by zero
            end

            step = error1 / J;
            L_mnc = L_mnc - step;

            % --- 4. Bounds Check ---
            if L_mnc < 0.01 
                L_mnc = 0.01;
            elseif L_mnc > 1.0
                L_mnc = 1.0;
            end
        end

        L_mn    = L_mnc;
        A_mt    = (F_ext_equil - F_tn * F_mo - V_mt * damp) / mass;
        V_mt    = V_mt + A_mt * dt;
        L_mt    = L_mt + V_mt * dt;

        F_m_AT(i)   = F_mo * (F_an + F_pn) * cos(Alpha);
    end

    valid_range     = round(time_len * 0.1) : time_len;
    sig_in          = u(valid_range);
    sig_out         = F_m_AT(valid_range);

    U_fft           = fft(sig_in);
    F_fft           = fft(sig_out);

    Fs              = 1 / dt;
    n_x             = length(sig_in);
    f               = (0 : n_x - 1) * (Fs / n_x);

    [~, idx]        = min(abs(f - freq));
    sin_f(k)        = freq;
    mag  (k)        = abs(F_fft(idx)) / abs(U_fft(idx));
    pha  (k)        = angle(F_fft(idx)) - angle(U_fft(idx));
    exc  (k)        = abs(U_fft(idx));

    fprintf('%d, %f, %f, %f, %.1f\n', k, freq, mag(k), pha(k), toc);
end

pha = unwrap(pha) * (180 / pi);

visual_result{1} = sin_f;
visual_result{2} = mag;
visual_result{3} = pha;
visual_result{4} = exc;

% Saving the simulation data into .mat file
% Need modify filename to adjust the simulation setting
% sol/gast (muscle type) | YB/OB (aging) | wod/d (damping)
% Format:: TMM_results_KMS_Lmn0_uo_sol_YB_wod.mat %

if exist('mass','var')
    if abs(mass - 3e8) < 1e-6     % 300000000 kg → isometric
        mass_label = 'isometric';
    else
        mass_label = [num2str(mass) 'kg'];
    end
else
    mass_label = 'unknownMass';
end

saveFolder = '.\TMM_results';
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end

fileName = [num2str(L_mn_target) '_' num2str(u0) '_' muscleName '_YB_wod_aF_n_' mass_label '.mat'];
fullPath = fullfile(saveFolder, fileName);
savingdata_freq = visual_result;
save(fullPath, 'savingdata_freq');

% Bode plots
figure();
sg              = sgtitle(sprintf('Bode plot w/ u_{0} = %.2f', u0), 'FontSize', 18);
cm              = hsv(1);

mag_values = 20 * log10(visual_result{2});

subplot(2, 1, 1);
semilogx(visual_result{1}, mag_values, 'o', 'MarkerSize', 1, ...
    'Color', cm, 'DisplayName', sprintf('L_{mn} = %.2f, V_{m} = %.2f', L_mn_target, V_m_target));
ylabel("Magnitude (dB)", 'FontSize', 13);
grid on;
xlim([0.1 100]);
ylim([10 100]);

subplot(2, 1, 2);
semilogx(visual_result{1}, visual_result{3}, 'o', 'MarkerSize', 1, ...
    'Color', cm, 'DisplayName', sprintf('L_{mn} = %.2f, V_{m} = %.2f', L_mn_target, V_m_target));
xlabel("Frequency (Hz)", 'FontSize', 15);
ylabel("Phase (deg)", 'FontSize', 13);
grid on;
xlim([0.1 100]);
ylim([-180 10]);


%% Functions
function input = sinwave(frq, t, a0, a, b)
    if b == 0
        phase_shift = 0;
    else
        val = (a0 - a) / b;
        % asin range protection
        if val > 1, val = 1; end
        if val < -1, val = -1; end
        phase_shift = asin(val);
    end

    input = a + sin(2 * pi * frq * t + phase_shift) * b;
end


function Tau = get_Tau(a, u, param)   
    if (u > a)
        Tau = param.Tau_a * (0.5 + 1.5 * a);
    else
        Tau = param.Tau_d / (0.5 + 1.5 * a);
    end
end


function eps_t = get_tendon_strain(F_tn, param)
    % tendon strain with exponential behavior
    eps_t1 = param.EPSttoe / param.K_toe * log(F_tn / param.F_ttoe * (exp(param.K_toe) - 1) + 1);
    % tendon strain with linear behavior
    eps_t2 = (F_tn - param.F_ttoe) / param.K_lin + param.EPSttoe;
    
    if (F_tn <= param.F_ttoe)
        eps_t = eps_t1;
    else
        eps_t = eps_t2;
    end
end


function F_pn = passive_muscle_force_normalized(L_mn, param)                                
    if (L_mn > 1.0)
        F_pn = (exp((param.K_pe * (L_mn - 1)) / param.EPSmo) - 1) / (exp(param.K_pe) - 1);
    else
        F_pn = 0;
    end
end


function f_v = active_muscle_force_velocity_multiplier(V_mn, a, param)
    v1 = V_mn / (0.25 + 0.75 * a);
   
    if (V_mn <= 0)
        f_v = (1 + v1) / (1 - v1 / param.A_f);
    else
        % When V is positive
        someth = v1 * (2 + 2 / param.A_f) / (param.F_mnlen - 1);
        f_v = (1 + param.F_mnlen * someth) / (1 + someth);
    end
end


function F_tn = tendon_force_normalized(eps_t, param)
    if (eps_t <= param.EPSttoe)
        F_tn = param.F_ttoe / (exp(param.K_toe) - 1) * (exp(param.K_toe*eps_t / param.EPSttoe) - 1);
    else
        F_tn = param.K_lin * (eps_t - param.EPSttoe) + param.F_ttoe;
    end 
end


function f_l = active_muscle_force_length_multiplier(L_mn, param)                                 
    f_l = exp(-((L_mn - 1).^2) / param.Gamma);
end


function TMM = init_TMM()
    TMM.A_f     = 0.3;                     % Force-velocity shape factor
    TMM.F_mnlen = 1.4;                     % Ratio of maximum lengthening muscle fiber force to isometric force
    TMM.Gamma   = 0.45;                    % Active force-length exponential shape factor
    TMM.K_pe    = 4.0;                     % Passive force-length exponential shape factor
    TMM.EPSmo   = 0.6;                     % Passive muscle strain due to maximum isometric force
    TMM.EPSto   = 0.033;                   % Tendon strain due to maximum isometric force
    TMM.EPSttoe = 0.609 * TMM.EPSto;       % Tendon strain above which the tendon exhibits linear behavior
    TMM.F_ttoe  = 0.33;                    % Normalized tendon force at the transition from nonlinear to linear behavior
    TMM.K_toe   = 3;                       % Tendon exponential shape factor
    TMM.K_lin   = 1.712 / TMM.EPSto;       % Tendon linear scale factor
    TMM.Tau_a   = 0.01;                    % Muscle fiber activation time constant 
    TMM.Tau_d   = 0.04;                    % Muscle fiber deactivation time constant
end