% By Minseung Kim and Seungbum Koo
% February 5, 2026

clear; clc;

% Sim configuration
L_mtn_input    = 1.05; % for passive test
% L_mn_input   = 1.0;
% V_m_input    = 0.0;
U_input        = 0.0;
Amp_input      = 0.005;
Mass_input     = 0.0;
Damping_input  = 0.0;
SimTime_input  = 120;
SimDt_input    = 0.001;
FreqLow_input  = 0.1;
FreqHigh_input = 100;
NumFreqSamples = 100;

% Muscle properties
muscleName = 'Soleus';
F_mo       = 3549.0;    % Isometric muscle force      (N)
L_mo       = 0.05;      % Optimal muscle fiber length (m)
L_ts       = 0.25;      % Tendon slack length         (m)
alphaOpt   = 0.4363;    % Pennation angle at optimal muscle length (rad)
V_mmax     = 10 * L_mo;

% External dynamics parameters
mass_ext      = Mass_input;
damping_ext   = Damping_input;

% Operating point
L_mtn_target  = L_mtn_input;

% Load Thelen-Muscle model parameters
TMM           = init_TMM();

L_mt0         = L_mo * cos(alphaOpt) + L_ts;
L_m_hight     = L_mo * sin(alphaOpt);

low_L_mn      = 0.0;
high_L_mn     = 2.0;
error1        = 1000;

while (abs(real(error1)) > 1e-11)
    L_mn   = 1/2 * (low_L_mn + high_L_mn);
    L_m    = L_mn * L_mo;
    Alpha  = asin(L_m_hight/L_m);
    L_t    = L_mtn_target * L_mt0 - L_m * cos(Alpha);
    eps_t  = L_t / L_ts - 1;

    F_tn   = tendon_force_normalized(eps_t, TMM);
    F_pn   = passive_muscle_force_normalized(L_mn, TMM);

    error1 = F_tn - F_pn * cos(Alpha);

    if error1 < 0
        high_L_mn = L_mn;
    else
        low_L_mn  = L_mn;
    end
end

% Operating point
L_mn0 = L_mn;
u0    = U_input;

totalTime   = SimTime_input;
dt          = SimDt_input;
time_array  = 0 : dt : (totalTime - dt);
time_len    = length(time_array);
steps       = NumFreqSamples;

freqlb      = FreqLow_input;
freqhb      = FreqHigh_input;
frequencies = logspace(log10(freqlb), log10(freqhb), steps);          
freq_len    = length(frequencies);

% Storage initialization
visual_result = cell(1, 4);
pha           = zeros(1, length(frequencies));
sin_f         = zeros(1, length(frequencies));
mag           = zeros(1, length(frequencies));
exc           = zeros(1, length(frequencies));

parfor k = 1 : length(frequencies)
    tic;
    freq          = frequencies(k);

    L_mn          = L_mn0;
    F_m_AT        = zeros(1, time_len);

    amp           = Amp_input * L_mt0 * L_mtn_target;
    L_mt_preset   = L_mt0 * L_mtn_target + amp * sin(2 * pi * freq * time_array);
    V_mt          = amp * 2 * pi * freq * cos(2 * pi * freq * time_array);
    A_mt          = -amp * (2 * pi * freq)^2 * sin(2 * pi * freq * time_array);

    for i = 2 : time_len
        L_mt        = L_mt_preset(i);

        % initial value
        L_mnc       = L_mn;

        % Solver parameters
        max_iter    = 50;
        iter_count  = 0;
        error1      = 1.0;
        delta       = 1e-7;
        tol         = 1e-8;

        % parfor temporary variable initialization
        Alpha       = 0;
        F_pn        = 0;
        F_an        = 0;
        a           = u0;
        while ( abs(real(error1)) > tol && iter_count < max_iter)
            iter_count = iter_count + 1;

            L_m   = L_mo * L_mnc;
            Alpha = asin(L_m_hight / L_m);
            L_t   = L_mt - (L_mnc * L_mo) * cos(Alpha);
            V_mn  = ((L_mnc - L_mn) * L_mo / dt) / V_mmax;
            eps_t = L_t / L_ts - 1;
            F_tn  = tendon_force_normalized(eps_t, TMM);
            f_l   = active_muscle_force_length_multiplier(L_mnc, TMM);
            f_v   = active_muscle_force_velocity_multiplier(V_mn, a, TMM);
            F_an  = a * f_l * f_v;
            F_pn  = passive_muscle_force_normalized(L_mnc, TMM);

            error1  = F_tn - (F_an + F_pn) * cos(Alpha);

            L_mnc_p = L_mnc + delta;
            L_m_p   = L_mo * L_mnc_p;
            Alpha_p = asin(L_m_hight / L_m_p);
            L_t_p   = L_mt - (L_mnc_p * L_mo) * cos(Alpha_p);
            V_mn_p  = ((L_mnc_p - L_mn) * L_mo / dt) / V_mmax;
            eps_t_p = L_t_p / L_ts - 1;
            F_tn_p  = tendon_force_normalized(eps_t_p, TMM);
            f_l_p   = active_muscle_force_length_multiplier(L_mnc_p, TMM);
            f_v_p   = active_muscle_force_velocity_multiplier(V_mn_p, a, TMM);
            F_an_p  = a * f_l_p * f_v_p;
            F_pn_p  = passive_muscle_force_normalized(L_mnc_p, TMM);
            
            error2  = F_tn_p - (F_an_p + F_pn_p) * cos(Alpha_p);

            J = (error2 - error1) / delta;
            
            if abs(J) < 1e-14
                J = 1e-14;
            end
            
            step  = error1 / J;
            L_mnc = L_mnc - step;
            
            if L_mnc < 0.01 
                L_mnc = 0.01;
            elseif L_mnc > 2.0
                L_mnc = 2.0;
            end
        end            

        if iter_count >= max_iter
             error('Warning: Newton solver did not converge at step %d\n', i);
        end

        L_mn         = L_mnc;
        
        if (L_mn * L_mo < L_mo)
            error("Fiber length is shorter than Lopt");
        end

        L_m   = L_mo * L_mnc;
        Alpha = asin(L_m_hight / L_m);
        L_t   = L_mt - (L_mnc * L_mo) * cos(Alpha);
        eps_t = L_t / L_ts - 1;
        F_tn  = tendon_force_normalized(eps_t, TMM);

        F_m_AT(i)   = F_mo * F_tn;
    end          

    valid_range     = round(time_len * 0.1) : time_len;
    sig_in          = L_mt_preset(valid_range);
    sig_out         = F_m_AT(valid_range);

    Fs              = 1 / dt;
    n_x             = length(sig_in);
    f               = (0 : n_x - 1) * (Fs / n_x);
    
    X_fft           = fft(sig_in);
    F_fft           = fft(sig_out);

    [~, idx]        = min(abs(f - freq));
    sin_f(k)        = freq;
    mag  (k)        = abs(F_fft(idx)) / abs(X_fft(idx));
    pha  (k)        = angle(F_fft(idx)) - angle(X_fft(idx));

    fprintf('%d, %f, %f, %f, %.1f\n', k, freq, mag(k), pha(k), toc);
end

pha = unwrap(pha) * (180 / pi);

visual_result{1} = sin_f;
visual_result{2} = mag;
visual_result{3} = pha;
visual_result{4} = exc;


%% Extracting simulation data into .mat file
% Need modify filename to adjust the simulation setting
% sol/gast (muscle type) | YB/OB (aging) | wod/d (damping)
% Format:: TMM_results_KMS_Lmn0_uo_sol_YB_wod.mat %

if exist('mass_ext','var')
    if abs(mass_ext - 3e9) < 1e-6     % 3e9 kg → static
        mass_label = 'isometric';
    else
        mass_label = [num2str(mass_ext) 'kg'];
    end
else
    mass_label = 'unknownMass';
end

saveFolder = '.\TMM_results';
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end

fileName        = [num2str(L_mtn_target) '_' num2str(u0) '_' muscleName '_YB_wod_aF_passive.mat'];
fullPath        = fullfile(saveFolder, fileName);
savingdata_freq = visual_result;
save(fullPath, 'savingdata_freq');

% Bode plots
figure();
sg       = sgtitle(sprintf('Bode plot of each conditions w/ u_{0} = %.2f', u0), 'FontSize', 18);
cm       = hsv(1);

mag_values = 20 * log10(visual_result{2});

subplot(2, 1, 1);
semilogx(visual_result{1}, mag_values, 'o', 'MarkerSize', 1, ...
    'Color', cm, 'DisplayName', sprintf('L_{mtn} = %.2f', L_mtn_target));
ylabel("Magnitude (dB)", 'FontSize', 13);
grid on;
xlim([0.1 100]);
ylim([10 100]);

subplot(2, 1, 2);
semilogx(visual_result{1}, visual_result{3}, 'o', 'MarkerSize', 1, ...
    'Color', cm, 'DisplayName', sprintf('L_{mtn} = %.2f', L_mtn_target));
xlabel("Frequency (Hz)", 'FontSize', 15);
ylabel("Phase (deg)", 'FontSize', 13);
grid on;
xlim([0.1 100]);
ylim([-180 10]);


%% Functions
function F_pn = passive_muscle_force_normalized(L_mn, param)                                
    if (L_mn > 1.0)
        F_pn = (exp((param.K_pe * (L_mn - 1)) / param.EPSmo) - 1) / (exp(param.K_pe) - 1);
    else
        F_pn = 0;
    end
end


function f_l = active_muscle_force_length_multiplier(L_mn, param)                                 
    f_l = exp(-((L_mn - 1).^2) / param.Gamma);
end


function f_v = active_muscle_force_velocity_multiplier(V_mn, a, param)
    v1 = V_mn / (0.25 + 0.75 * a);
   
    if (V_mn <= 0)
        f_v = (1 + v1) / (1 - v1 / param.A_f);
    else
        % When V is positive
        someth = v1 * (2 + 2 / param.A_f) / (param.F_mnlen - 1);
        f_v    = (1 + param.F_mnlen * someth) / (1 + someth);
    end
end


function F_tn = tendon_force_normalized(eps_t, param)
    if (eps_t <= param.EPSttoe)
        F_tn = param.F_ttoe / (exp(param.K_toe) - 1) * (exp(param.K_toe*eps_t / param.EPSttoe) - 1);
    else
        F_tn = param.K_lin * (eps_t - param.EPSttoe) + param.F_ttoe;
    end 
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