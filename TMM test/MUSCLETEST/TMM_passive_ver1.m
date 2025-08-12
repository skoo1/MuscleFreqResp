% By Minseung Kim, 2025-04-10
% Testing passive impedance of the muscle model

clear;
clc;
global F_mnlen K_pe Gamma A_f EPSmo EPSto EPSttoe F_ttoe K_toe K_lin Tau_a Tau_d V_mmax Alpha;

%% Problem Assumptions

% [ Properties of Soleus muscle ]

L_mo    = 0.05;                                                              % Optimal muscle fiber length       (m)     Previous one was 0.05
L_ts    = 0.25;                                                              % Tendon slack length               (m)     Pervious one was 0.25
F_mo    = 3549.0;                                                            % Maximum muscle isometric force    (N)     Previous one was 3549
% F_mo    = 2484.3;
Alpha   = 0.4363;                                                            % Pennationangle                    (rad)
% Alpha   = 0.3050;

% mass    = 3000000.0;                                                       % Mass of the muscle                (kg)
mass    = 30.0;
damp    = 1000.0;
      
% -----------------------------------------------------------

% [ Properties of Gastrocnemius medialis muscle ]

% L_mo    = 0.06;                                                              % Optimal muscle fiber length       (m)     Previous one was 0.05
% L_ts    = 0.39;                                                              % Tendon slack length               (m)     Pervious one was 0.25
% F_mo    = 1558.0;                                                            % Maximum muscle isometric force    (N)     Previous one was 3549
% Alpha   = 0.2967;                                                            % Pennationangle                    (rad)
% mass    = 30.0;                                                              % Mass of the muscle                (kg)
% damp    = 0.0;
% k_temp  = 0.0;                                                               % Stiffness                 

% [ Properties of Tibialis anterior muscle ]

% F_mo    = 905.0;        % (N)
% L_mo    = 0.0898;       % (m)
% Alpha   = 0.0873;       % (rad)
% L_ts    = 0.2043;       % (m)         
% mass    = 30.0;
% damp    = 1000.0;

% [ Properties of Biceps Femoris muscle ]

% F_mo    = 804.0;        % (N)
% L_mo    = 0.1730;       % (m)
% Alpha   = 0.4014;       % (rad)
% L_ts    = 0.0890;       % (m)         
% mass    = 30.0;
% damp    = 1000.0;

%% Constants 

A_f     = 0.3;                                                                % Force-velocity shape factor
% A_f     = 0.25;
F_mnlen = 1.4;                                                                % Ratio of maximum lengthening muscle fiber force to isometric force
% F_mnlen = 1.8;
Gamma   = 0.45;                                                               % Active force-length exponential shape factor
% Gamma   = 0.40;
K_pe    = 4.0;                                                                % Passive force-length exponential shape factor
% K_pe    = 3.5; 

EPSmo   = 0.6;                                                                % Passive muscle strain due to maximum isometric force
% EPSmo   = 0.5;
EPSto   = 0.033;                                                              % Tendon strain due to maximum isometric force
EPSttoe = 0.609 * EPSto;                                                      % Tendon strain above which the tendon exhibits linear behavior
F_ttoe  = 0.33;                                                               % Normalized tendon force at the transition from nonlinear to linear behavior
K_toe   = 3;                                                                  % Tendon exponential shape factor
K_lin   = 1.712 / EPSto;                                                      % Tendon linear scale factor

Tau_a   = 0.01;                                                               % Muscle fiber activation time constant 
Tau_d   = 0.04;                                                               % Muscle fiber deactivation time constant
% Tau_d   = 0.06;
V_mmax  = 10 * L_mo;                                                          % Maximum muscle contraction velocity expressed in optimal fiber lengths per second
% V_mmax  = 8 * L_mo;

%% Generating patch_info.txt

% filename        = 'patch_info.txt';
filename        = 'single_patch.txt';

fileID          = fopen(filename, 'w');
fprintf(fileID, 'L_mn0 V_m0\n');

% L_mn0_values    = [0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4];
% L_mn0_values    = [0.6, 0.8, 1.0, 1.2, 1.4];
L_mn0_values    = [1.0];
V_m0_values     = [0.0];
L_mt0_values    = [1.05];

for file_i = 1 : length(L_mn0_values)
    for file_j = 1 : length(V_m0_values)
        fprintf(fileID, '%.1f %.1f\n', L_mn0_values(file_i), V_m0_values(file_j));
    end
end

fclose(fileID);

if ~exist('patch_result', 'dir')
    mkdir('patch_result');
end

if exist('patch_result', 'dir')
    delete('patch_result/*.mat');  % 디스크에 남은 파일 초기화
end

%% Main Calculation

% patch_info = readtable('patch_info.txt');
patch_info = readtable('single_patch.txt');

index_L  = 1;
index_mt = 1;
index_V  = 1;
pnum     = 1;
 
correct_signal = 0;

L_range = 0.100;                                                                     % Muscle length +- range
V_range = 0.100;                                                                     % Velocity +- range

aim_Lnum = length(L_mn0_values);
aim_Ltnum = length(L_mt0_values);
aim_Vnum = length(V_m0_values);
aim_pnum = height(patch_info);

patch_result  = cell(aim_pnum, 1000, 1);
visual_result = cell(aim_pnum, 1);

while (index_mt <= aim_Ltnum)

    % const_b         = 0.01;
    const_b         = 0;
    cutpoint        = 1;                                                              % cutpoint = 1 means no data cutting

    L_mt0           = L_mo * cos(Alpha) + L_ts;                                                    % muscle-tendon length
    V_mt0           = 0;                                                              % muscle-tendon velocity
        
    mod_V_m0        = V_m0_values(index_V);
    mod_L_mn0       = L_mn0_values(index_L);

    mod_L_m0        = mod_L_mn0 * L_mo;

    mod_L_mt0       = L_mt0_values(index_mt); %%%%%%
   
    f_l0            = active_muscle_force_length_multiplier(mod_L_mn0);                % active force-length scale factor 
    fpe0            = passive_muscle_force_length_multiplier(mod_L_mn0);               % normalized muscle passive force

    u0              = 0.5;
    F_ext           = 0;

    F_n0            = F_ext / F_mo;                                                % normalized external force
    Eps_t0          = tendon_strain_decider2(F_n0);                                    % tendon strain by tendon_force region (linear, exponential)
    L_t0            = (1 + Eps_t0) * (L_ts);                                           % tendon length
    
    fse0            = tendon_force_multiplier(Eps_t0);

    F_mA0           = fse0 / cos(Alpha) - fpe0;

    start_mn        = 0.0;
    end_mn          = 2.0;
    pre_error       = 1000;

    while (abs(real(pre_error)) > 1e-11)
        L_mn0 = 1/2 * (start_mn + end_mn);
        L_m0  = L_mn0 * L_mo;
        L_t0  = mod_L_mt0 * L_mt0 - L_m0 * cos(Alpha);
        Eps_t0 = L_t0 / L_ts - 1;

        fse0  = tendon_force_multiplier(Eps_t0);
        fpe0  = passive_muscle_force_length_multiplier(L_mn0);

        pre_error = fse0 - fpe0 * cos(Alpha);

        if pre_error < 0
            end_mn = L_mn0;
        else
            start_mn = L_mn0;
        end
    end

    L_mt0 = L_mn0 * cos(Alpha) * L_mo + L_t0;
    fprintf("%.4f N force was generated\n\n", fse0 * F_mo);
    
    fprintf("Initial fiber length: %.4f m\n", L_mn0 * L_mo);
    fprintf("%.4f %% lengthened\n\n", (L_mn0 * L_mo - L_mo) / L_mo * 100);

    fprintf("Initial tendon length: %.4f m\n", L_t0);
    fprintf("%.4f %% lengthened\n\n", (L_t0 - L_ts) / L_ts * 100);

    % [ calculations ]

    totalTime    = 30;                                                                                                                                                                                                                                                                                   ;
    dt           = 0.001;
    time         = 0 : dt : (totalTime - dt);

    frequencies  = logspace(-1, 2, 100);                                              % Sine-wave frequency array <- To supply the input signal on a per-single-frequency basis

    freq_len     = length(frequencies);          
    
    time_len     = length(time);
    
    xArray       = zeros(freq_len, time_len);                                          % Displacement array
    vArray       = zeros(freq_len, time_len);                                          % Velocity array
    aArray       = zeros(freq_len, time_len);                                          % Acceleration array (NOT AN ACTIAVTION)
    
    x_ext        = zeros(freq_len, time_len); %%%
    v_ext        = zeros(freq_len, time_len); %%%
    a_ext        = zeros(freq_len, time_len); %%%

    L_mt_preset  = zeros(freq_len, time_len);
    L_mt_eff     = zeros(freq_len, time_len);

    % [ initialization ]

    L_t          = ones(freq_len, time_len) .* L_t0;                                   % tendon length
    L_mn         = ones(freq_len, time_len) .* L_mn0;                              % normalized muscle length
    Eps_t        = ones(freq_len, time_len) .* Eps_t0;                                 % tendon strain
    f_l          = ones(freq_len, time_len) .* f_l0;                                   % active force(length) multiplier
    fse          = ones(freq_len, time_len) .* fse0;                                   % tendon force
    fpe          = ones(freq_len, time_len) .* fpe0;                                   % muscle passive force
    F_mA         = ones(freq_len, time_len) .* F_mA0;                                  % muscle active force
    penn         = ones(freq_len, time_len) .* Alpha;
    Fmtot        = zeros(freq_len, time_len);
    Fmtot_       = zeros(freq_len, time_len);


    sig_in       = zeros(freq_len, time_len);  
    sig_in_a     = zeros(freq_len, time_len);  
    sig_out      = zeros(freq_len, time_len);         

    L_mt         = ones(freq_len, time_len) * L_mt0;                                          % muscle tendon length
    A_mt         = zeros(freq_len, time_len);                                          % muscle tendon accleration
    V_mt         = zeros(freq_len, time_len);                                          % muscle tendon velocity
    V_m          = zeros(freq_len, time_len);                                          % muscle velocity
    V_t          = zeros(freq_len, time_len);
    
    % patch_result = cell(length(frequencies), 1);
    
    idx          = zeros(1, length(frequencies));

    pha          = zeros(1, length(frequencies));
    sin_f        = zeros(1, length(frequencies));
    mag          = zeros(1, length(frequencies));
    exc          = zeros(1, length(frequencies));

    F_ext_pass   = zeros(freq_len, time_len);
    
    fprintf("-----------------------------\n\n");
    fprintf("[ Main Calculation started ] \n\n");
    fprintf("-----------------------------\n\n");

    THD_vals        = zeros(1, freq_len, 'single');
    
    for k = 1 : length(frequencies)
        tic;

        fprintf("<strong>[ Patch number: %d ]\n\n</strong>", pnum);
        fprintf("<strong>Initial length   : %f\nInitial Velocity : %f</strong>\n\n", mod_L_mn0, mod_V_m0);
        fprintf("<strong>Step number: %d</strong> \n\n", k);
    
        freq = frequencies(k);

        amp                 = 0.005 * L_mt0;
        L_mt_preset(k, :)   = L_mt0 + amp * sin(2 * pi * freq * time);
        V_mt(k, :)          = amp * 2 * pi * freq * cos(2 * pi * freq * time);
        A_mt(k, :)          = -amp * (2 * pi * freq)^2 * sin(2 * pi * freq * time);

        sig_in(k, :)        = L_mt_preset(k, :);
        u                   = 0.5;
              
        temp_result         = NaN(time_len, 7);
        temp_index          = 1;
    
        for i = 1 : time_len - 1                                                                            % Explicit method
            
            L_mt_eff(k, i) = L_mt_preset(k, i);

            start_      = 0.0;       
            error       = 1000;                                                                             % error should be less than 1e-10

            end_        = 1.5;
            
            while ( abs(real(error)) > 1e-11 )                                                          % --------------Bisection method---------------
                L_mm = 1/2 * (start_ + end_);                                                               % normalized new muscle length (midpoint)
                
                % penn_temp       = dAlpha(L_mt_eff(k, i), L_t(k, i));
                penn_temp = Alpha;

                L_t  (k, i+1)   = L_mt_eff(k, i) - (L_mm * L_mo) * cos(penn_temp);                           % new muscle velocity (explicit)
                
                V_m  (k, i+1)   = ((L_mm - L_mn(k, i)) * L_mo / dt) / V_mmax;
                Eps_t(k, i+1)   = L_t(k, i+1) / L_ts - 1;                                                   % new tendon strain
                fse  (k, i+1)   = tendon_force_multiplier(Eps_t(k, i+1));                                   % new tendon force
                f_l  (k, i+1)   = active_muscle_force_length_multiplier(L_mm);                              % new active muscle force(length) multiplier
                fpe  (k, i+1)   = passive_muscle_force_length_multiplier(L_mm);                             % new muscle passive force
                F_mA (k, i+1)   = active_force(V_m(k, i+1), f_l(k, i+1), u);
                % F_mA (k, i+1)   = 0;

                error           = fse(k, i+1) - (F_mA(k, i+1) + fpe(k, i+1)) * cos(penn_temp);              % error for force equilbrium (muscle, tendon)

                if (error < 0) 
                    end_    = L_mm;
                else  
                    start_  = L_mm; 
                end                                                                                         % setting new midpoint
            end            
                                                                                                            % -----------------iteration over----------------
            L_mn(k, i+1)    = L_mm;                                                                         % save Bisection result
            
            if (L_mn(k, i+1) * L_mo < L_mo)
                sprintf("Fiber length is shorter than Lopt");
            end

            % penn(k, i+1)    = Alpha;

            penn(k, i+1)    = penn_temp;

            Fmtot_(k, i+1)   = F_mo * (F_mA(k, i+1) + fpe(k, i+1)) * cos(penn(k, i+1));
            Fmtot(k, i+1)    = Fmtot_(k, i+1) + mass * A_mt(k, i+1) - damp * V_mt(k, i+1);
            
            temp_result(temp_index, :) = [L_mn(k, i+1), V_m(k, i+1), F_mA(k, i+1), fpe(k, i+1), Fmtot(k, i+1), fse(k, i+1), L_mt_preset(k, i+1)];
            temp_index = temp_index + 1;
            
        end          

        patch_filename = sprintf('patch_result/p%d_k%d.mat', pnum, k);
        save(patch_filename, 'temp_result', '-v7.3');
        
        loaded = load(patch_filename, 'temp_result');
        temp_result = loaded.temp_result;

        sig_out(k, cutpoint:end-1) = temp_result(cutpoint:end-1, 5);

        Fs              = 1 / dt;

        n_x             = length(sig_in(k, cutpoint:end-2));
        n_f             = length(sig_out(k, cutpoint:end-2));
        f               = (0 : n_x - 1) * (Fs / n_x);
        
        X_fft           = fft(sig_in(k, cutpoint:end-2));
        F_fft           = fft(sig_out(k, cutpoint:end-2));

        [~, idx(k)]     = min(abs(f - freq));

        % omega   = 2 * pi * f(idx(k));
        % Z_total = F_fft(idx(k)) / X_fft(idx(k));
        % Z_visco = Z_total - mass * (1i * omega)^2;

        sin_f(k)        = freq;
        % mag(k) = abs(Z_visco);
        % pha(k) = angle(Z_visco);
        mag  (k)        = abs(F_fft(idx(k))) / abs(X_fft(idx(k)));
        pha  (k)        = angle(F_fft(idx(k))) - angle(X_fft(idx(k)));
        % exc  (k)        = abs(A_fft(idx(k)));
        

        fprintf("%d step was totally done \n\n", k);
        toc;
 
        fprintf("\n-----------------------------\n\n");

        N_fft     = n_f;                                                % 실제 FFT 길이
        P_out     = abs(F_fft(1 : floor(N_fft/2) + 1)).^2;              % 출력신호 절반 구간 파워 스펙트럼       

    end

    pha = unwrap(pha) * (180 / pi);
    
    visual_result{pnum}{:, 1} = sin_f;
    visual_result{pnum}{:, 2} = mag;
    visual_result{pnum}{:, 3} = pha;
    visual_result{pnum}{:, 4} = exc;

    index_V = index_V + 1;

    % calculation ends

    if index_V > aim_Vnum
        index_V = 1;
        index_mt = index_mt + 1;
    end

    pnum = pnum + 1;
end

%% Extracting simulation data into .mat file

% Need modify filename to adjust the simulation setting
% sol/gast (muscle type) | YB/OB (aging) | wod/d (damping)

% Format:: TMM_results_KMS_Lmn0_uo_sol_YB_wod.mat %

saveFolder1 = 'G:\Research\TMM test\MUSCLETEST\TMM_result(nopenn)_KMS\beforeFFT';

if ~exist(saveFolder1, 'dir')
    mkdir(saveFolder1);
end

for file_idx = 1 : length(L_mn0_values)
    
    fileName = [num2str(L_mn0_values(file_idx)) '_' num2str(u0) '_sol_YB_wd1000_bF_passive_30kg.mat'];
    
    fullPath = fullfile(saveFolder1, fileName);

    savingdata_time = cell(1, freq_len);
    
    for k = 1 : freq_len
        patch_filename = sprintf('patch_result/p%d_k%d.mat', file_idx, k);
        loaded = load(patch_filename, 'temp_result');
        savingdata_time{k} = loaded.temp_result;
    end

    save(fullPath, 'savingdata_time', '-v7.3');
end

saveFolder2 = 'G:\Research\TMM test\MUSCLETEST\TMM_result(nopenn)_KMS\afterFFT';

if ~exist(saveFolder2, 'dir')
    mkdir(saveFolder2);
end

for file_idx = 1 : length(L_mn0_values)
    
    fileName = [num2str(L_mn0_values(file_idx)) '_' num2str(u0) '_sol_YB_wd1000_aF_passive_30kg.mat'];
    
    fullPath = fullfile(saveFolder2, fileName);

    savingdata_freq = visual_result{file_idx};

    save(fullPath, 'savingdata_freq');
end

%% Bode plots

figure();
% sg = sgtitle(sprintf('Midpoint: L_{mn} = %.2f, V_{m} = %.2f, at u_{0} = %.2f', mod_L_mn0, mod_V_m0, u0));
sg              = sgtitle(sprintf('Bode plot of each conditions w/ u_{0} = %.2f', u0), 'FontSize', 18);

cm              = hsv(pnum - 1);
cutoff_level    = -3;

nat_freqs       = zeros(pnum-1, 1);

for fnum = 1 : pnum - 1
    
    freq_values = visual_result{fnum}{:, 1};
    mag_values = 20 * log10(visual_result{fnum}{:, 2});
    valid_idx = find(freq_values > 1);

    [max_mag, rel_max_idx] = max(mag_values(valid_idx));
    nat_freq = freq_values(valid_idx(rel_max_idx));
    nat_freqs(fnum)    = nat_freq;

    subplot(2, 1, 1);
    semilogx(visual_result{fnum}{:, 1}, mag_values, 'o', 'MarkerSize', 1, ...
        'Color', cm(fnum, :), 'DisplayName', sprintf('L_{mn} = %.2f, V_{m} = %.2f', L_mn0_values(fnum), V_m0_values(1)));
    xlabel("Frequency (Hz)", 'FontSize', 15); ylabel("Magnitude (dB)", 'FontSize', 13); grid on;
    % ylim([0 100]);
    hold on;

    subplot(2, 1, 2);
    semilogx(visual_result{fnum}{:, 1}, visual_result{fnum}{:, 3}, 'o', 'MarkerSize', 1, ...
        'Color', cm(fnum, :), 'DisplayName', sprintf('L_{mn} = %.2f, V_{m} = %.2f', L_mn0_values(fnum), V_m0_values(1)));
    xlabel("Frequency (Hz)", 'FontSize', 15); ylabel("Phase (deg)", 'FontSize', 13); grid on;
    % ylim([-180 180]);
    hold on;
end

hold off;

%% Patching (Multi)
% TO VISUALIZE ALL THE PATCHES

% L_mn_patch_midval = L_mn0_values;
% V_m_patch_midval  = V_m0_values;
% 
% l_pnum  = length(L_mn_patch_midval);
% v_pnum  = length(V_m_patch_midval);
% 
% grid_interval = 100;
% 
% % fprintf("Total number of the patches: %d\n\n", l_pnum * v_pnum);
% 
% temp = cell(v_pnum, l_pnum);
% 
% for l = 1 : l_pnum
%     for v = 1 : v_pnum
%         L_mn_patch_temp_range = linspace(L_mn_patch_midval(l) - L_range, L_mn_patch_midval(l) + L_range, grid_interval);
%         V_m_patch_temp_range = linspace(V_m_patch_midval(v) - V_range, V_m_patch_midval(v) + V_range, grid_interval);
% 
%         [L_mn_patch_temp, V_m_patch_temp] = meshgrid(L_mn_patch_temp_range, V_m_patch_temp_range);
% 
%         F_patch_temp = zeros(grid_interval, grid_interval);
% 
%         for i = 1 : grid_interval
%             f_l_patch_temp = active_muscle_force_length_multiplier(L_mn_patch_temp(i));
%             F_pe_patch_temp = passive_muscle_force_length_multiplier(L_mn_patch_temp(i));
% 
%             for j = 1 : grid_interval
%                 F_patch_temp(i, j) = active_force(V_m_patch_temp(i), f_l_patch_temp, 0.5) + F_pe_patch_temp;
%             end
%         end
% 
%         temp{v, l} = {L_mn_patch_temp, V_m_patch_temp, F_patch_temp};
%     end
% end
% 
% figure();
% 
% for i = 1 : l_pnum
%     for j = 1 : v_pnum
%         temp_data = temp{j, i};
%         scatter3(temp_data{1}, temp_data{2}, temp_data{3}, 2, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'red', 'MarkerFaceAlpha', 0.1);
%         hold on;
%     end
% end
% 
% % [ Visualize the result when calculated f, l, v data satisfy the criteria (remaining in a patch) ]
% 
% colors = jet(freq_len);
% 
% count = 0;
% 
% for p_num_ = 1 : pnum - 1
%     for k = 1 : freq_len
%         color = colors(k, :);
% 
%         plot3(patch_result{p_num_, k}(cutpoint:end, 1), patch_result{p_num_, k}(cutpoint:end, 2), ...
%             patch_result{p_num_, k}(cutpoint:end, 3), 'o', 'Color', color, 'MarkerSize', 2);
%         count = count + 1;
% 
%         hold on;
%     end
% end
% 
% 
% L_mn_baseline_range = linspace(0.2, 1.8, 100);
% V_m_baseline_range = linspace(-2.0, 2.0, 100);
% [L_mn_grid, V_m_grid] = meshgrid(L_mn_baseline_range, V_m_baseline_range);
% 
% f_l_baseline_grid = zeros(size(L_mn_grid));
% F_pe_baseline_grid = zeros(size(L_mn_grid));
% 
% for i = 1 : numel(L_mn_grid)
%     f_l_baseline_grid(i) = active_muscle_force_length_multiplier(L_mn_grid(i));
%     F_pe_baseline_grid(i) = passive_muscle_force_length_multiplier(L_mn_grid(i));
% end
% 
% F_baseline_total_grid = zeros(size(L_mn_grid));
% for i = 1 : numel(L_mn_grid)
%     F_baseline_total_grid(i) = (active_force(V_m_grid(i), f_l_baseline_grid(i), 0.5) + F_pe_baseline_grid(i));
% end
% 
% figure();
% surf(L_mn_grid, V_m_grid, F_baseline_total_grid, 'FaceAlpha', 0.1);
% xlabel('Normalized Muscle Length (L/L_{mo})', 'FontSize', 18);
% ylabel('Muscle Velocity (V/V_{max})', 'FontSize', 18);
% zlabel('Muscle Force', 'FontSize', 18);
% title('Patch with base F-L-V Surface(act = 0.5)', 'FontSize', 18);
% 
% hold off;

%% Extra figure extracting
target_freq_num = 1;

figure();
subplot(2, 1, 1);
plot(time, savingdata_time{1, target_freq_num}(:, 4) * F_mo); title('Fiber Passive Force'); ylabel('(N)');

subplot(2, 1, 2);
plot(time, savingdata_time{1, target_freq_num}(:, 6) * F_mo); title('Tendon Force'); ylabel('(N)');


figure();
subplot(3, 1, 1);
plot(time, savingdata_time{1, target_freq_num}(:, 1) * L_mo); title('Fiber Length'); ylabel('(m)');

subplot(3, 1, 2);
plot(time, L_t(target_freq_num, :)'); title('Tendon Length'); ylabel('(m)');

subplot(3, 1, 3);
plot(time, savingdata_time{1, target_freq_num}(:, 7)); title('Total Length'); ylabel('(m)');

figure();
plot(time, V_m(target_freq_num, :)); title('Fiber velocity'); ylabel('(m/s)');

% figure();
% plot(time, ((patch_result{1, target_freq_num}(:, 3) + patch_result{1, target_freq_num}(:, 4)) * cos(Alpha) ...
%     - fse(target_freq_num, :)) * F_mo); title('Force Difference (N)');

%% Functions


function Tau = get_Tau(a, u)
    Tau_a   = 0.01;
    % Tau_a   = 0.005;
    Tau_d   = 0.04;                                                                         % Equation 2 of thelen(2003)
    % Tau_d   = 0.06;

    if (u > a)
        Tau = Tau_a * (0.5 + 1.5 * a);
    else
        Tau = Tau_d / (0.5 + 1.5 * a);
    end
end


function Eps_t = tendon_strain_decider2(F_n)
    EPSto   = 0.033;                                                                        % Tendon strain due to maximum isometric force
    EPSttoe = 0.609 * EPSto;                                                                % Tendon strain above which the tendon exhibits linear behavior
    F_ttoe  = 0.33;                                                                         % Normalized tendon force at the transition from nonlinear to linear behavior
    K_toe   = 3;                                                                            % Tendon exponential shape factor
    K_lin   = 1.712 / EPSto;                                                                % Equation 5 of thelen(2003)
    
    Eps_t1 = EPSttoe / K_toe * log(F_n / F_ttoe * (exp(K_toe) - 1) + 1);                    % tendon strain with exponential behavior
    Eps_t2 = (F_n - F_ttoe) / K_lin + EPSttoe;                                              % tendon strain with linear behavior
    
    if (F_n <= F_ttoe)
        Eps_t = Eps_t1;
    else
        Eps_t = Eps_t2;
    end
end


function F_pe = passive_muscle_force_length_multiplier(L_mn)                                
    K_pe    = 4.0;                                                                         
    EPSmo   = 0.6;        
    % EPSmo   = 0.5;

    if (L_mn > 1.0)
        F_pe = (exp((K_pe * (L_mn - 1)) / EPSmo) - 1) / (exp(K_pe) - 1);
    else
        F_pe = 0;
    end
end

function f_l0 = active_muscle_force_length_multiplier(L_mn)                                 
    Gamma = 0.45;                                                                           % Equation 4 of thelen (2003)
    f_l0 = exp(-((L_mn - 1).^2) / Gamma);
end


function F_mn = active_force(V_temp, f_l, a)                                                % Equation 6,7 of thelen (2003)
    A_f     = 0.3;                                                                          
    F_mnlen = 1.4;      
    % F_mnlen = 1.8;

    velocity_norm = V_temp / (0.25 + 0.75 * a);
   
    if (V_temp <= 0)
        % When V is negative
        velocity_factor = (1 + velocity_norm) / (1 - velocity_norm / A_f);                  % Same as force-velocity multiplier (fv)
    else
        % When V is positive
        someth = velocity_norm * (2 + 2 / A_f) / (F_mnlen - 1);
        velocity_factor = (1 + F_mnlen * someth) / (1 + someth);
    end
    F_mn = a * f_l * velocity_factor;
end


function F_t = tendon_force_multiplier(eps_t)
    % global EPSttoe F_ttoe K_toe K_lin;                                                    % Equation 5 of thelen (2003)

    EPSto   = 0.033;                                                                        % Tendon strain due to maximum isometric force
    EPSttoe = 0.609 * EPSto;                                                                % Tendon strain above which the tendon exhibits linear behavior
    F_ttoe  = 0.33;                                                                         % Normalized tendon force at the transition from nonlinear to linear behavior
    K_toe   = 3;                                                                            % Tendon exponential shape factor
    K_lin   = 1.712 / EPSto;       
   
    if (eps_t <= EPSttoe)
        F_t = F_ttoe / (exp(K_toe) - 1) * (exp(K_toe*eps_t / EPSttoe) - 1);
    else
        F_t = K_lin * (eps_t - EPSttoe) + F_ttoe;
    end 
end

function dPennationAngle = dAlpha(m_muscleLength, tendon_length)
    
    Alpha0  = 0.4363;               % Optimal pennation angle (rad)
    % Alpha0   = 0.3050;
    L_mo    = 0.05;                 % Optimal fiber length (m)

    m_const_muscle_height = L_mo * sin(Alpha0);

    cos_fiber_length = m_muscleLength - tendon_length;

    % m_muscleLength is the total muscle length
    if (m_const_muscle_height / cos_fiber_length) > 9.0
        dPennationAngle = atan(9.0);
    else
        dFiberLength = sqrt(cos_fiber_length * cos_fiber_length + m_const_muscle_height * m_const_muscle_height);
        dPennationAngle = asin(m_const_muscle_height / dFiberLength);
    end
end