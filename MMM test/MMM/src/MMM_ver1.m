% By Minseung Kim, 2024-12-27
% profile on
clc;
clear;

%% Problem Assumptions

% Gathering muscle info 

maximumNormalizedFiberVelocity          = 10;                % in units of norm fiber lengths/second
% maximumNormalizedFiberVelocity          = 8;                   % OB data
maximumPennationAngle                   = 89 * (pi/180); 
    
flag_useArnold2010SoleusArchitecture    = 1;
flag_useThelen2003SoleusArchitecture    = 1;
muscleAbbrArnold2010                    = 'soleus';

flag_plotNormMuscleCurves               = 0;
flag_updateCurves                       = 1;

muscleAbbr  = [];
if(flag_useArnold2010SoleusArchitecture == 1)
    muscleAbbr = muscleAbbrArnold2010;
else
    muscleAbbr = 'compBench';
end

tendonStrainAtOneNormForce = 0.033; %If this is greater than 0 this value will
                                 %be used to make the tendon-force-length
                                 %curve. Otherwise the default of 0.049 is
                                 %taken.

normMuscleCurves = createDefaultNormalizedMuscleCurves(muscleAbbr,...
                                        tendonStrainAtOneNormForce,...
                                        flag_updateCurves,...
                                        flag_plotNormMuscleCurves);

muscleName  = [];
fiso        = [];
lceOpt      = [];
alphaOpt    = [];
ltSlk       = [];

if(flag_useArnold2010SoleusArchitecture == 1)
    unitsMKSN = 1;
    arnold2010LegArch = getArnold2010LegMuscleArchitecture(unitsMKSN);

    idx_         =  getArnold2010MuscleIndex(muscleAbbrArnold2010,...
                          arnold2010LegArch.abbrevation);
    
    muscleName  = arnold2010LegArch.names{idx_};
    fiso        = arnold2010LegArch.peakForce(idx_);
    lceOpt      = arnold2010LegArch.optimalFiberLength(idx_);
    alphaOpt    = arnold2010LegArch.pennationAngle(idx_);
    ltSlk       = arnold2010LegArch.tendonSlackLength(idx_);

else 
    muscleName  = 'MMM';
    fiso        = 1;
    lceOpt      = 0.02;
    alphaOpt    = 30*(pi/180);
    ltSlk       = 0.20;
end

if(flag_useThelen2003SoleusArchitecture == 1)
    % muscleName      = 'Soleus';
    % fiso            = 3549.0;       % (N)
    % lceOpt          = 0.05;         % (m)
    % alphaOpt        = 0.4363;       % (rad)
    % ltSlk           = 0.25;         % (m)
    % ltSlk           = 0.20;         % (m)
    
    % muscleName      = 'Soleus_OB';
    % fiso            = 2484.3;       % (N)
    % lceOpt          = 0.05;         % (m)
    % alphaOpt        = 0.3050;       % (rad)
    % ltSlk           = 0.25;         % (m)
    % 
    % muscleName      = 'Gastroc';
    % fiso            = 1558.0;       % (N)
    % lceOpt          = 0.06;         % (m)
    % alphaOpt        = 0.2967;       % (rad)
    % ltSlk           = 0.39;         % (m)

    muscleName      = 'TA';
    fiso            = 905.0;        % (N)
    lceOpt          = 0.0898;       % (m)
    alphaOpt        = 0.0873;       % (rad)
    ltSlk           = 0.2043;       % (m) 

    % muscleName      = 'bifemsh';
    % fiso            = 804.0;        % (N)
    % lceOpt          = 0.1730;       % (m)
    % alphaOpt        = 0.4014;       % (rad)
    % ltSlk           = 0.0890;       % (m)
end

muscleArch = [];
muscleArch.name                             = muscleName;
muscleArch.abbr                             = muscleAbbr;
muscleArch.fiso                             = fiso;
muscleArch.optimalFiberLength               = lceOpt;
muscleArch.maximumNormalizedFiberVelocity   = maximumNormalizedFiberVelocity;
muscleArch.pennationAngle                   = alphaOpt;
muscleArch.tendonSlackLength                = ltSlk;
    
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

%% Generating patch_info.txt

% filename                = 'patch_info.txt';
filename                = 'single_patch.txt';

fileID                  = fopen(filename, 'w');
fprintf(fileID, 'L_mn0 V_m0\n');

% L_mn0_values            = [0.6, 0.8, 1.0, 1.2, 1.4];
L_mn0_values            = [1.0];
V_m0_values             = [0.0];

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

patch_info              = readtable('single_patch.txt');

index_L                 = 1;
index_V                 = 1;
pnum                    = 1;

L_range                 = 0.100;                                                                     % Muscle length +- range
V_range                 = 0.100;                                                                     % Velocity +- range

aim_Lnum                = length(L_mn0_values);
aim_Vnum                = length(V_m0_values);
aim_pnum                = height(patch_info);

modelConfig.useElasticTendon    = 1;
modelConfig.useFiberDamping     = 0;
modelConfig.damping             = 0.1;
modelConfig.minActivation       = 0.01;
modelConfig.iterMax             = 10000;
modelConfig.tol                 = 1e-10;
modelConfig.passiveOnlyMode     = false;

patch_result            = cell(aim_pnum, 100, 1);
visual_result           = cell(aim_pnum, 1);

while (index_L <= aim_Lnum)

    const_b             = 0.01;
    cutpoint            = 1;                                                              % cutpoint = 1 means no data cutting
        
    mod_V_m0            = V_m0_values(index_V);
    mod_L_mn0           = L_mn0_values(index_L);

    lceAT0              = mod_L_mn0 * lceOpt;
    mod_L_mt0           = lceAT0 * cos(alphaOpt) + ltSlk;
    fprintf('\n%.5f: mod_L_mt0 (desired L_mt0)\n', mod_L_mt0);

    lt0                 = mod_L_mt0 - lceAT0 * cos(alphaOpt);
    dlceAT0             = 0.0;

    V_mt0               = 0;                                                              % muscle-tendon velocity

    u0                  = 0.5;
    
    desiredLf           = lceAT0;
    F_ext_equil         = findExternalForceForFiberLength(desiredLf, u0, muscleArch, ...
        normMuscleCurves, modelConfig, alphaOpt, ltSlk);
    
    tendonForce_initial = F_ext_equil;
    
    fiberForce_initial  = tendonForce_initial / cos(alphaOpt);
    
    pre_L_mt0           = desiredLf * cos(alphaOpt) + ltSlk; 
    pre_dlceAT0         = 0.0;
    pre_lceAT0          = desiredLf;

    mass                = 3000000;           
    % mass                = 300;

    % [ calculations ]
    % Configuration for Proposal: totalTime = 120; frequencies = logspace(-1, 2, 400)

    totalTime           = 120;                                                                                                                                                                                                                                                                                   ;
    dt                  = 0.001;
    time                = 0 : dt : (totalTime - dt);

    frequencies         = logspace(-1, 2, 400);                                              
    
    freq_len            = length(frequencies);          
    time_len            = length(time);
    ext_damping         = 0.0;
    
    % [ initialization ]

    sig_in              = zeros(freq_len, time_len);  
    sig_in_a            = zeros(freq_len, time_len);  
    sig_out             = zeros(freq_len, time_len - cutpoint);         

    L_mt                = ones(freq_len, time_len) .* mod_L_mt0;                              % muscle tendon length
    A_mt                = zeros(freq_len, time_len);                                          % muscle tendon accleration
    V_mt                = zeros(freq_len, time_len);                                          % muscle tendon velocity
    
    fiberForce          = zeros(freq_len, time_len) .* fiberForce_initial;
    tendonForce         = ones(freq_len, time_len) .* tendonForce_initial;
    activefiberForce    = zeros(freq_len, time_len);
    passivefiberForce   = zeros(freq_len, time_len);

    force_diff          = zeros(freq_len, time_len);
    lceAT               = ones(freq_len, time_len) .* lceAT0;
    dlceAT              = zeros(freq_len, time_len);

    fiber_length        = zeros(freq_len, time_len);
    fiber_velocity      = zeros(freq_len, time_len);
    tendon_length       = ones(freq_len, time_len) .* ltSlk;
    tendon_velocity     = zeros(freq_len, time_len);
    alpha               = zeros(freq_len, time_len);

    idx                 = zeros(1, length(frequencies));
    pha                 = zeros(1, length(frequencies));
    sin_f               = zeros(1, length(frequencies));
    mag                 = zeros(1, length(frequencies));
    exc                 = zeros(1, length(frequencies));

    THD_vals            = zeros(1, freq_len, 'single');
    
    fprintf("-----------------------------\n\n");
    fprintf("[ Main Calculation started ] \n\n");
    fprintf("-----------------------------\n\n");
    
    for k = 1 : length(frequencies)
        tic;

        fprintf("<strong>[ Patch number: %d ]\n\n</strong>", pnum);
        fprintf("<strong>Initial length   : %f\nInitial Velocity : %f</strong>\n\n", mod_L_mn0, mod_V_m0);
        fprintf("<strong>Step number: %d</strong> \n\n", k);
    
        freq = frequencies(k);
        
        % u_Vec    = calcSinusoidState(time, u0, const_b, freq);
        % u        = u_Vec(2);
        u    = sinwave(freq, time, u0, u0, const_b);
        
        for u_temp = 1 : length(u)
            if u(u_temp) < 0; u(u_temp) = 0; end
            if u(u_temp) > 1.0; u(u_temp) = 1.0; end
        end

        a    = ones(1, length(u)) * u0;

        for act = 1 : length(u) - 1
            da_dt = (u(act) - a(act)) / get_Tau(a(act), u(act));
            a(act + 1) = a(act) + dt * da_dt;
        end
        
        sig_in(k, :)    = u;
        sig_in_a(k, :)  = a;
        
        temp_result     = NaN(time_len, 7);
        temp_index      = 1;

        L_mt(k, 1)      = pre_L_mt0;
        V_mt(k, 1)      = 0;
        lceAT(k, 1)     = pre_lceAT0 * cos(alphaOpt);
        dlceAT(k, 1)    = pre_dlceAT0;

        for i = 1 : time_len - 1                                                            % Explicit method

            pathState        = [V_mt(k, i); L_mt(k, i)];
            % muscleState      = [dlceAT(k, i); lceAT(k, i)];
            muscleState      = lceAT(k, i);

            mtInfo           = calcMillard2012DampedEquilibriumMuscleInfo(...
                sig_in_a(k, i), ...
                pathState, ...
                muscleState, ...
                muscleArch, ...
                normMuscleCurves, ...
                modelConfig);
            
            tendonForce(k, i+1)       = mtInfo.muscleDynamicsInfo.tendonForce;
            fiberForce(k, i+1)        = mtInfo.muscleDynamicsInfo.fiberForceAlongTendon;
            activefiberForce(k, i+1)  = mtInfo.muscleDynamicsInfo.activeFiberForce;
            passivefiberForce(k, i+1) = mtInfo.muscleDynamicsInfo.passiveFiberForce;

            % lceAT(k, i+1)             = mtInfo.state.value;
            dlceAT(k, i+1)            = mtInfo.state.derivative;
            lceAT(k, i+1)             = lceAT(k, i) + dlceAT(k, i+1) * dt;
            
            % lceAT(k, i+1)             = mtInfo.muscleLengthInfo.fiberLengthAlongTendon;
            % dlceAT(k, i+1)            = mtInfo.fiberVelocityInfo.fiberVelocityAlongTendon;
            alpha(k, i+1)             = mtInfo.muscleLengthInfo.pennationAngle;

            tendon_length(k, i+1)     = mtInfo.muscleLengthInfo.tendonLength;
            tendon_velocity(k, i+1)   = mtInfo.fiberVelocityInfo.tendonVelocity;

            A_mt(k, i+1)    = (F_ext_equil - tendonForce(k, i+1) - ext_damping * V_mt(k, i)) / mass;
            V_mt(k, i+1)    = V_mt(k, i) + A_mt(k, i+1) * dt;                       
            L_mt(k, i+1)    = L_mt(k, i) + V_mt(k, i+1) * dt;       

            temp_result(temp_index, :) = [lceAT(k, i+1), dlceAT(k, i+1), activefiberForce(k, i+1), ... 
                passivefiberForce(k, i+1), fiberForce(k, i+1), tendonForce(k, i+1), L_mt(k, i+1)];

            temp_index = temp_index + 1;
        end          

        % patch_result{pnum, k} = temp_result;
        patch_filename = sprintf('patch_result/p%d_k%d.mat', pnum, k);
        save(patch_filename, 'temp_result', '-v7.3');

        loaded = load(patch_filename, 'temp_result');
        temp_result = loaded.temp_result;
   
        % sig_out(k, :)   = patch_result{pnum, k}(cutpoint:end-1, 7);
        sig_out(k, :)   = temp_result(cutpoint:end-1, 5);
        
        A_fft           = fft(sig_in(k, cutpoint:end-1)); 
        F_fft           = fft(sig_out(k, cutpoint:end-1));
        
        Fs              = 1 / dt;
        n_x             = length(sig_in(k, cutpoint:end-1));
        n_f             = length(sig_out);
        f               = (0 : n_x - 1) * (Fs / n_x);
    
        [~, idx(k)]     = min(abs(f - freq));
    
        sin_f(k)        = freq;
        mag  (k)        = abs(F_fft(idx(k)) / A_fft(idx(k)));
        pha  (k)        = angle(F_fft(idx(k))) - angle(A_fft(idx(k)));
        exc  (k)        = abs(A_fft(idx(k)));

        fprintf("%d step was totally done \n\n", k);
        toc;
 
        fprintf("\n-----------------------------\n\n");

        % THD Calculation
        N_fft       = n_f;                                                % 실제 FFT 길이
        P_out       = abs(F_fft(1 : floor(N_fft/2) + 1)).^2;              % 출력신호 절반 구간 파워 스펙트럼
        
        basic_idx   = idx(k);                                             % 기본파 인덱스
        basic_power = P_out(basic_idx);                                   % 기본파 파워
        
        har_idx     = 2 * basic_idx : basic_idx : floor(N_fft/2) + 1; 
        har_idx     = har_idx(har_idx <= floor(N_fft/2) + 1);             % 범위 초과 방지
        
        if basic_power ~= 0
            har_power    = sum(P_out(har_idx));
            THD_vals(k)  = sqrt(har_power) / sqrt(basic_power);         % THD = sqrt(고조파 파워합) / sqrt(기본파 파워)
        else
            THD_vals(k)  = NaN;                                         % 기본파가 너무 작으면 계산 불가
        end
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
        index_L = index_L + 1;
    end

    pnum = pnum + 1;
end

%% Extracting simulation data into .mat file

% Need modify filename to adjust the simulation setting
% sol/gast (muscle type) | YB/OB (aging) | wod/d (damping)

% Format:: MMM_results_KMS_Lmn0_uo_sol_YB_wod.mat %

saveFolder1 = 'G:\Research\MMM test\MMM\src\MMM_result_KMS\beforeFFT';

if ~exist(saveFolder1, 'dir')
    mkdir(saveFolder1);
end

for file_idx = 1 : length(L_mn0_values)
    
    fileName = [num2str(L_mn0_values(file_idx)) '_' num2str(u0) '_TA_YB_wod_bF_sta.mat'];
    
    fullPath = fullfile(saveFolder1, fileName);

    % savingdata_time = patch_result{file_idx};
    savingdata_time = cell(1, freq_len);
    
    for k = 1 : freq_len
        patch_filename = sprintf('patch_result/p%d_k%d.mat', file_idx, k);
        loaded = load(patch_filename, 'temp_result');
        savingdata_time{k} = loaded.temp_result;
    end

    save(fullPath, 'savingdata_time', '-v7.3');
end

saveFolder2 = 'G:\Research\MMM test\MMM\src\MMM_result_KMS\afterFFT';

if ~exist(saveFolder2, 'dir')
    mkdir(saveFolder2);
end

for file_idx = 1 : length(L_mn0_values)
    
    fileName = [num2str(L_mn0_values(file_idx)) '_' num2str(u0) '_TA_YB_wod_aF_sta.mat'];
    
    fullPath = fullfile(saveFolder2, fileName);

    savingdata_freq = visual_result{file_idx};

    save(fullPath, 'savingdata_freq');
end


%% Debugging plot

% figure(); plot(time, lceAT(1, :) ./ cos(alpha(1, :))); title('Fiber length along Tendon');

% figure(); plot(time, tendon_length(1, :)); title('Tendon length');

figure(); plot(time, fiberForce(1, :)); title('Fiber force');

figure(); plot(time, tendonForce(1, :)); title('Tendon force');

%% Bode plots with fitted 2nd order system

figure();
sg = sgtitle(sprintf('Bode plot of each conditions w/ u_{0} = %.2f', u0));

cm = hsv(pnum - 1);

cutoff_level = -3;

nat_freqs = zeros(pnum-1, 1);

for fnum = 1 : pnum - 1

    bode_freq = visual_result{fnum}{:, 1}; % Frequency (Hz)
    bode_mag = visual_result{fnum}{:, 2}; % Magnitude (linear)
    bode_pha = visual_result{fnum}{:, 3}; % Phase (degrees)

    mag_dB = 20 * log10(bode_mag);
    omega = 2 * pi * bode_freq; % Frequency (rad/s)

    % Generate complex response for FRD
    bode_response = bode_mag .* (cosd(bode_pha) + 1j * sind(bode_pha));
    data = frd(bode_response, omega);

    sys = tfest(data, 2);
    [wn, zeta] = damp(sys);
    nat_freq_Hz = wn(1) / (2 * pi);

    subplot(2, 1, 1);
    semilogx(bode_freq, mag_dB, 'o', 'MarkerSize', 1, ...
        'Color', cm(fnum, :), 'DisplayName', sprintf('Original: L_{mn} = %.2f, V_{m} = %.2f', L_mn0_values(fnum), V_m0_values(1)));
    hold on;
    % plot(closest_freq, closest_mag_dB, 'x', 'MarkerSize', 10, ...
        % 'Color', cm(fnum, :), 'LineWidth', 2);
    % text(closest_freq, closest_mag_dB, sprintf(' f_n = %.2f Hz', closest_freq), ...
        % 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    grid on;

    % Plot phase
    subplot(2, 1, 2);
    semilogx(bode_freq, bode_pha, 'o', 'MarkerSize', 1, ...
        'Color', cm(fnum, :), 'DisplayName', sprintf('Original: L_{mn} = %.2f, V_{m} = %.2f', L_mn0_values(fnum), V_m0_values(1)));
    hold on;
    grid on;
end

% Add legends
% subplot(2, 1, 1);
% legend('Location', 'EastOutside', 'Box', 'off');
% 
% subplot(2, 1, 2);
% legend('Location', 'EastOutside', 'Box', 'off');

hold off;

% Display natural frequencies
disp('Natural frequencies for each patch: ');
disp(nat_freqs);

% Total Harmonic Distortion (THD) plot
figure();
semilogx(frequencies, 20*log10(THD_vals), 'k-o', 'MarkerSize', 4);
title('Total Harmonic Distortion (THD)');
xlabel('Frequency (Hz)');
ylabel('THD (dB)'); ylim([-60 0]);
grid on;

%% Patching (Multi)
% TO VISUALIZE ALL THE PATCHES

% L_mn_patch_midval = L_mn0_values;
% V_m_patch_midval  = V_m0_values;
% 
% l_pnum  = length(L_mn_patch_midval);
% v_pnum  = length(V_m_patch_midval);
% 
% grid_interval = 150;
% 
% fprintf("\nGenerating F-L-V patch data...\n");
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
%             for j = 1 : grid_interval
%                 lce_temp = L_mn_patch_temp(i, j) * lceOpt;
%                 dlce_temp = V_m_patch_temp(i, j) * maximumNormalizedFiberVelocity * lceOpt;
% 
%                 pathState = [dlce_temp; lce_temp * cos(alphaOpt) + ltSlk];
%                 muscleState = [dlce_temp; lce_temp];
% 
%                 mtInfo = calcMillard2012DampedEquilibriumMuscleInfo( ...
%                     0.5, ...
%                     pathState, ...
%                     muscleState, ...
%                     muscleArch, ...
%                     normMuscleCurves, ...
%                     modelConfig);
% 
%                 F_patch_temp(i, j) = mtInfo.muscleDynamicsInfo.fiberForceAlongTendon;
%             end
%         end
% 
%         temp{v, l} = {L_mn_patch_temp, V_m_patch_temp, F_patch_temp};
%     end
% end
% 
% figure();
% for i = 1 : l_pnum
%     for j = 1 : v_pnum
%         temp_data = temp{j, i};
%         scatter3(temp_data{1}, temp_data{2}, temp_data{3}, 2, 'MarkerEdgeColor', 'none', ...
%                  'MarkerFaceColor', 'blue', 'MarkerFaceAlpha', 0.1);
%         hold on;
%     end
% end
% 
% colors = jet(freq_len);
% count = 0; 
% 
% for p_num_ = 1 : pnum - 1
%     for k = 1 : freq_len
%         color = colors(k, :);
% 
%         plot3(patch_result{p_num_, k}(:, 1) / lceOpt, patch_result{p_num_, k}(:, 2), patch_result{p_num_, k}(:, 5), ...
%               'o', 'Color', color, 'MarkerSize', 2);
%         count = count + 1;
% 
%         hold on;
%     end
% end
% 
% L_mn_baseline_range = linspace(0.2, 1.8, 100); 
% V_m_baseline_range = linspace(-2.0, 2.0, 100); 
% 
% [L_mn_grid, V_m_grid] = meshgrid(L_mn_baseline_range, V_m_baseline_range);
% 
% F_baseline_total_grid = zeros(size(L_mn_grid));
% 
% fprintf("\nCalculating Baseline F-L-V Surface...\n");
% 
% for i = 1 : size(L_mn_grid, 1)
%     for j = 1 : size(L_mn_grid, 2)
%         lce_temp = L_mn_grid(i, j) * lceOpt;
%         dlce_temp = V_m_grid(i, j) * maximumNormalizedFiberVelocity * lceOpt;
% 
%         pathState = [dlce_temp; lce_temp * cos(alphaOpt) + ltSlk];
%         muscleState = [dlce_temp; lce_temp];
% 
%         mtInfo = calcMillard2012DampedEquilibriumMuscleInfo( ...
%             0.5, ...
%             pathState, ...
%             muscleState, ...
%             muscleArch, ...
%             normMuscleCurves, ...
%             modelConfig);
% 
%         F_baseline_total_grid(i, j) = mtInfo.muscleDynamicsInfo.fiberForceAlongTendon;
%     end
% end
% 
% surf(L_mn_grid, V_m_grid, F_baseline_total_grid, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% xlabel('Normalized Muscle Length (L/L_{mo})', 'FontSize', 14);
% ylabel('Muscle Velocity (V/V_{max})', 'FontSize', 14);
% zlabel('Muscle Force (N)', 'FontSize', 14);
% title('Baseline F-L-V Surface (Millard Model)', 'FontSize', 16);
% grid on;
% 
% hold off;

fs = 8000;
freqs = [440, 880, 660]; 
duration = 0.2;
for freq = freqs
    t = 0:1/fs:duration;
    y = sin(2 * pi * freq * t);
    sound(y, fs);
    pause(duration); 
end

% profile viewer

%% Functions

function input = sinwave(frq, t, a0, a, b)                                                        % Generating sine-wave fore acvitation
    % input = (sin(2 * pi * frq * t) + 1) / 2                                               % Starts at a0 when the F_ext given
    
    phase_shift = asin(2 * a0 - 1);
    
    input = (a + sin(2 * pi * frq * t + phase_shift) * b);
end

function Tau = get_Tau(a, u)
    Tau_a   = 0.01;                                                                    
    Tau_d   = 0.04;                                                                         % Equation 2 of thelen(2003)
    % Tau_d   = 0.06;

    if (u > a)
        Tau = Tau_a * (0.5 + 1.5 * a);
    else
        Tau = Tau_d / (0.5 + 1.5 * a);
    end
end 

function F_ext_equil = findExternalForceForFiberLength( ...
    desiredLf, ...   
    u0, ...         
    muscleArch, ...
    normMuscleCurves, ...
    modelConfig, ...
    alphaOpt, ...    
    ltSlk)     
       
    L_mt0  = desiredLf * cos(alphaOpt) + ltSlk; 

    pathState   = [0, L_mt0];
    muscleState = [0, desiredLf * cos(alphaOpt)];

    mtInfo = calcMillard2012DampedEquilibriumMuscleInfo( ...
        u0, ...
        pathState, ...
        muscleState, ...
        muscleArch, ...
        normMuscleCurves, ...
        modelConfig);

    currentFiberForce = mtInfo.muscleDynamicsInfo.fiberForceAlongTendon;
    
    F_ext_equil       = currentFiberForce;
end