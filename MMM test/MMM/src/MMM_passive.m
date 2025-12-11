% By Minseung Kim, 2024-12-27
% 2025-12-11

% clc;
% clear;

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
    muscleName      = 'Soleus';
    fiso            = 3549.0;       % (N)
    lceOpt          = 0.05;         % (m)
    alphaOpt        = 0.4363;       % (rad)
    ltSlk           = 0.25;         % (m)
    % ltSlk           = 0.20;         % (m)
    
    % muscleName      = 'Soleus_OB';
    % fiso            = 2484.3;       % (N)
    % lceOpt          = 0.05;         % (m)
    % alphaOpt        = 0.3050;       % (rad)
    % ltSlk           = 0.25;         % (m)

    % muscleName      = 'Gastroc';
    % fiso            = 1558.0;       % (N)
    % lceOpt          = 0.06;         % (m)
    % alphaOpt        = 0.2967;       % (rad)
    % ltSlk           = 0.39;         % (m)

    % muscleName      = 'TA';
    % fiso            = 905.0;        % (N)
    % lceOpt          = 0.0898;       % (m)
    % alphaOpt        = 0.0873;       % (rad)
    % ltSlk           = 0.2043;       % (m) 

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

filename                = 'patch_info.txt';

fileID                  = fopen(filename, 'w');
fprintf(fileID, 'L_mn0 V_m0\n');

if ~exist('L_mn0_values','var')
    L_mn0_values = [1.0];
end

if ~exist('V_m0_values','var')
    V_m0_values = [0.0];
end

L_mt0_values            = [1.05];

for file_i = 1 : length(L_mn0_values)
    for file_j = 1 : length(V_m0_values)
        fprintf(fileID, '%.1f %.1f\n', L_mn0_values(file_i), V_m0_values(file_j));
    end
end

fclose(fileID);

thisDir  = fileparts(mfilename('fullpath'));
patchDir = fullfile(thisDir, 'patch_result');

if ~exist(patchDir, 'dir')
    mkdir(patchDir);
end

delete(fullfile(patchDir, '*.mat'));

%% Main Calculation

patch_info              = readtable('patch_info.txt');

index_L                 = 1;
index_mt                = 1;
index_V                 = 1;
pnum                    = 1;

L_range                 = 0.050;                                                                     % Muscle length +- range
V_range                 = 0.050;                                                                     % Velocity +- range

aim_Lnum                = length(L_mn0_values);
aim_Ltnum               = length(L_mt0_values);
aim_Vnum                = length(V_m0_values);
aim_pnum                = height(patch_info);

modelConfig.useElasticTendon    = 1;
modelConfig.useFiberDamping     = 0;
modelConfig.damping             = 0.1;
modelConfig.minActivation       = 1e-10;
modelConfig.iterMax             = 10000;
modelConfig.tol                 = 1e-10;
modelConfig.passiveOnlyMode     = true;

patch_result            = cell(aim_pnum, 1000);
visual_result           = cell(aim_pnum, 1);

while (index_L <= aim_Lnum)

    const_b             = 0.0;
    cutpoint            = 1;                                                              % cutpoint = 1 means no data cutting
        
    mod_V_m0            = V_m0_values(index_V);
    mod_L_mn0           = L_mn0_values(index_L);
    mod_L_mt0           = L_mt0_values(index_mt);

    pre_L_mt0           = mod_L_mt0 * (ltSlk + lceOpt * cos(alphaOpt));

    Lf_lower            = 0.0;
    Lf_upper            = 3.0;

    err = 1000;

    while abs(real(err)) > 1e-5
        Lf_mid = 1/2 * (Lf_lower + Lf_upper);

        pathState   = [0; pre_L_mt0];
        muscleState = [0; Lf_mid];

        mtInfo = calcMillard2012DampedEquilibriumMuscleInfo( ...
            1e-10, ... % u = 0
            pathState, ...
            muscleState, ...
            muscleArch, ...
            normMuscleCurves, ...
            modelConfig);
        
        fpe  = mtInfo.muscleDynamicsInfo.fiberForceAlongTendon;
        fse  = mtInfo.muscleDynamicsInfo.tendonForce;

        err = fse - fpe;

        if (err < 0)
            Lf_upper = Lf_mid;
        else
            Lf_lower = Lf_mid;
        end
    end

    Lf0 = mtInfo.muscleLengthInfo.fiberLengthAlongTendon;
    Lt0 = mtInfo.muscleLengthInfo.tendonLength;
    alpha0 = mtInfo.muscleLengthInfo.pennationAngle;

    fprintf("%d\n", mtInfo.muscleDynamicsInfo.activeFiberForce);

    u0                  = 0;  

    pre_dlceAT0         = 0.0;
    pre_lceAT0          = Lf0;

    if ~exist('mass','var')
        mass = 0.0;                                                                 % Mass of the muscle                (kg)
    end

    if ~exist('damp','var')
        ext_damping = 0.0;
    end

    if ~exist('totalTime','var')
        totalTime = 50;
    end
    
    if ~exist('dt','var')
        dt = 0.001;
    end

    time                = 0 : dt : (totalTime - dt);

    if ~exist('steps','var')
        steps = 300;
    end                                          
    
    if ~exist('freqlb','var')
        freqlb = 0.1;                                                                  % Lower bound of the excitation frequency
    end

    if ~exist('freqhb','var')
        freqhb = 100;                                                                  % Upper bound of the excitation frequency
    end
    
    frequencies         = logspace(log10(freqlb), log10(freqhb), steps);        

    freq_len            = length(frequencies);          
    time_len            = length(time);
    
    % [ initialization ]

    sig_in              = zeros(freq_len, time_len);  
    sig_in_a            = zeros(freq_len, time_len);  
    sig_out             = zeros(freq_len, time_len - cutpoint);         

    L_mt_preset         = ones(freq_len, time_len) .* pre_L_mt0;                              % muscle tendon length
    L_mt                = zeros(freq_len, time_len);   
    A_mt                = zeros(freq_len, time_len);                                          % muscle tendon accleration
    V_mt                = zeros(freq_len, time_len);                                          % muscle tendon velocity
    
    fiberForce          = zeros(freq_len, time_len);
    tendonForce         = zeros(freq_len, time_len);
    activefiberForce    = zeros(freq_len, time_len);
    passivefiberForce   = zeros(freq_len, time_len);

    force_diff          = zeros(freq_len, time_len);
    lceAT               = ones(freq_len, time_len) .* pre_lceAT0;
    dlceAT              = zeros(freq_len, time_len);

    fiber_length        = zeros(freq_len, time_len);
    fiber_velocity      = zeros(freq_len, time_len);
    tendon_length       = ones(freq_len, time_len) .* Lt0;
    tendon_velocity     = zeros(freq_len, time_len);
    alpha               = zeros(freq_len, time_len);
    
    idx                 = zeros(1, length(frequencies));

    pha                 = zeros(1, length(frequencies));
    sin_f               = zeros(1, length(frequencies));
    mag                 = zeros(1, length(frequencies));
    exc                 = zeros(1, length(frequencies));

    x_ext               = zeros(freq_len, time_len); %%%
    v_ext               = zeros(freq_len, time_len); %%%
    a_ext               = zeros(freq_len, time_len); %%%
    mass_force          = zeros(freq_len, time_len);
    
    fprintf("-----------------------------\n\n");
    fprintf("[ Main Calculation started ] \n\n");
    fprintf("-----------------------------\n\n");
    
    for k = 1 : length(frequencies)
        
        % Check 'STOP' sign from GUI (MFR_main.m)
        if evalin('base','exist(''MFR_STOP'',''var'') && MFR_STOP')
            builtin('error','Simulation stopped by user via GUI.');
        end
        drawnow limitrate;

        tic;

        fprintf("<strong>[ Patch number: %d ]\n\n</strong>", pnum);
        fprintf("<strong>Initial length   : %f\nInitial Velocity : %f</strong>\n\n", mod_L_mn0, mod_V_m0);
        fprintf("<strong>Step number: %d</strong> \n\n", k);
    
        freq = frequencies(k);
              
        temp_result             = NaN(time_len, 7);
        temp_index              = 1;

        L_mt_preset(k, 1)       = pre_L_mt0;
        V_mt(k, 1)              = 0;
        lceAT(k, 1)             = pre_lceAT0;
        dlceAT(k, 1)            = pre_dlceAT0;

        amp                     = 0.005 * pre_L_mt0;
        L_mt_preset(k, :)       = pre_L_mt0 + amp * sin(2 * pi * freq * time);
        V_mt(k, :)              = amp * 2 * pi * freq * cos(2 * pi * freq * time);
        A_mt(k, :)              = -amp * (2 * pi * freq)^2 * sin(2 * pi * freq * time);
        
        sig_in(k, :)            = L_mt_preset(k, :);

        for i = 1 : time_len - 1                                                            % Explicit method
            
            L_mt(k, i) = L_mt_preset(k, i);
            
            start_  = 0.0; 
            end_    = 1.0; 
            error   = 1000; 

            while (abs(real(error)) > 1e-4)
                
                L_mm = 1/2 * (start_ + end_);
            
                mtInfo_loop = calcMillard2012DampedEquilibriumMuscleInfo( ...
                1e-10, ...
                [0; L_mt(k, i)], ...
                [0; L_mm], ...
                muscleArch, ...
                normMuscleCurves, ...
                modelConfig);
                
                fpe1 = mtInfo_loop.muscleDynamicsInfo.passiveFiberForce;
                ft1  = mtInfo_loop.muscleDynamicsInfo.tendonForce;
                error  = ft1 - fpe1;              % error for force equilbrium (muscle, tendon)
                
                if (error < 0) 
                    end_    = L_mm;
                else  
                    start_  = L_mm; 
                end 

                % if abs(error) < 1e-10 break; end
            end

            mtInfo_Eq = calcMillard2012DampedEquilibriumMuscleInfo( ...
                1e-10, ...
                [0; L_mt(k, i)], ...
                L_mm, ...
                muscleArch, ...
                normMuscleCurves, ...
                modelConfig);

            dlceAT(k, i+1)            = (L_mm - lceAT(k, i)) / dt;
            lceAT(k, i+1)             = L_mm;
            alpha(k, i+1)             = mtInfo_Eq.muscleLengthInfo.pennationAngle;
            tendon_length(k, i+1)     = mtInfo_Eq.muscleLengthInfo.tendonLength;
            tendon_velocity(k, i+1)   = mtInfo_Eq.fiberVelocityInfo.tendonVelocity;

            tendonForce(k, i+1)       = mtInfo_Eq.muscleDynamicsInfo.tendonForce;
            fiberForce(k, i+1)        = mtInfo_Eq.muscleDynamicsInfo.fiberForceAlongTendon;
            activefiberForce(k, i+1)  = mtInfo_Eq.muscleDynamicsInfo.activeFiberForce;
            passivefiberForce(k, i+1) = mtInfo_Eq.muscleDynamicsInfo.passiveFiberForce;

            % a_ext(k, i+1) = (tendonForce(k, i+1) - fse) / mass;
            v_ext(k, i+1) = v_ext(k, i) + dt * a_ext(k, i+1);
            x_ext(k, i+1) = x_ext(k, i) + dt * v_ext(k, i+1);

            mass_force(k, i+1) = fiberForce(k, i+1) + mass * A_mt(k, i+1);
                                  
            temp_result(temp_index, :) = [ lceAT(k, i+1), dlceAT(k, i+1), activefiberForce(k, i+1), passivefiberForce(k, i+1), ...
              mass_force(k, i+1), tendonForce(k, i+1), L_mt(k, i+1) ];

            temp_index = temp_index + 1;
            
        end          

        patch_filename = fullfile(patchDir, sprintf('p%d_k%d.mat', pnum, k));
        save(patch_filename, 'temp_result', '-v7.3');

        loaded = load(patch_filename, 'temp_result');
        temp_result = loaded.temp_result;
   
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

saveFolder1 = 'C:\Users\user\Desktop\MuscleFreqResp-main\MMM test\MMM\src\MMM_result\beforeFFT';

if ~exist(saveFolder1, 'dir')
    mkdir(saveFolder1);
end

for file_idx = 1 : length(L_mn0_values)
    
    fileName = [num2str(L_mn0_values(file_idx)) '_' num2str(u0) '_sol_YB_wod_bF_passive.mat'];
    
    fullPath = fullfile(saveFolder1, fileName);

    savingdata_time = cell(1, freq_len);
    
    for k = 1 : freq_len
        patch_filename = sprintf('patch_result/p%d_k%d.mat', file_idx, k);
        loaded = load(patch_filename, 'temp_result');
        savingdata_time{k} = loaded.temp_result;
    end

    save(fullPath, 'savingdata_time', '-v7.3');
end

saveFolder2 = 'C:\Users\user\Desktop\MuscleFreqResp-main\MMM test\MMM\src\MMM_result\afterFFT';

if ~exist(saveFolder2, 'dir')
    mkdir(saveFolder2);
end

for file_idx = 1 : length(L_mn0_values)
    
    fileName = [num2str(L_mn0_values(file_idx)) '_' num2str(u0) '_sol_YB_wod_aF_passive.mat'];
    
    fullPath = fullfile(saveFolder2, fileName);

    savingdata_freq = visual_result{file_idx};

    save(fullPath, 'savingdata_freq');
end


%% Debugging plot

% target_fnum = 1;
% 
% figure(); plot(time, lceAT(target_fnum, :) ./ cos(alpha(target_fnum, :))); title('Fiber length along Tendon');
% 
% figure(); plot(time, tendon_length(target_fnum, :)); title('Tendon length');
% 
% figure(); plot(time, fiberForce(target_fnum, :)); title('Fiber force');
% 
% figure(); plot(time, tendonForce(target_fnum, :)); title('Tendon force');

%% Bode plots with fitted 2nd order system

% figure();
% sg = sgtitle(sprintf('Bode plot of each conditions w/ u_{0} = %.2f', u0));
% 
% cm = hsv(pnum - 1);
% 
% cutoff_level = -3;
% 
% nat_freqs = zeros(pnum-1, 1);
% 
% for fnum = 1 : pnum - 1
% 
%     bode_freq = visual_result{fnum}{:, 1}; % Frequency (Hz)
%     bode_mag = visual_result{fnum}{:, 2}; % Magnitude (linear)
%     bode_pha = visual_result{fnum}{:, 3}; % Phase (degrees)
% 
%     mag_dB = 20 * log10(bode_mag);
%     omega = 2 * pi * bode_freq; % Frequency (rad/s)
% 
%     subplot(2, 1, 1);
%     semilogx(bode_freq, mag_dB, 'o', 'MarkerSize', 1, ...
%         'Color', cm(fnum, :), 'DisplayName', sprintf('Original: L_{mn} = %.2f, V_{m} = %.2f', L_mn0_values(fnum), V_m0_values(1)));
%     hold on;
%     grid on;
% 
%     subplot(2, 1, 2);
%     semilogx(bode_freq, bode_pha, 'o', 'MarkerSize', 1, ...
%         'Color', cm(fnum, :), 'DisplayName', sprintf('Original: L_{mn} = %.2f, V_{m} = %.2f', L_mn0_values(fnum), V_m0_values(1)));
%     hold on;
%     grid on;
% end
% 
% hold off;
