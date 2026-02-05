% By Minseung Kim and Seungbum Koo
% February 5, 2026

clc; clear;
addpath('.\MMM\src\');

% Sim configuration
L_mtn_input    = 1.05;
% L_mn_input   = 1.0; % for passive test
% V_m_input    = 0.0; % for passive test
U_input        = 0.0; % 1.49012e-08 Minimum in the Millard Muscle Library
Amp_input      = 0.005;
Mass_input     = 0.0;
Damping_input  = 0.0;
SimTime_input  = 100;
SimDt_input    = 0.001;
FreqLow_input  = 0.1;
FreqHigh_input = 100;
NumFreqSamples = 1000;

% Muscle properties
muscleName      = 'Soleus';
F_mo            = 3549.0;       % (N)
L_mo            = 0.05;         % (m)
L_ts            = 0.25;         % (m)
alphaOpt        = 0.4363;       % (rad)
V_mmax_norm     = 10;

% External dynamics parameters
mass         = Mass_input;
damp         = Damping_input;

% Operating point
L_mtn_target    = L_mtn_input;

% Load Millard-Muscle model parameters
[muscleArch, normMuscleCurves, modelConfig, MMM] = ...
    init_MMM(muscleName, muscleName, F_mo, L_mo, L_ts, alphaOpt, V_mmax_norm);

L_mt0               = L_ts + L_mo * cos(alphaOpt);

low_L_m_AT          = 0.0;
high_L_m_AT         = L_mt0;
error1              = 1000;

while abs(real(error1)) > 1e-8
    L_m_AT = 1/2 * (low_L_m_AT + high_L_m_AT);

    pathState   = [0; L_mt0 * L_mtn_target];
    muscleState = [0; L_m_AT];

    mtInfo = calcMillard2012DampedEquilibriumMuscleInfo( ...
        1e-10, pathState, muscleState, ...
        muscleArch, normMuscleCurves, modelConfig);
    
    F_m_AT = mtInfo.muscleDynamicsInfo.fiberForceAlongTendon;
    F_t    = mtInfo.muscleDynamicsInfo.tendonForce;

    error1 = F_t - F_m_AT;

    if (error1 < 0)
        high_L_m_AT = L_m_AT;
    else
        low_L_m_AT = L_m_AT;
    end
end

% Operating point
L_m_AT0      = L_m_AT;
u0           = U_input;  

totalTime    = SimTime_input;
dt           = SimDt_input;
time_array   = 0 : dt : (totalTime - dt);
time_len     = length(time_array);
steps        = NumFreqSamples;

freqlb       = FreqLow_input;
freqhb       = FreqHigh_input;
frequencies  = logspace(log10(freqlb), log10(freqhb), steps);        
freq_len     = length(frequencies);          

% Storage initialization
visual_result = cell(1, 4);
pha           = zeros(1, length(frequencies));
sin_f         = zeros(1, length(frequencies));
mag           = zeros(1, length(frequencies));
exc           = zeros(1, length(frequencies));

parfor k = 1 : length(frequencies)
    tic;
    F_m_AT    = zeros(1, time_len);

    L_m_AT    = L_m_AT0;
    freq      = frequencies(k);

    amp               = Amp_input * L_mt0 * L_mtn_target;
    L_mt_preset       = L_mt0 * L_mtn_target + amp * sin(2 * pi * freq * time_array);
    V_mt              = amp * 2 * pi * freq * cos(2 * pi * freq * time_array);
    A_mt              = -amp * (2 * pi * freq)^2 * sin(2 * pi * freq * time_array);
    sig_in            = L_mt_preset;

    for i = 1 : time_len
        L_mt = L_mt_preset(i);

        count    = 0;
        max_iter = 50;
        error1   = 1.0;
        delta    = 1e-7;

        while (abs(real(error1)) > 1e-8 && count < max_iter)
            count = count + 1;
        
            mtInfo1 = calcMillard2012DampedEquilibriumMuscleInfo( ...
                1e-10, [0; L_mt], [0; L_m_AT], ...
                muscleArch, normMuscleCurves, modelConfig);
            
            F_m_AT1 = mtInfo1.muscleDynamicsInfo.fiberForceAlongTendon;
            F_t1    = mtInfo1.muscleDynamicsInfo.tendonForce;
            error1  = F_t1 - F_m_AT1;        % error for force equilbrium (muscle, tendon)
            
            L_m_AT_perturb = L_m_AT + delta;
            mtInfo2 = calcMillard2012DampedEquilibriumMuscleInfo( ...
                1e-10, [0; L_mt], [0; L_m_AT_perturb], ...
                muscleArch, normMuscleCurves, modelConfig);

            F_m_AT2 = mtInfo2.muscleDynamicsInfo.fiberForceAlongTendon;
            F_t2    = mtInfo2.muscleDynamicsInfo.tendonForce;
            error2  = F_t2 - F_m_AT2;

            J       = (error2 - error1) / delta;
            if abs(J) < 1e-14
                error('Cannot update: stiffness is close to zero.');
            end

            L_m_AT  = L_m_AT - (error1 / J);
            if (L_m_AT < 0) || (L_m_AT > L_mt0)
                error('Newton solver failed: Solution %.4e is out of bounds [%.4e, %.4e].', ...
                    L_m_AT, low_L_m_AT, high_L_m_AT);
            end
        end

        mtInfo_Eq   = calcMillard2012DampedEquilibriumMuscleInfo( ...
            1e-10, [0; L_mt], L_m_AT, ...
            muscleArch, normMuscleCurves, modelConfig);

        F_m_AT(i)   = mtInfo_Eq.muscleDynamicsInfo.fiberForceAlongTendon;
    end          

    sig_out   = F_m_AT;

    X_fft           = fft(sig_in); 
    F_fft           = fft(sig_out);
    
    Fs              = 1 / dt;
    n_x             = length(sig_in);
    f               = (0 : n_x - 1) * (Fs / n_x);

    [~, idx]        = min(abs(f - freq));

    sin_f(k)        = freq;
    mag  (k)        = abs(F_fft(idx) / X_fft(idx));
    pha  (k)        = angle(F_fft(idx)) - angle(X_fft(idx));
    exc  (k)        = abs(X_fft(idx));

    fprintf('%d, %f, %f, %f, %.1f\n', k, freq, mag(k), pha(k), toc);

end

pha = unwrap(pha) * (180 / pi);

visual_result{1} = sin_f;
visual_result{2} = mag;
visual_result{3} = pha;
visual_result{4} = exc;

% Extracting simulation data into .mat file
% Need modify filename to adjust the simulation setting
% sol/gast (muscle type) | YB/OB (aging) | wod/d (damping)
% Format:: MMM_results_KMS_Lmn0_uo_sol_YB_wod.mat %

saveFolder = '.\MMM_results';
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end

fileName = [num2str(L_mtn_target) '_' num2str(u0) '_' muscleName '_YB_wod_aF_DEq_passive.mat'];
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
function [muscleArch, normMuscleCurves, modelConfig, MMM] = init_MMM(muscleName, muscleAbbr, F_mo, L_mo, L_ts, alphaOpt, V_mmax_norm)

    % 1. Default Parameters & Constants
    maximumPennationAngle = 89 * (pi/180); 
    
    % If this is greater than 0 this value will be used to make 
    % the tendon-force-length curve. Otherwise the default of 0.049 is taken.
    tendonStrainAtOneNormForce = 0.033; 
    flag_plotNormMuscleCurves = 0;
    flag_updateCurves = 1;

    % 2. Normalized Muscle Curves Creation (Millard Model)
    normMuscleCurves = createDefaultNormalizedMuscleCurves(muscleName,...
                                            tendonStrainAtOneNormForce,...
                                            flag_updateCurves,...
                                            flag_plotNormMuscleCurves);

    % 3. Muscle Architecture Setup (muscleArch)
    muscleArch = struct();
    muscleArch.name = muscleName;
    muscleArch.abbr = muscleAbbr;
    muscleArch.fiso = F_mo;
    muscleArch.optimalFiberLength = L_mo;
    muscleArch.maximumNormalizedFiberVelocity = V_mmax_norm;
    muscleArch.pennationAngle = alphaOpt;
    muscleArch.tendonSlackLength = L_ts;

    % 4. Kinematics Limit Calculation
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

    % 5. Model Configuration (modelConfig)
    modelConfig = struct();
    modelConfig.useElasticTendon    = 1;
    modelConfig.useFiberDamping     = 1;  
    modelConfig.damping             = 0.1;
    modelConfig.minActivation       = 1.49012e-08;  % Minimum in the Millard Muscle Library
    modelConfig.iterMax             = 10000;
    modelConfig.tol                 = 1e-10;

    % 6. Additional Constants (MMM)
    MMM = struct();
    MMM.Tau_a = 0.01;
    MMM.Tau_d = 0.04;
end