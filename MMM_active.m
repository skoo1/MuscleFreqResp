% By Minseung Kim and Seungbum Koo
% January 26, 2026

clc; clear;
addpath('.\MMM\src\');

% Sim configuration
% L_mtn_intput = 1.05; % for passive test
L_mn_input     = 1.0;
V_m_input      = 0.0;
U_input        = 0.5;
Amp_input      = 0.01;
Mass_input     = 3e9;
Damping_imput  = 0.0;
SimTime_input  = 120;
SimDt_input    = 0.001;
FreqLow_input  = 0.1;
FreqHigh_input = 100;
NumFreqSamples = 100;

% Muscle properties
muscleName      = 'Soleus';
F_mo            = 3549.0;       % (N)
L_mo            = 0.05;         % (m)
L_ts            = 0.25;         % (m)
alphaOpt        = 0.4363;       % (rad)
V_mmax_norm     = 10;

% External dynamics parameters
mass = Mass_input;
damp = Damping_imput;

% Operating point
L_mn_target = L_mn_input; 
V_m_target  = V_m_input;
u0          = U_input;

% Load Millard-Muscle model parameters
[muscleArch, normMuscleCurves, modelConfig, MMM] = ...
    init_MMM(muscleName, muscleName, F_mo, L_mo, L_ts, alphaOpt, V_mmax_norm);

% Muscle parameter initialization
L_mn0           = L_mn_target;
L_m0            = L_mn_target * L_mo;
L_mt0           = L_m0 * cos(alphaOpt) + L_ts;
V_mt0           = V_m_target;
V_m0            = V_m_target;
L_m_AT0         = L_m0 * cos(alphaOpt);

% Calculate the equilibrium external force at initial state
F_ext_equil = findExternalForceForFiberLength( ...
    L_m_AT0, L_ts, u0, muscleArch, normMuscleCurves, modelConfig );

% Dynamics simulation parameters
totalTime  = SimTime_input;
dt         = SimDt_input;
time_array = 0 : dt : (totalTime - dt);
time_len   = length(time_array);
steps      = NumFreqSamples;

% Frequency sampling parameters
freqlb      = FreqLow_input;
freqhb      = FreqHigh_input;
frequencies = logspace(log10(freqlb), log10(freqhb), steps);                                              
freq_len    = length(frequencies);          

% Storage initialization
visual_result       = cell(1, 4);
pha                 = zeros(1, length(frequencies));
sin_f               = zeros(1, length(frequencies));
mag                 = zeros(1, length(frequencies));
exc                 = zeros(1, length(frequencies));

parfor k = 1 : length(frequencies)       
    tic;
    freq = frequencies(k);

    F_m_AT = zeros(1, time_len);

    amp  = Amp_input;
    u    = sinwave(freq, time_array, u0, u0, amp);
    u(u < 0) = 0; 
    u(u > 1) = 1;

    a    = ones(1, length(u)) * u0;
    a(1) = u(1);
    
    L_mt      = L_mt0;
    V_mt      = 0;
    L_m_AT    = L_m_AT0;
    
    % Calculate the force measured at the sensor or F_m_AT(i)
    % as a response to the input signal u(i)
    for i = 1 : time_len
        if (i>=2)
            da_dt = (u(i) - a(i-1)) / get_Tau(a(i-1), u(i), MMM);
            a(i) = a(i-1) + dt * da_dt;
        end

        pathState        = [V_mt; L_mt];
        muscleState      = L_m_AT;
        mtInfo           = calcMillard2012DampedEquilibriumMuscleInfo(...
                                a(i), pathState, muscleState, ...
                                muscleArch, normMuscleCurves, modelConfig);
        
        F_t             = mtInfo.muscleDynamicsInfo.tendonForce;
        F_m_AT(i)        = mtInfo.muscleDynamicsInfo.fiberForceAlongTendon;

        V_m_AT          = mtInfo.state.derivative;
        L_m_AT          = L_m_AT + V_m_AT * dt;

        A_mt    = (F_ext_equil - F_t - damp * V_mt) / mass;
        V_mt    = V_mt + A_mt * dt;
        L_mt    = L_mt + V_mt * dt;
    end

    valid_range = round(time_len * 0.1) : time_len-1;
    sig_in  = u(valid_range);
    sig_out = F_m_AT(valid_range);
    
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
% sol/gast (muscle type) | YB/OB (aging) | wod/d (damping)
% Format:: MMM_results_KMS_Lmn0_uo_sol_YB_wod.mat %

if exist('mass','var')
    if abs(mass - 3e8) < 1e-6     % 300000000 kg → isometric
        mass_label = 'isometric';
    else
        mass_label = [num2str(mass) 'kg'];
    end
else
    mass_label = 'unknownMass';
end

saveFolder = '.\MMM_results';
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


function F_ext_equil = findExternalForceForFiberLength( ...
    L_m_AT, L_ts, u0, muscleArch, normMuscleCurves, modelConfig )

    L_mt  = L_m_AT + L_ts; 
    pathState   = [0, L_mt];
    muscleState = [0, L_m_AT];

    mtInfo = calcMillard2012DampedEquilibriumMuscleInfo( ...
        u0, pathState, muscleState, ...
        muscleArch, normMuscleCurves, modelConfig);

    currentFiberForce = mtInfo.muscleDynamicsInfo.fiberForceAlongTendon;
    F_ext_equil       = currentFiberForce;
end


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
    modelConfig.useFiberDamping     = 0;  
    modelConfig.damping             = 0.1;
    modelConfig.minActivation       = 1e-10;
    modelConfig.iterMax             = 10000;
    modelConfig.tol                 = 1e-10;

    % 6. Additional Constants (MMM)
    MMM = struct();
    MMM.Tau_a = 0.01;
    MMM.Tau_d = 0.04;
end
