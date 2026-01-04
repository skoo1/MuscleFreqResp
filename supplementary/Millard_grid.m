clear; clc;

global flag_updateCurves flag_plotNormMuscleCurves
flag_updateCurves = 1;
flag_plotNormMuscleCurves = 0;

thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(thisDir);

addpath(genpath(fullfile(repoRoot, 'MMM test', 'MMM')));
addpath(genpath(fullfile(repoRoot, 'MMM test', 'MMM', 'src')));

assignin('base','flag_updateCurves',flag_updateCurves);
assignin('base','flag_plotNormMuscleCurves',flag_plotNormMuscleCurves);

oldDir = pwd;
cd(thisDir);

muscleAbbr = 'soleus';
tendonStrainAtOneNormForce = 0.049;

fiso   = 3549.0;
lceOpt = 0.05;

Lts    = 0.25;
alpha0 = 0.0;

maximumNormalizedFiberVelocity = 10;

n_len = 101;
n_vel = 101;

normL_grid = linspace(0.5, 1.5, n_len);
normV_grid = linspace(-1.0, 1.0, n_vel);
a_list     = [0.5, 1.0];

normMuscleCurves = createDefaultNormalizedMuscleCurves( ...
    muscleAbbr, tendonStrainAtOneNormForce, flag_updateCurves, 0);

muscleArchitecture = struct();
muscleArchitecture.name  = muscleAbbr;
muscleArchitecture.abbr  = muscleAbbr;

muscleArchitecture.fiso  = fiso;
muscleArchitecture.optimalFiberLength = lceOpt;

muscleArchitecture.maximumNormalizedFiberVelocity = maximumNormalizedFiberVelocity;
muscleArchitecture.pennationAngle    = alpha0;
muscleArchitecture.tendonSlackLength = Lts;

muscleArchitecture.minimumFiberLength = 0.01;
muscleArchitecture.minimumFiberLengthAlongTendon = 0.01;
muscleArchitecture.pennationAngleAtMinumumFiberLength = alpha0;

modelConfig = struct();
modelConfig.useElasticTendon = 0;
modelConfig.useFiberDamping  = 0;
modelConfig.damping          = 0.1;
modelConfig.iterMax          = 100;
modelConfig.tol              = 1e-10;
modelConfig.minActivation    = 1e-8;
modelConfig.passiveOnlyMode  = 0;

nRows = numel(a_list) * n_len * n_vel;

activation          = zeros(nRows,1);
norm_fiber_length   = zeros(nRows,1);
norm_fiber_velocity = zeros(nRows,1);

activeFiberForce_N  = zeros(nRows,1);
passiveFiberForce_N = zeros(nRows,1);
totalFiberForce_N   = zeros(nRows,1);

activeAT_N  = zeros(nRows,1);
passiveAT_N = zeros(nRows,1);
totalAT_N   = zeros(nRows,1);

% --- Normalized versions (force / fiso)
activeAT_norm  = zeros(nRows,1);
passiveAT_norm = zeros(nRows,1);
totalAT_norm   = zeros(nRows,1);

row = 0;

for a = a_list
    for i = 1:n_len
        lceN = normL_grid(i);
        lce  = lceN * lceOpt;

        for j = 1:n_vel
            row = row + 1;
            dlceN = normV_grid(j);
            dlce  = dlceN * (lceOpt * maximumNormalizedFiberVelocity);

            fkAT = calcFixedWidthPennatedFiberKinematicsAlongTendon( ...
                lce, dlce, lceOpt, alpha0);

            lceAT  = fkAT.fiberLengthAlongTendon;
            dlceAT = fkAT.fiberVelocityAlongTendon;

            lp  = lceAT + Lts;
            dlp = dlceAT;

            activationState = a;
            pathState = [dlp; lp];
            muscleState = [];

            mtInfo = calcMillard2012DampedEquilibriumMuscleInfo( ...
                activationState, pathState, muscleState, ...
                muscleArchitecture, normMuscleCurves, modelConfig);

            activation(row)          = a;
            norm_fiber_length(row)   = lceN;
            norm_fiber_velocity(row) = dlceN;

            activeFiberForce_N(row)  = mtInfo.muscleDynamicsInfo.activeFiberForce;
            passiveFiberForce_N(row) = mtInfo.muscleDynamicsInfo.passiveFiberForce;
            totalFiberForce_N(row)   = mtInfo.muscleDynamicsInfo.fiberForce;

            activeAT_N(row)  = mtInfo.muscleDynamicsInfo.activeFiberForce * mtInfo.muscleLengthInfo.cosPennationAngle;
            passiveAT_N(row) = mtInfo.muscleDynamicsInfo.passiveFiberForce * mtInfo.muscleLengthInfo.cosPennationAngle;
            totalAT_N(row)   = mtInfo.muscleDynamicsInfo.fiberForceAlongTendon;

            activeAT_norm(row)  = activeAT_N(row)  / fiso;
            passiveAT_norm(row) = passiveAT_N(row) / fiso;
            totalAT_norm(row)   = totalAT_N(row)   / fiso;
        end
    end
end

T = table(activation, norm_fiber_length, norm_fiber_velocity, ...
          activeAT_N, passiveAT_N, totalAT_N, ...
          activeAT_norm, passiveAT_norm, totalAT_norm, ...
          activeFiberForce_N, passiveFiberForce_N, totalFiberForce_N);

outCsv = fullfile(thisDir, 'Millard_grid.csv');
writetable(T, outCsv);

fprintf('[OK] Wrote: %s\n', outCsv);
fprintf('[OK] Rows : %d\n', height(T));

cd(oldDir);
