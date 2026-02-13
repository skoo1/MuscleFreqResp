% main.m - Run Analysis for All Simulations
clear; clc; close all;

% 공통 설정
L_mn_target = 1.0; 

%% 1. MMM Active Simulation
% 폴더명은 C++ 코드에서 생성한 이름과 일치해야 합니다.
try
    csv_folder  = './MMM_Active_results_csv/';
    save_folder = './MMM_results/';
    u0          = 0.5;
    analyze_simulation_results(csv_folder, save_folder, 'MMM', 'Active', u0, L_mn_target);
catch
    fprintf('MMM_Active Error!\n')
end

%% 2. MMM Passive Simulation
try
    csv_folder  = './MMM_Passive_results_csv/';
    save_folder = './MMM_results/';
    u0          = 0.0; % Passive는 u=0
    analyze_simulation_results(csv_folder, save_folder, 'MMM', 'Passive', u0, L_mn_target);
catch
    fprintf('MMM_Passive Error!\n')
end

%% 3. TMM Active Simulation
try
    csv_folder  = './TMM_Active_results_csv/';
    save_folder = './TMM_results/';
    u0          = 0.5;
    analyze_simulation_results(csv_folder, save_folder, 'TMM', 'Active', u0, L_mn_target);
catch
    fprintf('TMM_Active Error!\n')
end

%% 4. TMM Passive Simulation
try
    csv_folder  = './TMM_Passive_results_csv/';
    save_folder = './TMM_results/';
    u0          = 0.0;
    analyze_simulation_results(csv_folder, save_folder, 'TMM', 'Passive', u0, L_mn_target);
catch
    fprintf('TMM_Passive Error!\n')
end

fprintf('\nAll analyses completed successfully.\n');