% By Minseung Kim, Seungwoo Yoon and Seungbum Koo
% KAIST, Daejeon, South Korea
% February 23, 2026

% main.m - Run Analysis for All Simulations
clear; clc; close all;

% L_mn_target = 1.0; 

%% 1. MMM Active Simulation
try
    csv_folder  = './MMM_Classic_Active_csv_0.600000/';
    save_folder = './MMM_Classic_Active_result_0.600000/';
    u0          = 0.5;
    analyze_simulation_results(csv_folder, save_folder, 'MMM-Classic', 'Active', u0, 0.6);
catch
    fprintf('MMM_Classic_Active_0.600000 Error!\n')
end

try
    csv_folder  = './MMM_Classic_Active_csv_0.800000/';
    save_folder = './MMM_Classic_Active_result_0.800000/';
    u0          = 0.5;
    analyze_simulation_results(csv_folder, save_folder, 'MMM-Classic', 'Active', u0, 0.8);
catch
    fprintf('MMM_Classic_Active_0.800000 Error!\n')
end

try
    csv_folder  = './MMM_Classic_Active_csv_1.000000/';
    save_folder = './MMM_Classic_Active_result_1.000000/';
    u0          = 0.5;
    analyze_simulation_results(csv_folder, save_folder, 'MMM-Classic', 'Active', u0, 1.0);
catch
    fprintf('MMM_Classic_Active_1.000000 Error!\n')
end

try
    csv_folder  = './MMM_Classic_Active_csv_1.200000/';
    save_folder = './MMM_Classic_Active_result_1.200000/';
    u0          = 0.5;
    analyze_simulation_results(csv_folder, save_folder, 'MMM-Classic', 'Active', u0, 1.2);
catch
    fprintf('MMM_Classic_Active_1.200000 Error!\n')
end

try
    csv_folder  = './MMM_Classic_Active_csv_1.400000/';
    save_folder = './MMM_Classic_Active_result_1.400000/';
    u0          = 0.5;
    analyze_simulation_results(csv_folder, save_folder, 'MMM-Classic', 'Active', u0, 1.4);
catch
    fprintf('MMM_Classic_Active_1.400000 Error!\n')
end

try
    csv_folder  = './MMM_DEq_Active_csv_0.600000/';
    save_folder = './MMM_DEq_Active_result_0.600000/';
    u0          = 0.5;
    analyze_simulation_results(csv_folder, save_folder, 'MMM-DEq', 'Active', u0, 0.6);
catch
    fprintf('MMM_DEq_Active_0.600000 Error!\n')
end

try
    csv_folder  = './MMM_DEq_Active_csv_0.800000/';
    save_folder = './MMM_DEq_Active_result_0.800000/';
    u0          = 0.5;
    analyze_simulation_results(csv_folder, save_folder, 'MMM-DEq', 'Active', u0, 0.8);
catch
    fprintf('MMM_DEq_Active_0.800000 Error!\n')
end

try
    csv_folder  = './MMM_DEq_Active_csv_1.000000/';
    save_folder = './MMM_DEq_Active_result_1.000000/';
    u0          = 0.5;
    analyze_simulation_results(csv_folder, save_folder, 'MMM-DEq', 'Active', u0, 1.0);
catch
    fprintf('MMM_DEq_Active_1.000000 Error!\n')
end

try
    csv_folder  = './MMM_DEq_Active_csv_1.200000/';
    save_folder = './MMM_DEq_Active_result_1.200000/';
    u0          = 0.5;
    analyze_simulation_results(csv_folder, save_folder, 'MMM-DEq', 'Active', u0, 1.2);
catch
    fprintf('MMM_DEq_Active_1.200000 Error!\n')
end

try
    csv_folder  = './MMM_DEq_Active_csv_1.400000/';
    save_folder = './MMM_DEq_Active_result_1.400000/';
    u0          = 0.5;
    analyze_simulation_results(csv_folder, save_folder, 'MMM-DEq', 'Active', u0, 1.4);
catch
    fprintf('MMM_DEq_Active_1.400000 Error!\n')
end

try
    csv_folder  = './MMM_Rigid_Active_csv_0.600000/';
    save_folder = './MMM_Rigid_Active_result_0.600000/';
    u0          = 0.5;
    analyze_simulation_results(csv_folder, save_folder, 'MMM-Rigid', 'Active', u0, 0.6);
catch
    fprintf('MMM_Rigid_Active_0.600000 Error!\n')
end

try
    csv_folder  = './MMM_Rigid_Active_csv_0.800000/';
    save_folder = './MMM_Rigid_Active_result_0.800000/';
    u0          = 0.5;
    analyze_simulation_results(csv_folder, save_folder, 'MMM-Rigid', 'Active', u0, 0.8);
catch
    fprintf('MMM_Rigid_Active_0.800000 Error!\n')
end

try
    csv_folder  = './MMM_Rigid_Active_csv_1.000000/';
    save_folder = './MMM_Rigid_Active_result_1.000000/';
    u0          = 0.5;
    analyze_simulation_results(csv_folder, save_folder, 'MMM-Rigid', 'Active', u0, 1.0);
catch
    fprintf('MMM_Rigid_Active_1.000000 Error!\n')
end

try
    csv_folder  = './MMM_Rigid_Active_csv_1.200000/';
    save_folder = './MMM_Rigid_Active_result_1.200000/';
    u0          = 0.5;
    analyze_simulation_results(csv_folder, save_folder, 'MMM-Rigid', 'Active', u0, 1.2);
catch
    fprintf('MMM_Rigid_Active_1.200000 Error!\n')
end

try
    csv_folder  = './MMM_Rigid_Active_csv_1.400000/';
    save_folder = './MMM_Rigid_Active_result_1.400000/';
    u0          = 0.5;
    analyze_simulation_results(csv_folder, save_folder, 'MMM-Rigid', 'Active', u0, 1.4);
catch
    fprintf('MMM_Rigid_Active_1.400000 Error!\n')
end

%% 2. MMM Passive Simulation
try
    csv_folder  = './MMM_Classic_Passive_csv/';
    save_folder = './MMM_Classic_Passive_result/';
    u0          = 0.0;
    analyze_simulation_results(csv_folder, save_folder, 'MMM-Classic', 'Passive', u0, 1.05);
catch
    fprintf('MMM_Classic_Passive Error!\n')
end

try
    csv_folder  = './MMM_DEq_Passive_csv/';
    save_folder = './MMM_DEq_Passive_result/';
    u0          = 0.0;
    analyze_simulation_results(csv_folder, save_folder, 'MMM-DEq', 'Passive', u0, 1.05);
catch
    fprintf('MMM_DEq_Passive Error!\n')
end

try
    csv_folder  = './MMM_Rigid_Passive_csv/';
    save_folder = './MMM_Rigid_Passive_result/';
    u0          = 0.0;
    analyze_simulation_results(csv_folder, save_folder, 'MMM-Rigid', 'Passive', u0, 1.05);
catch
    fprintf('MMM_Rigid_Passive Error!\n')
end


%% 3. TMM Active Simulation
try
    csv_folder  = './TMM_Active_result_csv_0.600000/';
    save_folder = './TMM_Active_result_0.600000/';
    u0          = 0.5;
    analyze_simulation_results(csv_folder, save_folder, 'TMM', 'Active', u0, 0.6);
catch
    fprintf('TMM_Active_0.600000 Error!\n')
end

try
    csv_folder  = './TMM_Active_result_csv_0.800000/';
    save_folder = './TMM_Active_result_0.800000/';
    u0          = 0.5;
    analyze_simulation_results(csv_folder, save_folder, 'TMM', 'Active', u0, 0.8);
catch
    fprintf('TMM_Active_0.800000 Error!\n')
end

try
    csv_folder  = './TMM_Active_result_csv_1.000000/';
    save_folder = './TMM_Active_result_1.000000/';
    u0          = 0.5;
    analyze_simulation_results(csv_folder, save_folder, 'TMM', 'Active', u0, 1.0);
catch
    fprintf('TMM_Active_1.000000 Error!\n')
end

try
    csv_folder  = './TMM_Active_result_csv_1.200000/';
    save_folder = './TMM_Active_result_1.200000/';
    u0          = 0.5;
    analyze_simulation_results(csv_folder, save_folder, 'TMM', 'Active', u0, 1.2);
catch
    fprintf('TMM_Active_1.200000 Error!\n')
end

try
    csv_folder  = './TMM_Active_result_csv_1.400000/';
    save_folder = './TMM_Active_result_1.400000/';
    u0          = 0.5;
    analyze_simulation_results(csv_folder, save_folder, 'TMM', 'Active', u0, 1.4);
catch
    fprintf('TMM_Active_1.400000 Error!\n')
end

%% 4. TMM Passive Simulation
try
    csv_folder  = './TMM_Passive_result_csv/';
    save_folder = './TMM_Passive_result/';
    u0          = 0.0;
    analyze_simulation_results(csv_folder, save_folder, 'TMM', 'Passive', u0, 1.05);
catch
    fprintf('TMM_Passive Error!\n')
end

fprintf('\nAll analyses completed successfully.\n');