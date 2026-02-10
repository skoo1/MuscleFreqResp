% By Minseung Kim and Seungbum Koo
% KAIST, Daejeon, South Korea
% February 8, 2026

clear; clc; close all;

%% Thelen Muscle Model
fprintf('--- Testing Thelen Model ---\n');

% Thelen Muscle Instantiation
% ThelenMuscle(Muscle_name, F_iso, L_m_opt, L_t_slack, Alpha_opt, Mass_not_used, Damping_not_used)
thelen = ThelenMuscle('Soleus', 3549, 0.05, 0.25, 0.4363, 0.0, 0.0, "None");

% Thelen Muscle Tests
results_thelen_active  = run_active_test(thelen);
results_thelen_passive = run_passive_test(thelen);

%% Millard Muscle Model
fprintf('--- Testing Millard Models ---\n');

% Include Millard Muscle Library
% % https://github.com/mjhmilla/Millard2012EquilibriumMuscleMatlabPort
addpath('./MMM/src/');

fprintf('--- Testing Millard Classic Model ---\n');
% MillardMuscle(Muscle_name, F_iso, L_m_opt, L_t_slack, Alpha_opt, Mass_not_used, Damping_not_used)
millard_Classic = MillardMuscle('Soleus', 3549, 0.05, 0.25, 0.4363, 0.0, 0.0, "Classic");

% Millard Muscle Tests
results_millard_Classic_active  = run_active_test(millard_Classic);
results_millard_Classic_passive = run_passive_test(millard_Classic);

fprintf('--- Testing Millard DEq Model ---\n');
millard_DEq = MillardMuscle('Soleus', 3549, 0.05, 0.25, 0.4363, 0.0, 0.0, "DEq");

% Millard Muscle Tests
results_millard_DEq_active  = run_active_test(millard_DEq);
results_millard_DEq_passive = run_passive_test(millard_DEq);

fprintf('--- Testing Millard Rigid Model ---\n');
millard_Rigid = MillardMuscle('Soleus', 3549, 0.05, 0.25, 0.4363, 0.0, 0.0, "Rigid");

% Millard Muscle Tests
results_millard_Rigid_active  = run_active_test(millard_Rigid);
results_millard_Rigid_passive = run_passive_test(millard_Rigid);


%% Results comparison
% Active tests ======
figure('Name', 'Model Comparison', 'Color', 'w');
% Magnitude Plot
subplot(2,1,1);
semilogx(results_thelen_active{1}, 20*log10(results_thelen_active{2}), 'r-', 'LineWidth', 1.5); hold on;
semilogx(results_millard_Classic_active{1}, 20*log10(results_millard_Classic_active{2}), 'g-', 'LineWidth', 1.5);
semilogx(results_millard_DEq_active{1}, 20*log10(results_millard_DEq_active{2}), 'b-', 'LineWidth', 1.5);
semilogx(results_millard_Rigid_active{1}, 20*log10(results_millard_Rigid_active{2}), 'y-', 'LineWidth', 1.5);
title('Active Response: Magnitude');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
legend('Thelen', 'Millard Classic', 'Millard DEq', 'Millard Rigid'); grid on; xlim([0.1 100]); ylim([10 100]);

% Phase Plot
subplot(2,1,2);
semilogx(results_thelen_active{1}, results_thelen_active{3}, 'r-', 'LineWidth', 1.5); hold on;
semilogx(results_millard_Classic_active{1}, results_millard_Classic_active{3}, 'g-', 'LineWidth', 1.5);
semilogx(results_millard_DEq_active{1}, results_millard_DEq_active{3}, 'b-', 'LineWidth', 1.5);
semilogx(results_millard_Rigid_active{1}, results_millard_Rigid_active{3}, 'y-', 'LineWidth', 1.5);

title('Active Response: Phase');
xlabel('Frequency (Hz)'); ylabel('Phase (deg)');
grid on; xlim([0.1 100]); ylim([-180 10]);
% Passive tests ======
figure('Name', 'Model Comparison', 'Color', 'w');
% Magnitude Plot
subplot(2,1,1);
semilogx(results_thelen_passive{1}, 20*log10(results_thelen_passive{2}), 'r-', 'LineWidth', 1.5); hold on;
semilogx(results_millard_Classic_passive{1}, 20*log10(results_millard_Classic_passive{2}), 'g-', 'LineWidth', 1.5);
semilogx(results_millard_DEq_passive{1}, 20*log10(results_millard_DEq_passive{2}), 'b-', 'LineWidth', 1.5);
semilogx(results_millard_Rigid_passive{1}, 20*log10(results_millard_Rigid_passive{2}), 'y-', 'LineWidth', 1.5);
title('Passive Response: Magnitude');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
legend('Thelen', 'Millard Classic', 'Millard DEq', 'Millard Rigid'); grid on; xlim([0.1 100]); ylim([10 100]);

% Phase Plot
subplot(2,1,2);
semilogx(results_thelen_passive{1}, results_thelen_passive{3}, 'r-', 'LineWidth', 1.5); hold on;
semilogx(results_millard_Classic_passive{1}, results_millard_Classic_passive{3}, 'g-', 'LineWidth', 1.5);
semilogx(results_millard_DEq_passive{1}, results_millard_DEq_passive{3}, 'b-', 'LineWidth', 1.5);
semilogx(results_millard_Rigid_passive{1}, results_millard_Rigid_passive{3}, 'y-', 'LineWidth', 1.5);
title('Passive Response: Phase');
xlabel('Frequency (Hz)'); ylabel('Phase (deg)');
grid on; xlim([0.1 100]); ylim([-180 10]);