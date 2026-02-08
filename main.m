% By Minseung Kim and Seungbum Koo
% KAIST, Daejeon, South Korea
% February 8, 2026

clear; clc; close all;

%% Thelen Muscle Model
fprintf('--- Testing Thelen Model ---\n');

% Thelen Muscle Instantiation
% ThelenMuscle(Muscle_name, F_iso, L_m_opt, L_t_slack, Alpha_opt, Mass, Damping)
thelen = ThelenMuscle('Soleus', 3549, 0.05, 0.25, 0.4363, 1, 0.0);

% Thelen Muscle Tests
results_thelen_active  = run_active_test(thelen);
results_thelen_passive = run_passive_test(thelen);

%% Millard Muscle Model
fprintf('--- Testing Millard Model ---\n');

% Include Millard Muscle Library
% % https://github.com/mjhmilla/Millard2012EquilibriumMuscleMatlabPort
addpath('./MMM/src/');

% Millard Muscle Instantiation
% MillardMuscle(Muscle_name, F_iso, L_m_opt, L_t_slack, Alpha_opt, Mass, Damping)
millard = MillardMuscle('Soleus', 3549, 0.05, 0.25, 0.4363, 1, 0.0);

% Millard Muscle Tests
results_millard_active  = run_active_test(millard);
results_millard_passive = run_passive_test(millard);


%% Results comparison
% Active tests ======
figure('Name', 'Model Comparison', 'Color', 'w');
% Magnitude Plot
subplot(2,1,1);
semilogx(results_thelen_active{1}, 20*log10(results_thelen_active{2}), 'r-', 'LineWidth', 1.5); hold on;
semilogx(results_millard_active{1}, 20*log10(results_millard_active{2}), 'b--', 'LineWidth', 1.5);
title('Active Response: Magnitude');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
legend('Thelen', 'Millard'); grid on; xlim([0.1 100]); ylim([10 100]);

% Phase Plot
subplot(2,1,2);
semilogx(results_thelen_active{1}, results_thelen_active{3}, 'r-', 'LineWidth', 1.5); hold on;
semilogx(results_millard_active{1}, results_millard_active{3}, 'b--', 'LineWidth', 1.5);
title('Active Response: Phase');
xlabel('Frequency (Hz)'); ylabel('Phase (deg)');
grid on; xlim([0.1 100]); ylim([-180 10]);
% Passive tests ======
figure('Name', 'Model Comparison', 'Color', 'w');
% Magnitude Plot
subplot(2,1,1);
semilogx(results_thelen_passive{1}, 20*log10(results_thelen_passive{2}), 'r-', 'LineWidth', 1.5); hold on;
semilogx(results_millard_passive{1}, 20*log10(results_millard_passive{2}), 'b--', 'LineWidth', 1.5);
title('Passive Response: Magnitude');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
legend('Thelen', 'Millard'); grid on; xlim([0.1 100]); ylim([10 100]);

% Phase Plot
subplot(2,1,2);
semilogx(results_thelen_passive{1}, results_thelen_passive{3}, 'r-', 'LineWidth', 1.5); hold on;
semilogx(results_millard_passive{1}, results_millard_passive{3}, 'b--', 'LineWidth', 1.5);
title('Passive Response: Phase');
xlabel('Frequency (Hz)'); ylabel('Phase (deg)');
grid on; xlim([0.1 100]); ylim([-180 10]);