%% Descriptions 
% By Minseung Kim, 2024-11-01
% Another revision at 2025-11-27

% sol/gast (muscle type) | YB/OB (aging) | wod/d (damping)
% Format:: TMM_results_KMS_Lmn0_uo_sol_YB_wod.mat %

clc; 

% clear;

if evalin('base','exist(''GUI_Bode_files'',''var'')')
    GUI_Bode_files = evalin('base','GUI_Bode_files');

    if ischar(GUI_Bode_files)
        GUI_Bode_files = {GUI_Bode_files};
    end

    n         = numel(GUI_Bode_files);
    datapaths = cell(n,1);
    fileNames = cell(n,1);

    for i = 1:n
        [dp, fn, ext] = fileparts(GUI_Bode_files{i});
        datapaths{i} = dp;
        fileNames{i} = [fn ext];
    end

else
    n           = input('Enter the number of files to load: ');
    datapaths   = cell(n, 1); 
    fileNames   = cell(n, 1);
    
    % Example:: datapath1 = 'C:\Users\Minseung Kim\Desktop\Github\TMM test\TMM\TMM_result\afterFFT';
    % Example:: fileName1 = '1_0.5_sol_YB_wod_aF_n.mat'; 
    
    for i = 1 : n
        datapaths{i} = input(['Enter path for file ' num2str(i) ': '], 's');
        fileNames{i} = input(['Enter name for file ' num2str(i) ': '], 's');
    end
end
%%

sin_f_data      = cell(n, 1);
mag_data        = cell(n, 1);
phase_data      = cell(n, 1);

parts           = split(fileNames{1}, '_');
u0_str          = parts{2};
u0              = str2double(u0_str);

L_mn0_values    = zeros(1, n);
nat_freqs       = zeros(1, n);
variable_names  = cell(n, 1);
legend_entries  = cell(n, 1);
muscle_names    = cell(n, 1); 
age_types       = cell(n, 1);
last_part       = cell(n, 1);
damping_types   = cell(n, 1);
MMM_variations  = cell(n, 1);
Model_names     = {'TMM', 'MMM', 'MMM-DEq', 'MMM-Rigid'};

for i = 1 : n
    filePath            = fullfile(datapaths{i}, fileNames{i});
    data                = load(filePath);

    savingdata_freq     = data.savingdata_freq;
    
    if istable(savingdata_freq)
        sin_f_data{i}   = savingdata_freq{:, 1};   % freq
        mag_data{i}     = savingdata_freq{:, 2};   % mag
        phase_data{i}   = savingdata_freq{:, 3};   % phase
        
    elseif iscell(savingdata_freq)
        sin_f_data{i}   = savingdata_freq{1};
        mag_data{i}     = savingdata_freq{2};
        phase_data{i}   = savingdata_freq{3};
        
    elseif isnumeric(savingdata_freq)
        if size(savingdata_freq,2) < 3
            error('savingdata_freq (numeric) must have at least 3 columns, but size is %dx%d.', ...
                  size(savingdata_freq,1), size(savingdata_freq,2));
        end

        sin_f_data{i}  = savingdata_freq(:, 1);
        mag_data{i}    = savingdata_freq(:, 2);
        phase_data{i}  = savingdata_freq(:, 3);

    else
        error('Unsupported savingdata_freq type: %s', class(savingdata_freq));
    end
    
    parts               = split(fileNames{i}, '_');
    L_mn0_values(i)     = str2double(parts{1});
    muscle_name         = {'Sol', 'Gas', 'TA'};
    age_types{i}        = parts{4};
    damping_types{i}    = parts{5};
    mass_types          = {'Isometric', '30 kg', '300 kg'};

    last_part           = parts{end};
    variable_name       = extractBefore(last_part, '.mat');
    variable_names{i}   = variable_name;
end

%% Bode plots
% figure();
% 
% markerStyles = {'o', 's', 'd', '^', 'v', 'p', 'h'};
% lineStyles = {'-', '--', ':', '-.'};
% 
% cm              = hsv(n);
% cutoff_level    = -3;
% 
% for i = 1 : n
%     freq_values = sin_f_data{i};
%     mag_values = mag_data{i};
%     pha_values = phase_data{i};
% 
%     mag_dB = 20 * log10(mag_values);
% 
%     % Reference magnitude at 0.1 Hz (first index)
%     ref_mag_dB = mag_dB(1); % First index corresponds to 0.1 Hz
% 
%     cutoff_mag = ref_mag_dB - 3; % Define -3 dB cutoff level
%     lower_cutoff_idx = find(mag_dB <= cutoff_mag, 1, 'first');
% 
%     if ~isempty(lower_cutoff_idx)
%         lower_cutoff_freq = freq_values(lower_cutoff_idx);
%     else
%         lower_cutoff_freq = NaN;
%     end
% 
%     % Determine natural frequency as the peak magnitude frequency
%     [max_mag, max_idx] = max(mag_dB);
%     nat_freq_Hz = freq_values(max_idx);
% 
%     % Magnitude plot
%     subplot(2, 1, 1);
%     semilogx(freq_values, mag_dB, 'o', 'MarkerSize', 1, 'Color', cm(i, :), ...
%         'DisplayName', sprintf('%s (Original)', L_mn0_values(i)));
%     grid on;
%     hold on;
% 
%     xlabel("Frequency (Hz)", 'FontSize', 18, 'FontName', 'Times New Roman'); 
%     ylabel("|F/X|", 'FontSize', 18, 'FontName', 'Times New Roman'); grid on;
%     ylim([85 95]);
% 
%     peak_mag = mag_dB(max_idx); % Magnitude at natural frequency
%     % plot(freq_values(max_idx), peak_mag, 'x', 'MarkerSize', 10, 'Color', cm(i, :), 'LineWidth', 2);
%     % text(freq_values(max_idx), peak_mag + 8, sprintf('$\\mathbf{w_n \\approx %.2f}$ Hz', nat_freq_Hz), ...
%          % 'Color', cm(i, :), 'FontSize', 20, 'HorizontalAlignment', 'center', 'Interpreter', 'latex', 'FontName', 'Times New Roman');
%     % lower_cutoff_freq = 0;
%     % Highlight -3 dB point
%     if ~isnan(lower_cutoff_freq)
%         % plot(lower_cutoff_freq, cutoff_mag, 'x', 'MarkerSize', 10, 'Color', 'black', 'LineWidth', 2);
%         % text(lower_cutoff_freq - 15, cutoff_mag - 10, '-3dB Point', 'Color', cm(i, :), ...
%              % 'FontSize', 15, 'HorizontalAlignment', 'center', 'FontName', 'Times New Roman');
%     end
% 
%     % Define model name manually (e.g., 'TMM', 'MMM', etc.)
%     model_name = 'MMM-DEq'; % Type manually
% 
%     % Bandwidth arrow using line
%     if ~isnan(lower_cutoff_freq)
% 
%         bandwidth_y = cutoff_mag - 20; % Adjust height below cutoff_mag
%         % line([0.1, lower_cutoff_freq], [bandwidth_y-1, bandwidth_y-1], 'Color', 'blue', 'LineWidth', 1.5);
% 
%         % Left arrowhead (< slightly inside the plot at 0.1 Hz)
%         % text(0.1 * 1.25, bandwidth_y, '<', 'Color', 'blue', 'FontSize', 18, 'HorizontalAlignment', 'right', 'FontWeight', 'bold');
% 
%         % Right arrowhead (> slightly inside the plot at lower_cutoff_freq)
%         % text(lower_cutoff_freq * 0.85, bandwidth_y, '>', 'Color', 'blue', 'FontSize', 18, 'HorizontalAlignment', 'left', 'FontWeight', 'bold');
% 
%         bandwidth_value = lower_cutoff_freq - 0.1; % Bandwidth in Hz
%         mid_freq = sqrt(0.1 * lower_cutoff_freq);
%         if isempty(model_name)
%             warning('Model name is not specified. Please assign a value to "model_name".');
%         end
%         % text(mid_freq, bandwidth_y - 8, sprintf('Bandwidth: %.2f Hz (%s)', bandwidth_value, model_name), ...
%              % 'Color', 'blue', 'FontSize', 17, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontName', 'Times New Roman');
%     end
% 
%     % Phase plot
%     subplot(2, 1, 2);
%     % plot(freq_values, pha_values, 'o', 'MarkerSize', 1, 'Color', cm(i, :), ...
%     %     'DisplayName', sprintf('%s (Original)', variable_names{i}));
%     semilogx(freq_values, pha_values, 'o', 'MarkerSize', 1, 'Color', cm(i, :), ...
%         'DisplayName', sprintf('%s (Original)', damping_types{i}));
%     grid on;
%     hold on;
% 
%     xlabel("Frequency (Hz)", 'FontSize', 18, 'FontName', 'Times New Roman'); 
%     ylabel('\angleF/X', 'FontSize', 18, 'FontName', 'Times New Roman');
%     grid on;
%     ylim([-50 50]);
% 
%     legend_entries{i} = sprintf('L^M_0 = %.2f, w_n ≈ %.2f Hz', L_mn0_values(i), nat_freqs(i));
% end
% 
% % Add legends
% subplot(2, 1, 1);
% hold on;
% plot_handles = [];
% for i = 1:n
%     plot_handles = [plot_handles; semilogx(NaN, NaN, 'o', 'MarkerSize', 13, 'Color', cm(i, :), ...
%         'MarkerFaceColor', cm(i, :), 'DisplayName', sprintf('%s', L_mn0_values(i)))];
% end
% legend(plot_handles, 'Location', 'EastOutside', 'Box', 'off', 'FontSize', 15);
% 
% subplot(2, 1, 2);
% hold on;
% plot_handles_phase = [];
% for i = 1:n
%     plot_handles_phase = [plot_handles_phase; semilogx(NaN, NaN, 'o', 'MarkerSize', 13, 'Color', cm(i, :), ...
%         'MarkerFaceColor', cm(i, :), 'DisplayName', sprintf('%s', L_mn0_values(i)))];
% end
% legend(plot_handles_phase, 'Location', 'EastOutside', 'Box', 'off', 'FontSize', 15);
% 
% disp('Natural frequencies for each patch: ');
% disp(nat_freqs);


%% For generating Single Figure
figure();

% sg_title_text = sprintf('TMM Results of %s muscle, initial excitation = %.2f', ...
%                         muscle_names{1}, u0);
% sg = sgtitle(sg_title_text, 'FontSize', 23, 'FontName', 'Times New Roman');

cm = hsv(n);

nat_freqs = zeros(n, 1);
lower_cutoff_freq = zeros(n, 1);

hold on;
for i = 1:n
    idx = i;
    freq_values = sin_f_data{idx};
    mag_values  = mag_data{idx};

    mag_dB = 20 * log10(mag_values);

    peak_search_range = (freq_values >= 1.0);
    if nnz(peak_search_range) >= 2
        [~, relative_peak_idx] = max(mag_dB(peak_search_range));
        peak_candidates = find(peak_search_range);
        max_idx = peak_candidates(relative_peak_idx);
    else
        [~, max_idx] = max(mag_dB);
    end
    nat_freq_Hz   = freq_values(max_idx);
    nat_freqs(idx) = nat_freq_Hz;

    plateau_range = (freq_values >= 0.1) & (freq_values <= 1.3);
    if nnz(plateau_range) >= 5
        ref_mag_dB = mean(mag_dB(plateau_range));  
    else
        ref_mag_dB = mag_dB(1);                    
    end
    cutoff_mag = ref_mag_dB - 3;              

    search_idx = max_idx:length(mag_dB);         
    lower_cutoff_idx_local = find(mag_dB(search_idx) <= cutoff_mag, 1, 'first');

    if ~isempty(lower_cutoff_idx_local)
        lower_cutoff_idx   = search_idx(lower_cutoff_idx_local);
        lower_cutoff_freq(idx) = freq_values(lower_cutoff_idx);
    else
        lower_cutoff_freq(idx) = NaN;
    end

end

[~, sorted_idx] = sort(lower_cutoff_freq, 'ascend');

hold on;
for i = 1 : n
    idx = sorted_idx(i);
    freq_values = sin_f_data{idx};
    mag_values  = mag_data{idx};
    mag_dB = 20 * log10(mag_values);

    str = sprintf('%s, Bandwidth: %.2f Hz', L_mn0_values(idx), lower_cutoff_freq(idx));

    semilogx(freq_values, mag_dB, 'o-', 'MarkerSize', 4, 'LineWidth', 1.5, ...
             'Color', cm(idx, :), ...
             'DisplayName', str);
    ylim([10 80]);
    
    bw_value = lower_cutoff_freq(idx);
    end_mag = interp1(freq_values, mag_dB, bw_value);
    plot(bw_value, end_mag, 'X', 'Color', cm(idx,:), ...
        'MarkerSize', 20, 'LineWidth', 2, 'HandleVisibility', 'off');

    y_start = 53;      
    delta_y = 10;      
    y_band  = y_start - delta_y*(i-1);
    x1 = 0.1;
    x2 = bw_value;
    x_mid = 10^((log10(x1)+log10(x2))/2);
    
    plot([x2 x2], [end_mag y_band], '--', ...
         'Color', cm(idx,:), 'LineWidth', 2, 'HandleVisibility','off');

    arrow_abs_length = min(0.05, 0.5 * (log10(x2) - log10(x1)));
    arrow_abs_height = 3.0;
    
    x2_tail = 10^(log10(x2) - arrow_abs_length);
    
    plot([x1 x2_tail], [y_band y_band], '-', ...
         'LineWidth', 8, 'Color', cm(idx,:), 'HandleVisibility','off');

    x_head = [x2_tail, x2_tail, x2];
    y_head = [y_band - arrow_abs_height/2, ...
              y_band + arrow_abs_height/2, y_band];
    patch(x_head, y_head, cm(idx,:), ...
          'EdgeColor', 'none', 'HandleVisibility','off');

    text(x_mid, y_band + 3.5, ...
         sprintf('L_{mn0} = %.2f, Bandwidth: %.2f Hz', L_mn0_values(idx), x2), ...
         'FontSize', 27, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'center', 'FontName', 'Times New Roman', ...
         'Color', 'k', 'HandleVisibility','off');

end
hold off;

xl = [0.1 100];
yl = ylim(gca);

x_pos = 10^(log10(xl(1)) + 0.98 * (log10(xl(2)) - log10(xl(1))));
y_pos = yl(1) + 0.95 * (yl(2) - yl(1));

text(x_pos, y_pos, '$\mathbf{\times}$ marker: $-3\ \mathrm{dB}$ cutoff point', ...
    'Interpreter', 'latex', ...
    'HorizontalAlignment', 'right', ...
    'FontSize', 25, 'Color', 'k', ...
    'FontName', 'Times New Roman', ...
    'EdgeColor', 'k', 'BackgroundColor', 'white');

hold off;

xlabel("Frequency (Hz)", 'FontSize', 30, 'FontName', 'Times New Roman');
y1 = ylabel({'Magnitude (dB)'}, ...
            'FontSize', 30, 'FontName', 'Times New Roman', 'Rotation', 0);
grid on;

set(gca, 'XScale', 'log');
set(gca, 'FontSize', 25);
