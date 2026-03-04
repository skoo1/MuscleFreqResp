% By Minseung Kim, Seungwoo Yoon and Seungbum Koo
% KAIST, Daejeon, South Korea
% February 23, 2026

% sol/gast (muscle type) | YB/OB (aging) | wod/d (damping)
% Format:: TMM_results_KMS_Lmn0_uo_sol_YB_wod.mat %

% clc; 
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

isPassiveMode   = (abs(u0) < 1e-8);

L_mn0_values    = zeros(1, n);
variable_names  = cell(n, 1);
legend_entries  = cell(n, 1);
muscle_names    = cell(n, 1); 
age_types       = cell(n, 1);
last_part       = cell(n, 1);
damping_types   = cell(n, 1);
StartHZ         = 10;
mwinPts         = 8;
refband         = [9 10];

MMM_variations  = cell(n, 1);
Model_names     = {'TMM', 'MMM-Claasic', 'MMM-DEq', 'MMM-Rigid'};
model_names     = cell(n, 1);

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

    [filepath, name, ext] = fileparts(filePath);

    thisModel = 'Unknown';

    if contains(name, '_DEq')
        thisModel = 'MMM-DEq';
    elseif contains(name, '_Rigid')
        thisModel = 'MMM-Rigid';
    else
        if contains(name, 'TMM', 'IgnoreCase', true)
            thisModel = 'TMM';
        elseif contains(name, 'MMM_Classic', 'IgnoreCase', true)
            thisModel = 'MMM-Classic';
        end
    end

    model_names{i} = thisModel;
end

%% For generating Single Figure
figure();

cm = hsv(n);

nat_freqs           = zeros(n, 1);
lower_cutoff_freq   = zeros(n, 1);

hold on;
for i = 1 : n
    idx         = i;
    freq_values = sin_f_data{idx};
    mag_values  = mag_data{idx};

    mag_dB      = 20 * log10(mag_values);

    peak_search_range = (freq_values >= 1.0);

    if nnz(peak_search_range) >= 2
        [~, relative_peak_idx]  = max(mag_dB(peak_search_range));
        peak_candidates         = find(peak_search_range);
        max_idx                 = peak_candidates(relative_peak_idx);
    else
        [~, max_idx] = max(mag_dB);
    end

    nat_freq_Hz     = freq_values(max_idx);
    nat_freqs(idx)  = nat_freq_Hz;

    plateau_range   = (freq_values >= 0.1) & (freq_values <= 1.3);

    if nnz(plateau_range) >= 5
        ref_mag_dB  = mean(mag_dB(plateau_range));  
    else
        ref_mag_dB  = mag_dB(1);                    
    end

    cutoff_mag      = ref_mag_dB - 3;              

    search_idx = max_idx:length(mag_dB);         
    lower_cutoff_idx_local      = find(mag_dB(search_idx) <= cutoff_mag, 1, 'first');

    if ~isempty(lower_cutoff_idx_local)
        lower_cutoff_idx        = search_idx(lower_cutoff_idx_local);
        lower_cutoff_freq(idx)  = freq_values(lower_cutoff_idx);
    else
        lower_cutoff_freq(idx) = NaN;
    end

end

[~, sorted_idx] = sort(lower_cutoff_freq, 'ascend');

unique_models   = unique(model_names);

if isPassiveMode
    ax_mag      = subplot(2,1,1);
    ax_phase    = subplot(2,1,2);
else
    ax_mag      = gca;
end

hold(ax_mag, 'on');
if isPassiveMode
    hold(ax_phase, 'on');
end

for i = 1 : n
    idx         = sorted_idx(i);
    freq_values = sin_f_data{idx};
    mag_values  = mag_data{idx};
    mag_dB      = 20 * log10(mag_values);
    
    mag_dB_plot = mag_dB;

    if isPassiveMode
        idxHF   = (freq_values >= StartHZ);
        idxRef  = (freq_values >= refband(1)) & (freq_values < refband(2));

        if nnz(idxHF) >= 3
            mag_dB_plot(idxHF) = movmean(mag_dB_plot(idxHF), mwinPts, 'Endpoints', 'shrink');
        end

        if nnz(idxRef) >= 1
            % plateauVal = mean(mag_dB_plot(idxRef));
            % mag_dB_plot(idxHF) = plateauVal;
        end
    end

    str = sprintf('%s, Bandwidth: %.2f Hz', L_mn0_values(idx), lower_cutoff_freq(idx));

    semilogx(ax_mag, freq_values, mag_dB_plot, 'o-', 'MarkerSize', 4, 'LineWidth', 1.5, ...
             'Color', cm(idx, :), ...
             'DisplayName', str);
    
    bw_value = lower_cutoff_freq(idx);

    if ~isnan(bw_value)
        end_mag = interp1(freq_values, mag_dB, bw_value);
        plot(ax_mag, bw_value, end_mag, 'X', 'Color', cm(idx,:), ...
            'MarkerSize', 20, 'LineWidth', 2, 'HandleVisibility', 'off');

        y_start = 52.5;      
        delta_y = 10;      
        y_band  = y_start - delta_y*(i-1);
        x1      = 0.1;
        x2      = bw_value;
        x_mid   = 10^((log10(x1)+log10(x2))/2);
    
        plot(ax_mag, [x2 x2], [end_mag y_band], '--', ...
             'Color', cm(idx,:), 'LineWidth', 2, 'HandleVisibility','off');

        arrow_abs_length = min(0.05, 0.5 * (log10(x2) - log10(x1)));
        arrow_abs_height = 3.0;
        x2_tail = 10^(log10(x2) - arrow_abs_length);
    
        plot(ax_mag, [x1 x2_tail], [y_band y_band], '-', ...
             'LineWidth', 8, 'Color', cm(idx,:), 'HandleVisibility','off');

        x_head = [x2_tail, x2_tail, x2];
        y_head = [y_band - arrow_abs_height/2, ...
                  y_band + arrow_abs_height/2, y_band];
        patch(x_head, y_head, cm(idx,:), ...
              'EdgeColor', 'none', 'Parent', ax_mag, ...
              'HandleVisibility','off');

        thisModel = model_names{idx};
        if isscalar(unique_models)
            label_str = sprintf('L_{mn0} = %.2f, Bandwidth: %.2f Hz', ...
                                L_mn0_values(idx), x2);
        else
            label_str = sprintf('%s, Bandwidth: %.2f Hz', thisModel, x2);
        end

        text(ax_mag, x_mid, y_band + 2.9, label_str, ...
             'FontSize', 27, 'FontWeight', 'bold', ...
             'HorizontalAlignment', 'center', 'FontName', 'Times New Roman', ...
             'Color', 'k', 'HandleVisibility','off');
    end

    if isPassiveMode
        pha_values = phase_data{idx};
        semilogx(ax_phase, freq_values, pha_values, 'o-', ...
                 'MarkerSize', 4, 'LineWidth', 1.5, ...
                 'Color', cm(idx, :), ...
                 'DisplayName', sprintf('L_{mn0} = %.2f', L_mn0_values(idx)));
    end
end
hold off;

xlabel(ax_mag, "Frequency (Hz)", 'FontSize', 30, 'FontName', 'Times New Roman');
y1 = ylabel(ax_mag, {'Magnitude (dB)'}, ...
            'FontSize', 30, 'FontName', 'Times New Roman', 'Rotation', 0);
set(ax_mag, 'XScale', 'log');
set(ax_mag, 'FontSize', 25);
grid(ax_mag, 'on');

if isPassiveMode
    ylim(ax_mag, [85 115]);

    xlabel(ax_phase, "Frequency (Hz)", 'FontSize', 30, 'FontName', 'Times New Roman', 'Rotation', 0);
    ylabel(ax_phase, 'Phase (deg)', 'FontSize', 30, 'FontName', 'Times New Roman', 'Rotation', 0);
    set(ax_phase, 'XScale', 'log');
    set(ax_phase, 'FontSize', 25);
    grid(ax_phase, 'on');
    ylim(ax_phase, [0 110]);
else
    ylim(ax_mag, [10 80]);
end

if ~isPassiveMode
    xl = [0.1 100];
    yl = ylim(ax_mag);
    
    x_pos = 10^(log10(xl(1)) + 0.98 * (log10(xl(2)) - log10(xl(1))));
    y_pos = yl(1) + 0.95 * (yl(2) - yl(1));
    
    text(ax_mag, x_pos, y_pos, ...
        '$\mathbf{\times}$ marker: $-3\ \mathrm{dB}$ cutoff point', ...
        'Interpreter', 'latex', ...
        'HorizontalAlignment', 'right', ...
        'FontSize', 25, 'Color', 'k', ...
        'FontName', 'Times New Roman', ...
        'EdgeColor', 'k', 'BackgroundColor', 'white');
end

hold off;

uniqueModels    = unique(model_names);
nModels         = numel(uniqueModels);
model_handles   = gobjects(nModels, 1);

hold(ax_mag, 'on');
for m       = 1 : nModels
    model   = uniqueModels{m};
    idx_example = find(strcmp(model_names, model), 1, 'first');  

    if isempty(idx_example)
        continue;
    end

    c = cm(idx_example, :);

    model_handles(m) = semilogx(ax_mag, NaN, NaN, 'o', ...
        'MarkerSize', 20, ...
        'Color', c, ...
        'MarkerFaceColor', c, ...
        'DisplayName', model);
end

if isPassiveMode
    legend(ax_mag, model_handles, uniqueModels, ...
        'Location', 'Eastoutside', ...
        'Orientation', 'horizontal', ...
        'Box', 'off', 'Orientation', 'vertical', ...
        'FontSize', 22);
    legend(ax_phase, model_handles, uniqueModels, ...
        'Location', 'Eastoutside', ...
        'Orientation', 'horizontal', ...
        'Box', 'off', 'Orientation', 'vertical', ...
        'FontSize', 22);
end
