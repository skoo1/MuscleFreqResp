% By Minseung Kim, Seungwoo Yoon and Seungbum Koo
% KAIST, Daejeon, South Korea
% February 23, 2026

function analyze_simulation_results(csv_folder, save_folder, model_name, sim_type, u0, L_mn)

    fprintf('Analyzing %s - %s Simulation...\n', model_name, sim_type);

    freq_file = fullfile(csv_folder, 'frequency_list.csv');
    if ~isfile(freq_file)
        error('Frequency list file not found: %s', freq_file);
    end
    freq_list = readmatrix(freq_file);
    sin_f = freq_list(:, 2)'; 
    numFreq = length(sin_f);

    dt = 0.001; 
    Fs = 1/dt;
    muscleName = 'Soleus';
    
    mag = zeros(1, numFreq);
    pha = zeros(1, numFreq);

    for k = 1:numFreq
        csv_name = sprintf('freq_res_%d.csv', k-1);
        file_path = fullfile(csv_folder, csv_name);
        
        if ~isfile(file_path)
            warning('File not found: %s', file_path);
            continue;
        end
        
        try
            data = readmatrix(file_path);
        catch
            data = csvread(file_path, 1, 0); 
        end
        
        % Col 1: Time
        % Col 2: Input (Active='u', Passive='L_mt')
        % Col 3: Output (Force 'F_m_AT')
        u_sig = data(:, 2); 
        f_sig = data(:, 3);
        
        valid_range = round(length(u_sig) * 0.1) : length(u_sig);
        sig_in = u_sig(valid_range);
        sig_out = f_sig(valid_range);
        
        U_fft = fft(sig_in); 
        F_fft = fft(sig_out);
        
        n_x = length(sig_in);
        f_axis = (0 : n_x - 1) * (Fs / n_x);
        
        target_freq = sin_f(k);
        [~, idx] = min(abs(f_axis - target_freq));
        
        % Active: Force / Activation (Gain)
        % Passive: Force / Length (Stiffness)
        mag(k) = abs(F_fft(idx)) / abs(U_fft(idx));
        pha(k) = angle(F_fft(idx)) - angle(U_fft(idx));
        
        if mod(k, 20) == 0
            fprintf('  Processing %d/%d: Freq = %.2f Hz\n', k, numFreq, target_freq);
        end
    end

    pha = unwrap(pha) * (180 / pi);

    if ~exist(save_folder, 'dir')
        mkdir(save_folder);
    end
    
    % format: [L_mn]_[u0]_[muscle]_[age]_[damping]_[mass].mat
    mass_label = 'isometric'; 
    mat_name = sprintf('%.1f_%.1f_%s_YB_wod_%s_%s.mat', ...
                       L_mn, u0, muscleName, model_name, mass_label); 
    
    fullPath = fullfile(save_folder, mat_name);
    
    visual_result = {sin_f, mag, pha, []}; 
    savingdata_freq = visual_result;
    
    save(fullPath, 'savingdata_freq');
    fprintf('  Saved results to: %s\n', fullPath);

    figure('Name', sprintf('%s %s Simulation', model_name, sim_type));
    sgtitle(sprintf('%s %s Bode Plot (u_{0}=%.2f)', model_name, sim_type, u0), 'FontSize', 16);
    cm = hsv(1);
    
    % Magnitude
    subplot(2, 1, 1);
    mag_db = 20 * log10(mag);
    semilogx(sin_f, mag_db, 'o-', 'MarkerSize', 3, 'LineWidth', 1, 'Color', cm);
    grid on;
    xlim([0.1 100]);
    ylim([10 150]);
    if strcmp(sim_type, 'Active')
        ylabel("Magnitude (dB) [N/\Delta u]", 'FontSize', 12);
    else
        ylabel("Magnitude (dB) [N/m]", 'FontSize', 12); % Stiffness
    end
    title('Magnitude Response');

    % Phase
    subplot(2, 1, 2);
    semilogx(sin_f, pha, 'o-', 'MarkerSize', 3, 'LineWidth', 1, 'Color', cm);
    grid on;
    xlim([0.1 100]);
    ylim([-180 50]);
    xlabel("Frequency (Hz)", 'FontSize', 14);
    ylabel("Phase (deg)", 'FontSize', 12);
    title('Phase Response');
    
    drawnow;
end