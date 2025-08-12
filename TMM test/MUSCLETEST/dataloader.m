%% Descriptions
% By Minseung Kim, 2024-11-01

% sol/gast (muscle type) | YB/OB (aging) | wod/d (damping)
% Format:: TMM_results_KMS_Lmn0_uo_sol_YB_wod.mat %

clc; clear;

n           = input('Enter the number of files to load: ');
datapaths   = cell(n, 1); 
fileNames   = cell(n, 1);

% Example:: datapath1 = 'G:\Research\TMM test\MUSCLETEST\TMM_result(nopenn)_KMS\afterFFT';
% Example:: datapath2 = 'G:\Research\MMM test\MMM\src\MMM_result_KMS\afterFFT';
% Example:: datapath3 = 'G:\Research\NMM test\NMM_result_KMS\afterFFT';
% Example:: fileName1 = '1_0.5_sol_YB_wd1000_aF_n.mat'; 
% 1_0_sol_YB_wod_aF_passive_0kg.mat

for i = 1 : n
    datapaths{i} = input(['Enter path for file ' num2str(i) ': '], 's');
    fileNames{i} = input(['Enter name for file ' num2str(i) ': '], 's');
end

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
    sin_f_data{i}       = savingdata_freq{1};
    mag_data{i}         = savingdata_freq{2};
    phase_data{i}       = savingdata_freq{3};
    
    parts               = split(fileNames{i}, '_');
    L_mn0_values(i)     = str2double(parts{1});
    % muscle_names{i}     = parts{3}; % Muscle name
    muscle_name         = {'Sol', 'Gas', 'TA'};
    age_types{i}        = parts{4};
    damping_types{i}    = parts{5};
    % MMM_variations{i}   = parts{8};
    mass_types          = {'Isometric', '30 kg', '300 kg'};

    last_part           = parts{end};
    variable_name       = extractBefore(last_part, '.mat');
    variable_names{i}   = variable_name; % Save for legend
end

%% Bode plots
% 
figure();

% sg_title_text = sprintf('Results of %s muscle, initial u_o = %.2f, L^M_o = %.2f', ...
%                         muscle_names{1}, u0, L_mn0_values(1));

% sg_title_text = sprintf('MMM-Rigid Results of %s muscle, initial u_o = %.2f, L^M_o = %.2f', ...
%                         muscle_names{1}, u0, L_mn0_values(1));

% sg_title_text = sprintf('TMM Results of %s muscle, initial u_o = %.2f', ...
%     muscle_names{1}, u0);

% sg_title_text = sprintf('MMM-Classic Results of %s muscle, initial u_o = %.2f', ...
%     muscle_names{1}, u0);

% sg_title_text = sprintf('TMM Results of initial u_o = %.2f, L^M_o = %.2f', ...
%       u0, L_mn0_values(1));


% sg = sgtitle(sg_title_text, 'FontSize', 20, 'FontName', 'Times New Roman');

markerStyles = {'o', 's', 'd', '^', 'v', 'p', 'h'}; % n ≤ 7 이라고 가정
lineStyles = {'-', '--', ':', '-.'}; % 선 스타일도 가능하면

cm              = hsv(n);
cutoff_level    = -3;

% nat_freqs = zeros(pnum-1, 1);
for i = 1 : n
    freq_values = sin_f_data{i};
    mag_values = mag_data{i};
    pha_values = phase_data{i};

    mag_dB = 20 * log10(mag_values);

    % % 초기 평균값 구간 (예: 0.1–1.0Hz)
    % plateau_mask = (freq_values >= 0.1) & (freq_values <= 3.0);
    % if nnz(plateau_mask) >= 3
    %     plateau_avg = mean(mag_dB(plateau_mask));
    % else
    %     plateau_avg = mag_dB(1);  % fallback
    % end
    % 
    % % 설정: ±band_dB 만큼 허용
    % band_dB = 0.05;  % 예: ±3dB band
    % clip_min = plateau_avg - band_dB;
    % clip_max = plateau_avg + band_dB;
    % 
    % % 적용할 주파수 범위 (예: 20Hz 이상에서만 적용)
    % clip_mask = (freq_values >= 10);
    % mag_dB(clip_mask) = min(max(mag_dB(clip_mask), clip_min), clip_max);
    % % 
    % % 
    % omega = 2 * pi * freq_values; % Frequency in rad/s

    % Reference magnitude at 0.1 Hz (first index)
    ref_mag_dB = mag_dB(1); % First index corresponds to 0.1 Hz

    cutoff_mag = ref_mag_dB - 3; % Define -3 dB cutoff level
    lower_cutoff_idx = find(mag_dB <= cutoff_mag, 1, 'first');

    if ~isempty(lower_cutoff_idx)
        lower_cutoff_freq = freq_values(lower_cutoff_idx);
    else
        lower_cutoff_freq = NaN;
    end

    % Determine natural frequency as the peak magnitude frequency
    [max_mag, max_idx] = max(mag_dB);
    nat_freq_Hz = freq_values(max_idx);

    % Magnitude plot
    subplot(2, 1, 1);
    semilogx(freq_values, mag_dB, 'o', 'MarkerSize', 1, 'Color', cm(i, :), ...
        'DisplayName', sprintf('%s (Original)', L_mn0_values(i)));
    grid on;
    hold on;

    xlabel("Frequency (Hz)", 'FontSize', 18, 'FontName', 'Times New Roman'); 
    ylabel("|F/X|", 'FontSize', 18, 'FontName', 'Times New Roman'); grid on;
    ylim([85 95]);

    peak_mag = mag_dB(max_idx); % Magnitude at natural frequency
    % plot(freq_values(max_idx), peak_mag, 'x', 'MarkerSize', 10, 'Color', cm(i, :), 'LineWidth', 2);
    % text(freq_values(max_idx), peak_mag + 8, sprintf('$\\mathbf{w_n \\approx %.2f}$ Hz', nat_freq_Hz), ...
         % 'Color', cm(i, :), 'FontSize', 20, 'HorizontalAlignment', 'center', 'Interpreter', 'latex', 'FontName', 'Times New Roman');
    % lower_cutoff_freq = 0;
    % Highlight -3 dB point
    if ~isnan(lower_cutoff_freq)
        % plot(lower_cutoff_freq, cutoff_mag, 'x', 'MarkerSize', 10, 'Color', 'black', 'LineWidth', 2);
        % text(lower_cutoff_freq - 15, cutoff_mag - 10, '-3dB Point', 'Color', cm(i, :), ...
             % 'FontSize', 15, 'HorizontalAlignment', 'center', 'FontName', 'Times New Roman');
    end

    % Define model name manually (e.g., 'TMM', 'MMM', etc.)
    model_name = 'MMM-DEq'; % 수동으로 'TMM' 또는 'MMM'을 입력하세요

    % Bandwidth arrow using line
    if ~isnan(lower_cutoff_freq)
        % Bandwidth line coordinates
        bandwidth_y = cutoff_mag - 20; % Adjust height below cutoff_mag
        % line([0.1, lower_cutoff_freq], [bandwidth_y-1, bandwidth_y-1], 'Color', 'blue', 'LineWidth', 1.5);

        % Left arrowhead (< slightly inside the plot at 0.1 Hz)
        % text(0.1 * 1.25, bandwidth_y, '<', 'Color', 'blue', 'FontSize', 18, 'HorizontalAlignment', 'right', 'FontWeight', 'bold');

        % Right arrowhead (> slightly inside the plot at lower_cutoff_freq)
        % text(lower_cutoff_freq * 0.85, bandwidth_y, '>', 'Color', 'blue', 'FontSize', 18, 'HorizontalAlignment', 'left', 'FontWeight', 'bold');

        bandwidth_value = lower_cutoff_freq - 0.1; % Bandwidth in Hz
        mid_freq = sqrt(0.1 * lower_cutoff_freq);
        if isempty(model_name)
            warning('Model name is not specified. Please assign a value to "model_name".');
        end
        % text(mid_freq, bandwidth_y - 8, sprintf('Bandwidth: %.2f Hz (%s)', bandwidth_value, model_name), ...
             % 'Color', 'blue', 'FontSize', 17, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    end

    % Phase plot
    subplot(2, 1, 2);
    % plot(freq_values, pha_values, 'o', 'MarkerSize', 1, 'Color', cm(i, :), ...
    %     'DisplayName', sprintf('%s (Original)', variable_names{i}));
    semilogx(freq_values, pha_values, 'o', 'MarkerSize', 1, 'Color', cm(i, :), ...
        'DisplayName', sprintf('%s (Original)', damping_types{i}));
    grid on;
    hold on;

    xlabel("Frequency (Hz)", 'FontSize', 18, 'FontName', 'Times New Roman'); 
    ylabel('\angleF/X', 'FontSize', 18, 'FontName', 'Times New Roman');
    grid on;
    ylim([-50 50]);

    % nat_freqs(i) = nat_freq_Hz;
    % 
    % bode_response = mag_values .* (cosd(pha_values) + 1j * sind(pha_values));
    % data = frd(bode_response, omega);
    % 
    % sys = tfest(data, 2);
    % [wn, zeta] = damp(sys);
    % nat_freq_Hz = wn(1) / (2 * pi); % Convert to Hz
    % 
    % % Compute Bode response of fitted system
    % [mag_fit, pha_fit, omega_fit] = bode(sys, omega);
    % mag_fit_dB = 20 * log10(squeeze(mag_fit));
    % pha_fit_deg = squeeze(pha_fit);
    % 
    % subplot(2, 1, 1);
    % semilogx(freq_values, mag_dB, 'o', 'MarkerSize', 1, 'Color', cm(i, :), ...
    %     'DisplayName', sprintf('%s (Original)', variable_names{i}));
    % grid on;
    % hold on;
    % 
    % % semilogx(freq_values, mag_fit_dB, '-', 'LineWidth', 2, 'Color', cm(i, :), ...
    % %     'DisplayName', sprintf('%s (Fitted, f_n = %.2f Hz)', variable_names{i}, nat_freq_Hz));
    % xlabel("Frequency (Hz)", 'FontSize', 10); ylabel("Magnitude (dB)", 'FontSize', 14); grid on;
    % ylim([10 60]);
    % 
    % [~, closest_idx] = min(abs(freq_values - nat_freq_Hz)); % Closest index to natural frequency
    % peak_mag = mag_dB(closest_idx); % Magnitude at natural frequency
    % % plot(freq_values(closest_idx), peak_mag, 'x', 'MarkerSize', 10, 'Color', cm(i, :), 'LineWidth', 2);
    % % text(freq_values(closest_idx), peak_mag + 6, sprintf('$\\mathbf{w_n \\approx %.2f}$ Hz', nat_freq_Hz), ...
    % %  'Color', cm(i, :), 'FontSize', 13, 'HorizontalAlignment', 'center', 'Interpreter', 'latex');
    % 
    % subplot(2, 1, 2);
    % semilogx(freq_values, pha_values, 'o', 'MarkerSize', 1, 'Color', cm(i, :), ...
    %     'DisplayName', sprintf('%s (Original)', variable_names{i}));
    % grid on;
    % hold on;
    % 
    % % semilogx(freq_values, pha_fit_deg, '-', 'LineWidth', 2, 'Color', cm(i, :), ...
    % %     'DisplayName', sprintf('%s (Fitted)', va5riable_names{i}));
    % xlabel("Frequency (Hz)", 'FontSize', 10); ylabel("Phase (deg)", 'FontSize', 14); grid on;
    % ylim([-180 50]);

    % nat_freqs(i) = nat_freq_Hz;

    % legend_entries{i} = sprintf('L^M_0 = %.2f, w_n ≈ %.2f Hz, %s', L_mn0_values(i), nat_freq_Hz, variable_names{i});
    legend_entries{i} = sprintf('L^M_0 = %.2f, w_n ≈ %.2f Hz', L_mn0_values(i), nat_freqs(i));
    % legend_entries{i} = sprintf('%s, w_n ≈ %.2f Hz', muscle_names{i}, nat_freqs(i));
end

% Add legends
subplot(2, 1, 1);
hold on;
plot_handles = [];
for i = 1:n
    plot_handles = [plot_handles; semilogx(NaN, NaN, 'o', 'MarkerSize', 13, 'Color', cm(i, :), ...
        'MarkerFaceColor', cm(i, :), 'DisplayName', sprintf('%s', L_mn0_values(i)))];
end
legend(plot_handles, 'Location', 'EastOutside', 'Box', 'off', 'FontSize', 15);

subplot(2, 1, 2);
hold on;
plot_handles_phase = [];
for i = 1:n
    plot_handles_phase = [plot_handles_phase; semilogx(NaN, NaN, 'o', 'MarkerSize', 13, 'Color', cm(i, :), ...
        'MarkerFaceColor', cm(i, :), 'DisplayName', sprintf('%s', L_mn0_values(i)))];
end
legend(plot_handles_phase, 'Location', 'EastOutside', 'Box', 'off', 'FontSize', 15);

disp('Natural frequencies for each patch: ');
disp(nat_freqs);


%% 단일 Figure 생성
figure();

% 제목 텍스트 설정
% sg_title_text = sprintf('TMM Results of %s muscle, initial excitation = %.2f', ...
%                         muscle_names{1}, u0);
% sg = sgtitle(sg_title_text, 'FontSize', 23, 'FontName', 'Times New Roman');

% 색상 맵 설정
cm = hsv(n);  % n개의 고유한 색상 생성

% 자연 주파수 저장 배열 초기화
nat_freqs = zeros(n, 1);
lower_cutoff_freq = zeros(n, 1);

% Magnitude Plot 생성
hold on;
for i = 1:n
    % 데이터 가져오기
    idx = i;
    freq_values = sin_f_data{idx};
    mag_values  = mag_data{idx};

    % 크기를 dB 단위로 변환
    mag_dB = 20 * log10(mag_values);

    % Reference magnitude at 0.1 Hz (first index)
    ref_mag_dB = mag_dB(1); % First index corresponds to 0.1 Hz
    cutoff_mag = ref_mag_dB - 3; % Define -3 dB cutoff level
    lower_cutoff_idx = find(mag_dB <= cutoff_mag, 1, 'first');

    if ~isempty(lower_cutoff_idx)
        lower_cutoff_freq(idx) = freq_values(lower_cutoff_idx);
    else
        lower_cutoff_freq(idx) = NaN;
    end

    % Step 1: 저주파 plateau 구간에서 평균값 계산 (예: 0.1–1.0 Hz)
    plateau_range = (freq_values >= 0.1) & (freq_values <= 1.3);
    if nnz(plateau_range) >= 3
        mag_plateau = mag_dB(plateau_range);
        std_plateau = std(mag_plateau);

        if std_plateau > 0.5  % 노이즈가 많은 경우에만 평균값 사용
            ref_mag_dB = mean(mag_plateau);
            cutoff_mag = ref_mag_dB - 3;
        else
            ref_mag_dB = mag_dB(find(plateau_range, 1));  % 첫 값 사용
            cutoff_mag = ref_mag_dB - 3;
        end
    else
        ref_mag_dB = mag_dB(1);  % fallback
    end

    % Step 2: 시스템 peak은 plateau 이후에서 탐색 (예: freq >= 2.0Hz)
    peak_search_range = (freq_values >= 1.0);
    if nnz(peak_search_range) >= 1.3
        [~, relative_peak_idx] = max(mag_dB(peak_search_range));
        peak_candidates = find(peak_search_range);
        max_idx = peak_candidates(relative_peak_idx);
    else
        [~, max_idx] = max(mag_dB);  % fallback
    end
    nat_freq_Hz = freq_values(max_idx);
    nat_freqs(idx) = nat_freq_Hz;

    % Step 3: -3dB cutoff는 peak 이후에서 탐색
    if ~isnan(cutoff_mag)
        search_idx = max_idx:length(mag_dB);
        lower_cutoff_idx_local = find(mag_dB(search_idx) <= cutoff_mag, 1, 'first');
        if ~isempty(lower_cutoff_idx_local)
            lower_cutoff_idx = search_idx(lower_cutoff_idx_local);
            lower_cutoff_freq(idx) = freq_values(lower_cutoff_idx);
        else
            lower_cutoff_freq(idx) = NaN;
        end
    else
        lower_cutoff_freq(idx) = NaN;
    end

    % 최대 크기 위치에서 자연 주파수 계산
    [~, max_idx] = max(mag_dB);
    nat_freq_Hz = freq_values(max_idx);
    nat_freqs(idx) = nat_freq_Hz;
end

[~, sorted_idx] = sort(lower_cutoff_freq, 'ascend');

hold on;
for i = 1 : n
    idx = sorted_idx(i);
    freq_values = sin_f_data{idx};
    mag_values  = mag_data{idx};
    mag_dB = 20 * log10(mag_values);
    % Magnitude Plot 추가
    % semilogx(freq_values, mag_dB, 'o-', 'MarkerSize', 4, 'LineWidth', 1.5, ...
    %          'Color', cm(i, :), ...
    %          'DisplayName', sprintf('$L^M_o = %.2f$, Bandwidth: %.2f Hz', L_mn0_values(i), lower_cutoff_freq(i)));
    
    % str = sprintf('$L_{mn0} = %.2f$, Bandwidth: %.2f Hz', L_mn0_values(i), lower_cutoff_freq(i));
    str = sprintf('%s, Bandwidth: %.2f Hz', L_mn0_values(idx), lower_cutoff_freq(idx));

    semilogx(freq_values, mag_dB, 'o-', 'MarkerSize', 4, 'LineWidth', 1.5, ...
             'Color', cm(idx, :), ...
             'DisplayName', str);
    ylim([10 80]);
    
    % === 기본 계산 ===
    bw_value = lower_cutoff_freq(idx);
    end_mag = interp1(freq_values, mag_dB, bw_value);
    plot(bw_value, end_mag, 'X', 'Color', cm(idx,:), ...
        'MarkerSize', 20, 'LineWidth', 2, 'HandleVisibility', 'off');

    % y 위치: 선은 아래로 충분히 떨어뜨리기
    y_start = 53;      % 첫 화살표 y 위치
    delta_y = 10;       % 화살표 간 간격을 더 크게
    y_band  = y_start - delta_y*(i-1);
    x1 = 0.1;
    x2 = bw_value;
    x_mid = 10^((log10(x1)+log10(x2))/2);
    
    % 1) 수직 점선 (data units)
    plot([x2 x2], [end_mag y_band], '--', ...
         'Color', cm(idx,:), 'LineWidth', 2, 'HandleVisibility','off');

    % 가변 삼각형 길이
    arrow_abs_length = min(0.05, 0.5 * (log10(x2) - log10(x1)));
    arrow_abs_height = 3.0;
    
    % 삼각형 시작점
    x2_tail = 10^(log10(x2) - arrow_abs_length);
    
    % 화살표 몸통
    plot([x1 x2_tail], [y_band y_band], '-', ...
         'LineWidth', 8, 'Color', cm(idx,:), 'HandleVisibility','off');
    
    % 화살촉
    x_head = [x2_tail, x2_tail, x2];
    y_head = [y_band - arrow_abs_height/2, ...
              y_band + arrow_abs_height/2, y_band];
    patch(x_head, y_head, cm(idx,:), ...
          'EdgeColor', 'none', 'HandleVisibility','off');

    % 3) 텍스트
    text(x_mid, y_band + 3.5, ...
         sprintf('L_{mn0} = %.2f, Bandwidth: %.2f Hz', L_mn0_values(idx), x2), ...
         'FontSize', 27, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'center', 'FontName', 'Times New Roman', ...
         'Color', 'k', 'HandleVisibility','off');
    % ────────────────────
    
    % % === annotation 화살표 그리기 ===
    % xlim_vals = [0.1 100]; 
    % ylim_vals = [10 80];
    % ax_pos = ax.Position;
    % 
    % % 로그 스케일 → 정규화 (0~1)
    % x1_norm = (log10(x1) - log10(xlim_vals(1))) / (log10(xlim_vals(2)) - log10(xlim_vals(1)));
    % x2_norm = (log10(x2) - log10(xlim_vals(1))) / (log10(xlim_vals(2)) - log10(xlim_vals(1)));
    % y_norm = (y_band - ylim_vals(1)) / (ylim_vals(2) - ylim_vals(1));
    % 
    % % figure 내 정규화 좌표
    % x_fig1 = ax_pos(1) + x1_norm * ax_pos(3);
    % x_fig2 = ax_pos(1) + x2_norm * ax_pos(3);
    % y_fig  = ax_pos(2) + y_norm  * ax_pos(4);

end
hold off;

% 안내 문구 위치 재계산
xl = [0.1 100];
yl = ylim(gca);

% X는 log scale이므로 중간 위치 80%쯤에 배치
x_pos = 10^(log10(xl(1)) + 0.98 * (log10(xl(2)) - log10(xl(1))));
y_pos = yl(1) + 0.95 * (yl(2) - yl(1));  % 거의 최상단

text(x_pos, y_pos, '$\mathbf{\times}$ marker: $-3\ \mathrm{dB}$ cutoff point', ...
    'Interpreter', 'latex', ...
    'HorizontalAlignment', 'right', ...
    'FontSize', 25, 'Color', 'k', ...
    'FontName', 'Times New Roman', ...
    'EdgeColor', 'k', 'BackgroundColor', 'white');

hold off;

% 축 레이블 및 스타일 설정
xlabel("Frequency (Hz)", 'FontSize', 30, 'FontName', 'Times New Roman');
y1 = ylabel({'Magnitude (dB)'}, ...
            'FontSize', 30, 'FontName', 'Times New Roman', 'Rotation', 0);
% ylim([0 60]); % y-축 범위 설정
grid on;

% x축을 로그 스케일로 설정
set(gca, 'XScale', 'log');
set(gca, 'FontSize', 25);

% 범례 추가 (그래프 내부 오른쪽 위에 위치, 줄바꿈)
% legend('Location', 'south', 'Orientation', 'vertical', ...
%        'FontSize', 25, 'Interpreter', 'latex', ...
%        'Box', 'on', 'FontName', 'Times New Roman'); % 범례 박스를 표시

% % 자연 주파수 출력
% disp('Natural frequencies for each patch: ');
% disp(nat_freqs);
