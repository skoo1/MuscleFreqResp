function analyze_simulation_results(csv_folder, save_folder, model_name, sim_type, u0, L_mn)
    % analyze_simulation_results
    % C++ 시뮬레이션 결과(CSV)를 읽어 Bode Plot을 그리고 .mat 파일로 저장합니다.
    %
    % Inputs:
    %   csv_folder  : CSV 파일들이 있는 폴더 경로 (예: './MMM_Active_results_csv/')
    %   save_folder : 결과를 저장할 폴더 경로 (예: './MMM_results')
    %   model_name  : 모델 이름 ('MMM' or 'TMM')
    %   sim_type    : 시뮬레이션 타입 ('Active' or 'Passive')
    %   u0          : 초기 활성도 (Active=0.5, Passive=0.0)
    %   L_mn        : 정규화된 근섬유 길이 타겟 (보통 1.0)

    fprintf('Analyzing %s - %s Simulation...\n', model_name, sim_type);

    % 1. 주파수 리스트 읽기
    freq_file = fullfile(csv_folder, 'frequency_list.csv');
    if ~isfile(freq_file)
        error('Frequency list file not found: %s', freq_file);
    end
    freq_list = readmatrix(freq_file);
    % C++에서 "index, freq" 형식으로 저장했으므로 2번째 컬럼 사용
    sin_f = freq_list(:, 2)'; 
    numFreq = length(sin_f);

    % 2. 초기 설정
    dt = 0.001; 
    Fs = 1/dt;
    muscleName = 'Soleus';
    
    mag = zeros(1, numFreq);
    pha = zeros(1, numFreq);

    % 3. 주파수별 데이터 처리 루프
    for k = 1:numFreq
        % 파일 읽기 (k-1은 C++ 0-based index)
        csv_name = sprintf('freq_res_%d.csv', k-1);
        file_path = fullfile(csv_folder, csv_name);
        
        if ~isfile(file_path)
            warning('File not found: %s', file_path);
            continue;
        end
        
        % CSV 읽기 (헤더가 있는 경우 readmatrix가 자동으로 숫자만 추출하거나 옵션 필요)
        % 안전하게 NumHeaderLines 옵션 사용 권장하나, 최신 버전은 자동 감지
        try
            data = readmatrix(file_path);
        catch
            % 구버전 호환용
            data = csvread(file_path, 1, 0); 
        end
        
        % 데이터 추출
        % Col 1: Time
        % Col 2: Input (Active='u', Passive='L_mt')
        % Col 3: Output (Force 'F_m_AT')
        u_sig = data(:, 2); 
        f_sig = data(:, 3);
        
        % FFT 분석 범위 설정 (초기 과도 응답 제외, 뒤쪽 90% 사용)
        valid_range = round(length(u_sig) * 0.1) : length(u_sig);
        sig_in = u_sig(valid_range);
        sig_out = f_sig(valid_range);
        
        % FFT 계산
        U_fft = fft(sig_in); 
        F_fft = fft(sig_out);
        
        n_x = length(sig_in);
        f_axis = (0 : n_x - 1) * (Fs / n_x);
        
        % 현재 주파수 인덱스 찾기
        target_freq = sin_f(k);
        [~, idx] = min(abs(f_axis - target_freq));
        
        % Magnitude & Phase 계산
        % Active: Force / Activation (Gain)
        % Passive: Force / Length (Stiffness)
        mag(k) = abs(F_fft(idx)) / abs(U_fft(idx));
        pha(k) = angle(F_fft(idx)) - angle(U_fft(idx));
        
        if mod(k, 20) == 0
            fprintf('  Processing %d/%d: Freq = %.2f Hz\n', k, numFreq, target_freq);
        end
    end

    % 4. 위상 보정 (Unwrap & Degree 변환)
    pha = unwrap(pha) * (180 / pi);

    % 5. 결과 저장 (.mat)
    if ~exist(save_folder, 'dir')
        mkdir(save_folder);
    end
    
    % 파일명 생성 규칙 준수
    % 형식: [L_mn]_[u0]_[muscle]_[age]_[damping]_[mass].mat
    mass_label = 'isometric'; 
    mat_name = sprintf('%.1f_%.1f_%s_YB_wod_%s_%s.mat', ...
                       L_mn, u0, muscleName, model_name, mass_label); 
    % 참고: 파일명에 시뮬레이션 타입(Active/Passive) 구분을 위해 model_name 뒤에 붙임
    % 기존 포맷을 유지하고 싶다면 인자를 조정하세요. 여기서는 식별을 위해 model_name에 포함시킴.
    
    fullPath = fullfile(save_folder, mat_name);
    
    visual_result = {sin_f, mag, pha, []}; 
    savingdata_freq = visual_result;
    
    save(fullPath, 'savingdata_freq');
    fprintf('  Saved results to: %s\n', fullPath);

    % 6. Bode Plot 그리기
    figure('Name', sprintf('%s %s Simulation', model_name, sim_type));
    sgtitle(sprintf('%s %s Bode Plot (u_{0}=%.2f)', model_name, sim_type, u0), 'FontSize', 16);
    cm = hsv(1);
    
    % Magnitude
    subplot(2, 1, 1);
    % Passive(Stiffness)와 Active(Gain)의 물리적 단위가 다르므로 Y축 라벨 구분
    mag_db = 20 * log10(mag);
    semilogx(sin_f, mag_db, 'o-', 'MarkerSize', 3, 'LineWidth', 1, 'Color', cm);
    grid on;
    xlim([0.1 100]);
    ylim([10 100]);
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
    ylim([-180 10]);
    xlabel("Frequency (Hz)", 'FontSize', 14);
    ylabel("Phase (deg)", 'FontSize', 12);
    title('Phase Response');
    
    drawnow;
end