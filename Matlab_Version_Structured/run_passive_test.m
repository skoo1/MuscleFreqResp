% By Minseung Kim and Seungbum Koo
% KAIST, Daejeon, South Korea
% February 8, 2026

function results = run_passive_test(muscleModel)
    % Sim configuration
    L_mtn_input    = 1.05;
    % L_mn_input     = 1.0;
    % V_m_input      = 0.0;  % static condition
    % V_mt_input     = 0.0;  % static condition
    U_input        = 0.0;
    Amp_input      = 0.005;
    Mass_input     = 0.0;
    Damping_imput  = 0.0;
    SimTime_input  = 120;
    SimDt_input    = 0.001;
    FreqLow_input  = 0.1;
    FreqHigh_input = 100;
    NumFreqSamples = 10;

    % 1. Simulation Configuration
    SimTime    = SimTime_input;
    dt         = SimDt_input;
    time_array = 0 : dt : (SimTime - dt);
    time_len   = length(time_array);

    % Frequency Sweep Parameters
    freq_low  = FreqLow_input;
    freq_high = FreqHigh_input;
    num_freqs = NumFreqSamples;
    freqs     = logspace(log10(freq_low), log10(freq_high), num_freqs);

    % External dynamics parameters (this test is kinematics driven)
    % mass_ext    = Mass_input;
    % damping_ext = Damping_imput;

    % Operating Point
    L_mtn_target = L_mtn_input;
    u0           = U_input;

    % Initialize the muscle state to equilibrium at t=0
    L_mt0 = muscleModel.L_mo * cos(muscleModel.AlphaOpt) + muscleModel.L_ts;
    a0 = u0;
    [~, F_t0, ~] = muscleModel.initialize_static_given_Lmt(a0, L_mtn_target * L_mt0);
    muscleModel0 = muscleModel;

    % Prepare Storage
    results = cell(4, 1);
    sin_f = zeros(num_freqs, 1);
    mag   = zeros(num_freqs, 1);
    pha   = zeros(num_freqs, 1);
    exc   = zeros(num_freqs, 1);

    % Main Frequency Loop
    fprintf('Running Passive Test on %s...\n', muscleModel.MuscleName);

    parfor k = 1 : num_freqs
        tic;
        freq = freqs(k);

        % Data Logging
        F_m_AT = zeros(1, time_len);
        F_m_AT(1) = F_t0;

        muscleModel = muscleModel0; % reset

        amp           = Amp_input * L_mt0 * L_mtn_target;
        L_mt_preset   = L_mt0 * L_mtn_target + amp * sin(2 * pi * freq * time_array);        
        V_mt_preset   = amp * 2 * pi * freq * cos(2 * pi * freq * time_array);
        % A_mt          = -amp * (2 * pi * freq)^2 * sin(2 * pi * freq * time_array);
        % starts from static condition
        % muscleModel.V_m  = 0; % from initialization

        % --- Time Integration Loop ---
        for i = 2 : time_len
        % Update Muscle Dynamics
        % Kinematically driven by given L_mt.
        % Finds the normalized fiber length (L_mn) that satisfies the force balance (F_m == F_t).
        % Quasi-static: While V_m is computed and could be used for the force-velocity multiplier (f_v),
        % F_a remains zero in this scenario since u = 0 and a = 0.
            F_t = muscleModel.updateDynamicsQuasiStatic(dt, u0, L_mt_preset(i), V_mt_preset(i));

            % Record Sensor Force
            F_m_AT(i) = F_t;
        end

        % FFT Analysis
        % Remove Transient (First 10%)
        valid_range = round(time_len * 0.1) : time_len;
        sig_in  = L_mt_preset(valid_range);
        sig_out = F_m_AT(valid_range);

        U_fft = fft(sig_in);
        F_fft = fft(sig_out);

        Fs     = 1 / dt;
        n_x    = length(sig_in);
        f_axis = (0 : n_x - 1) * (Fs / n_x);

        [~, idx] = min(abs(f_axis - freq));

        sin_f(k) = freq;
        mag(k)   = abs(F_fft(idx) / U_fft(idx));
        pha(k)   = angle(F_fft(idx)) - angle(U_fft(idx));
        exc(k)   = abs(U_fft(idx));

        fprintf('%d, %f, %f, %f, %.1f\n', k, freq, mag(k), pha(k), toc);
    end

    % Finalize Results
    % Unwrap Phase
    pha = unwrap(pha) * (180 / pi);

    results{1} = sin_f;
    results{2} = mag;
    results{3} = pha;
    results{4} = exc;

    % Plotting (Graph Drawing)
    % Create a new figure for the Bode Plot
    figure('Name', ['Passive Test - ' muscleModel.MuscleName]);
    sgtitle(sprintf('Bode plot w/ u_{0} = %.2f (%s)', u0, muscleModel.MuscleName), 'FontSize', 16);

    mag_db = 20 * log10(mag); % Convert magnitude to dB

    % Subplot 1: Magnitude Response
    subplot(2, 1, 1);
    semilogx(sin_f, mag_db, '-o', 'MarkerSize', 3, 'LineWidth', 1.5);
    ylabel("Magnitude (dB)", 'FontSize', 12); 
    grid on; 
    xlim([0.1 100]);
    ylim([10 100]);
    title('Magnitude Response');

    % Subplot 2: Phase Response
    subplot(2, 1, 2);
    semilogx(sin_f, pha, '-o', 'MarkerSize', 3, 'LineWidth', 1.5);
    xlabel("Frequency (Hz)", 'FontSize', 12); 
    ylabel("Phase (deg)", 'FontSize', 12); 
    grid on; 
    xlim([0.1 100]);
    ylim([-180 10]);
    title('Phase Response');
end