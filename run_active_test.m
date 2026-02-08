% By Minseung Kim and Seungbum Koo
% KAIST, Daejeon, South Korea
% February 8, 2026

function results = run_active_test(muscleModel)
    % Sim configuration
    L_mn_input     = 1.0;
    V_m_input      = 0.0;  % static condition
    V_mt_input     = 0.0;  % static condition
    U_input        = 0.5;
    Amp_input      = 0.01;
    Mass_input     = 3e9;
    Damping_imput  = 0.0;
    SimTime_input  = 120;
    SimDt_input    = 0.001;
    FreqLow_input  = 0.1;
    FreqHigh_input = 100;
    NumFreqSamples = 10;

    % Simulation Configuration
    SimTime    = SimTime_input;
    dt         = SimDt_input;
    time_array = 0 : dt : (SimTime - dt);
    time_len   = length(time_array);

    % Frequency Sweep Parameters
    freq_low  = FreqLow_input;
    freq_high = FreqHigh_input;
    num_freqs = NumFreqSamples;
    freqs     = logspace(log10(freq_low), log10(freq_high), num_freqs);

    % External dynamics parameters
    mass_ext    = Mass_input;
    damping_ext = Damping_imput;

    % Operating Point
    L_mn_target = L_mn_input;
    V_m_target  = V_m_input;
    V_mt_target = V_mt_input;
    u0 = U_input;

    % Muscle parameter initialization
    L_mn0 = L_mn_target;
    L_m0  = L_mn0 * muscleModel.L_mo;
    V_mt0 = V_mt_target;
    V_m0  = V_m_target;

    % Initialize the muscle state to equilibrium at t=0
    a0 = u0;
    [L_mt0, F_ext_equil, Alpha0] = muscleModel.initialize_static_given_Lm(a0, L_m0);
    muscleModel0 = muscleModel;

    % Prepare Storage
    results = cell(4, 1);
    sin_f = zeros(num_freqs, 1);
    mag   = zeros(num_freqs, 1);
    pha   = zeros(num_freqs, 1);
    exc   = zeros(num_freqs, 1);

    % Main Frequency Loop
    fprintf('Running Active Test on %s...\n', muscleModel.MuscleName);

    parfor k = 1 : num_freqs
        tic;
        freq = freqs(k);

        % Data Logging
        F_m_AT = zeros(1, time_len);

        % Generate Input Signal (Sine Wave)
        amp = Amp_input;
        u_mean = u0;
        u_init = u0;
        u = sinwave(freq, time_array, u_init, u_mean, amp);

        % Clamp Activation (0 ~ 1)
        u(u < 0) = 0; 
        u(u > 1) = 1;

        muscleModel = muscleModel0; % reset
        muscleModel.a    = u(1);
        muscleModel.L_m  = L_mn0 * muscleModel.L_mo;
        muscleModel.V_m  = V_m0;
        muscleModel.L_mt = L_mt0;
        muscleModel.V_mt = V_mt0;
        muscleModel.L_mn = L_mn0;
        muscleModel.F_t  = F_ext_equil;
        muscleModel.F_m  = F_ext_equil/cos(Alpha0);

        % --- Time Integration Loop ---
        for i = 2 : time_len
            % Update Muscle Dynamics
            % Interaction with external environment: F_ext_equil, mass_ext, damping_ext
            % Compute the fiber force along the tendon (F_m_AT) or tendon force (F_t) based on u.
            % Internal states, including fiber length (L_m) and tendon length (L_t), are updated internally.
            F_t = muscleModel.updateDynamics(dt, u(i), F_ext_equil, mass_ext, damping_ext);

            % Record Sensor Force
            F_m_AT(i) = F_t;
        end

        % FFT Analysis
        % Remove Transient (First 10%)
        valid_range = round(time_len * 0.1) : time_len;
        
        sig_in  = u(valid_range);
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
    figure('Name', ['Active Test - ' muscleModel.MuscleName]);
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

function input = sinwave(frq, t, a0, a, b)
    if b == 0
        phase_shift = 0;
    else
        val = (a0 - a) / b;
        % asin range protection
        if val > 1, val = 1; end
        if val < -1, val = -1; end
        phase_shift = asin(val);
    end

    input = a + sin(2 * pi * frq * t + phase_shift) * b;
end