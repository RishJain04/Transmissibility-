
%% FREQUENCY-DEPENDENT DAMPING IDENTIFICATION
% LOWESS smoothing + Analytical damping identification
%
% Transmissibility:
% T = Xout / Xin
% T_dB = 20*log10(Xout/Xin)
%
% Outputs:
% 1. Transmissibility curve
% 2. Frequency-dependent damping curve
% 3. LOWESS-smoothed damping curve
% 4. Apparant fitting

close all;
clear;
clc;

%% =========================================================
% 1. DEFINE EXPERIMENTAL DATASETS
% ==========================================================

datasets = {
'Low Density aerogel', 'in.xlsx', 'out.xlsx'    
'High Density aerogel', 'foam input.xlsx', 'foam output.xlsx';
    'Medium Density aerogel', 'i_f_m.xlsx', 'o_f_m.xlsx';
    
};

% LOWESS smoothing parameters
transmissibility_span = 0.01;
damping_span = 0.01;

% Maximum frequency for plotting and analysis
max_frequency = 500;

% Minimum input FFT amplitude threshold
input_threshold_factor = 0.01;

% Maximum physically accepted damping ratio
max_damping = 2;

% Summary storage
Summary = table();

%% =========================================================
% 2. PROCESS EACH EXPERIMENTAL DATASET
% ==========================================================

for dataset_index = 1:size(datasets,1)

    group_name = datasets{dataset_index,1};

    input_file = datasets{dataset_index,2};
    output_file = datasets{dataset_index,3};

    fprintf('\n====================================\n');
    fprintf('Processing: %s\n',group_name);
    fprintf('====================================\n');

    %% 2.1 READ EXPERIMENTAL DATA

    Data_input = readmatrix(input_file);
    Data_output = readmatrix(output_file);

    time_in = Data_input(:,1);
    disp_in = Data_input(:,2);

    time_out = Data_output(:,1);
    disp_out = Data_output(:,2);

    % Check data lengths

    if length(time_in) ~= length(time_out)
        warning('Input/output lengths do not match for %s.', ...
            group_name);
        continue;
    end

    % Check time vectors

    if max(abs(time_in-time_out)) > 1e-8
        warning('Input/output time vectors do not match for %s.', ...
            group_name);
        continue;
    end

    % Check sampling interval

    dt_vector = diff(time_in);
    dt = mean(dt_vector);

    %if dt <= 0 || max(abs(dt_vector-dt)) > 1e-6*dt
        %warning('Invalid or nonuniform sampling for %s.', ...
           % group_name);
        %continue;
    %end

    Fs = 1/dt;

    %% =====================================================
    % 3. PREPROCESS EXPERIMENTAL SIGNALS
    % ======================================================

    % Remove mean and linear trend

    disp_in = detrend(disp_in);
    disp_out = detrend(disp_out);

    N = length(disp_in);

    %% =====================================================
    % 4. FFT ANALYSIS
    % ======================================================

    Y_in = fft(disp_in);
    Y_out = fft(disp_out);

    n = floor(N/2)+1;

    freq = (0:n-1)'*Fs/N;

    % Single-sided FFT amplitudes

    FFT_in = abs(Y_in(1:n))/N;
    FFT_out = abs(Y_out(1:n))/N;

    % Correct single-sided amplitudes

    if rem(N,2) == 0

        FFT_in(2:end-1) = 2*FFT_in(2:end-1);
        FFT_out(2:end-1) = 2*FFT_out(2:end-1);

    else

        FFT_in(2:end) = 2*FFT_in(2:end);
        FFT_out(2:end) = 2*FFT_out(2:end);

    end

    %% =====================================================
    % 5. CALCULATE LINEAR TRANSMISSIBILITY
    % ======================================================

    TR = NaN(size(freq));

    threshold = input_threshold_factor*max(FFT_in);

    valid_input = FFT_in > threshold;

    TR(valid_input) = ...
        FFT_out(valid_input)./FFT_in(valid_input);

    %% =====================================================
    % 6. CONVERT TRANSMISSIBILITY TO dB
    % ======================================================

    TR_dB = NaN(size(freq));

    valid_TR = isfinite(TR) & TR > 0;

    TR_dB(valid_TR) = 20*log10(TR(valid_TR));

    %% =====================================================
    % 7. LOWESS SMOOTHING OF TRANSMISSIBILITY
    % ======================================================

    valid_smooth = isfinite(TR) & freq > 0;

    freq_valid = freq(valid_smooth);
    TR_valid = TR(valid_smooth);

    if length(freq_valid) < 5

        warning('Insufficient valid data for %s.',group_name);
        continue;

    end

    % LOWESS smoothing in linear transmissibility

    TR_smooth = smooth(freq_valid,TR_valid, ...
        transmissibility_span,'lowess');

    freq_smooth = freq_valid;

    % Convert smoothed transmissibility to dB

    TR_smooth_dB = NaN(size(TR_smooth));

    positive_smooth = TR_smooth > 0;

    TR_smooth_dB(positive_smooth) = ...
        20*log10(TR_smooth(positive_smooth));

%% =====================================================
% 8. IDENTIFY NATURAL FREQUENCY
% ======================================================

valid_peak = isfinite(TR_smooth) & TR_smooth > 0;

if ~any(valid_peak)

    warning('No valid transmissibility peak for %s.', ...
        group_name);
    continue;

end

valid_freq = freq_smooth(valid_peak);
valid_trans = TR_smooth(valid_peak);

% Find all local peaks
[peak_values, peak_locations] = findpeaks(valid_trans);

peak_frequencies = valid_freq(peak_locations);

% Sort peaks from highest to lowest amplitude
[~, sort_index] = sort(peak_values,'descend');

peak_values = peak_values(sort_index);
peak_frequencies = peak_frequencies(sort_index);

% Select the highest peak initially
peak_index = 1;

% If the highest peak is below 10 Hz,
% select the next highest peak above 10 Hz

if peak_frequencies(1) < 10

    next_peak = find(peak_frequencies >= 10,1,'first');

    if isempty(next_peak)

        warning('No peak above 10 Hz found for %s.',group_name);
        continue;

    end

    peak_index = next_peak;

end

% Identified natural frequency
fn = peak_frequencies(peak_index);

peak_value = peak_values(peak_index);

fprintf('Selected Natural Frequency = %.4f Hz\n',fn);
fprintf('Selected Peak Transmissibility = %.4f\n',peak_value);
    %% =====================================================
    % 9. FREQUENCY-DEPENDENT DAMPING IDENTIFICATION
    % ======================================================

    zeta_freq = NaN(size(freq_smooth));

    for k = 1:length(freq_smooth)

        if ~isfinite(TR_smooth(k)) || ...
                TR_smooth(k) <= 0 || freq_smooth(k) <= 0

            continue;

        end

        % Frequency ratio

        r = freq_smooth(k)/fn;

        Texp = TR_smooth(k);

        % Analytical inversion

        A = Texp^2*(1-r^2)^2 - 1;

        B = 4*r^2*(Texp^2-1);

        if abs(B) < 1e-12
            continue;
        end

        zeta_squared = -A/B;

        if isfinite(zeta_squared) && zeta_squared >= 0

            zeta_freq(k) = sqrt(zeta_squared);

        end

    end

    %% =====================================================
    % 10. FILTER UNRELIABLE DAMPING VALUES
    % ======================================================

    zeta_freq(zeta_freq > max_damping) = NaN;

    zeta_freq(zeta_freq < 0) = NaN;

    %% =====================================================
    % 11. LOWESS SMOOTHING OF DAMPING CURVE
    % ======================================================

    valid_zeta = isfinite(zeta_freq);

    zeta_smooth = NaN(size(zeta_freq));

    if nnz(valid_zeta) > 5

        zeta_smooth(valid_zeta) = smooth( ...
            freq_smooth(valid_zeta), ...
            zeta_freq(valid_zeta), ...
            damping_span, ...
            'lowess');

    end

    %% =====================================================
    % 12. DAMPING AT NATURAL FREQUENCY
    % ======================================================

    if any(valid_zeta)

        damping_at_fn = interp1( ...
            freq_smooth(valid_zeta), ...
            zeta_freq(valid_zeta), ...
            fn,'nearest','extrap');

    else

        damping_at_fn = NaN;

    end

    fprintf('Damping at natural frequency = %.5f\n', ...
        damping_at_fn);


%% =====================================================
% 13. PLOT ACTUAL, LOWESS AND APPARENT FITTING CURVE
% ======================================================

figure('Color','w','Position',[100 100 1000 600]);

% 1. Actual experimental transmissibility

valid_actual = isfinite(TR_dB);

h1 = plot(freq(valid_actual),TR_dB(valid_actual), ...
    'Color',[0.65 0.75 0.90], ...
    'LineWidth',0.9);

hold on;

% 2. LOWESS-smoothed transmissibility

h2 = plot(freq_smooth,TR_smooth_dB, ...
    'r','LineWidth',2);

% 3. Apparent theoretical fitting curve

h3 = gobjects(1);

if isfinite(fn) && isfinite(damping_at_fn) && damping_at_fn > 0

    f_fit = linspace(0.1,max_frequency,2000);

    T_fit_dB = transmissibility_model( ...
        f_fit,fn,damping_at_fn);

    h3 = plot(f_fit,T_fit_dB, ...
        'k--','LineWidth',2);

end

% 4. Identified natural frequency

h4 = xline(fn,'b--', ...
    sprintf('Natural Frequency = %.2f Hz',fn), ...
    'LineWidth',1.5);

% Labels

xlabel('Frequency (Hz)','FontSize',12);

ylabel('Transmissibility (dB)','FontSize',12);

title(['Experimental Transmissibility - ',group_name], ...
    'FontSize',14);

% Legend using explicit handles

if isgraphics(h3)

    legend([h1,h2,h3,h4], ...
        {'Actual Experimental Transmissibility', ...
         'LOWESS Smoothed Transmissibility', ...
         'Apparent Theoretical Fitting Curve', ...
         'Identified Natural Frequency'}, ...
         'Location','best');

else

    legend([h1,h2,h4], ...
        {'Actual Experimental Transmissibility', ...
         'LOWESS Smoothed Transmissibility', ...
         'Identified Natural Frequency'}, ...
         'Location','best');

end

grid on;
box on;

xlim([0 max_frequency]);

set(gca,'FontSize',11,'LineWidth',1);

hold off;


    %% =====================================================
    % 14. PLOT FREQUENCY-DEPENDENT DAMPING
    % ======================================================

    figure('Color','w');

    plot(freq_smooth,zeta_freq, ...
        'Color',[0.7 0.7 0.7], ...
        'LineWidth',1);

    hold on;

    plot(freq_smooth,zeta_smooth, ...
        'r','LineWidth',1.8);

    xline(fn,'b--', ...
        sprintf('Natural frequency = %.2f Hz',fn), ...
        'LineWidth',1.5);

    xlabel('Frequency (Hz)');

    ylabel('Equivalent Damping Ratio, \zeta');

    title(['Frequency-Dependent Damping - ',group_name]);

    legend('Calculated Damping', ...
        'LOWESS Smoothed Damping', ...
        'Natural Frequency', ...
        'Location','best');

    grid on;
    xlim([0 max_frequency]);

    %% =====================================================
    % 15. EXPORT INDIVIDUAL RESULTS TO EXCEL
    % ======================================================

    Results = table( ...
        freq_smooth, ...
        TR_smooth, ...
        TR_smooth_dB, ...
        zeta_freq, ...
        zeta_smooth, ...
        'VariableNames',{ ...
        'Frequency_Hz', ...
        'Transmissibility', ...
        'Transmissibility_dB', ...
        'Damping_Ratio', ...
        'Smoothed_Damping_Ratio'});

    output_excel = [group_name,'_Damping_Results.xlsx'];

    writetable(Results,output_excel);

    fprintf('Results saved to %s\n',output_excel);

    %% =====================================================
    % 16. STORE SUMMARY RESULTS
    % ======================================================

    NewRow = table( ...
        string(group_name), ...
        fn, ...
        damping_at_fn, ...
        'VariableNames',{ ...
        'Group', ...
        'Natural_Frequency_Hz', ...
        'Damping_at_Natural_Frequency'});

    Summary = [Summary;NewRow];

end

%% =========================================================
% 17. FINAL SUMMARY
% ==========================================================

fprintf('\n====================================\n');
fprintf('FINAL DAMPING IDENTIFICATION RESULTS\n');
fprintf('====================================\n');

disp(Summary);

writetable(Summary,'Final_Damping_Summary.xlsx');

fprintf('\nAll processing completed.\n');
fprintf('Final summary saved to Final_Damping_Summary.xlsx\n');
%% =========================================================
% 18. THEORETICAL TRANSMISSIBILITY MODEL
% ==========================================================

function T_fit_dB = transmissibility_model(f, fn, zeta)

    r = f ./ fn;

    num = sqrt(1 + (2*zeta.*r).^2);

    den = sqrt((1-r.^2).^2 + (2*zeta.*r).^2);

    T_abs = num ./ den;

    T_fit_dB = 20*log10(T_abs + 1e-12);

end
