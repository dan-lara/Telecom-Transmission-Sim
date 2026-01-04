clc; clear; close all;

%% ==========================================
%% 1. CONFIGURATION OF PARAMETERS FOR WEBLAB
%% ==========================================

% Fixed sampling frequency of WebLab (200 MHz)
fs_weblab = 200e6;

% Signal length: 10,000 samples (sufficient for DPD)
target_length = 10000;

% Choice of modulation - 16QAM recommended for DPD characterization
modulation_type = '16QAM';

% RRC (Root Raised Cosine) filter parameters
rolloff = 0.35;      % Roll-off factor (standard LTE)
span = 10;           % Filter span in symbols
sps = 16;            % Samples per symbol (200 MHz / 8 = 25 MHz)

%% ==========================================
%% 2. DIGITAL SIGNAL GENERATION
%% ==========================================

% 2.1 Configuration of modulation parameters
switch upper(modulation_type)
    case 'BPSK'
        M = 2;  k = 1;   % 2 states, 1 bit per symbol
    case 'QPSK'
        M = 4;  k = 2;   % 4 states, 2 bits per symbol
    case '16QAM'
        M = 16; k = 4;   % 16 states, 4 bits per symbol
    otherwise
        error('Modulation not supported. Use BPSK, QPSK or 16QAM.');
end

% 2.2 Calculate the number of bits needed
% We add 2*span to compensate for edge effects of filtering
num_symbols_raw = ceil(target_length / sps) + 2*span;
num_bits_needed = num_symbols_raw * k;

% 2.3 Generation of random bits
fprintf('Generating %d random bits...\n', num_bits_needed);
tx_bits = randi([0 1], num_bits_needed, 1);

% 2.4 Baseband modulation
fprintf('%s modulation in progress...\n', modulation_type);
tx_symbols = qammod(tx_bits, M, 'InputType', 'bit', 'UnitAveragePower', true);

% Check for NaN or Inf
if any(isnan(tx_symbols(:))) || any(isinf(tx_symbols(:)))
    error('ERROR: Modulation produced NaN or Inf values!');
end

% 2.5 Design of RRC filter
fprintf('Designing RRC filter (rolloff=%.2f, span=%d, sps=%d)...\n', rolloff, span, sps);
h_rrc = rcosdesign(rolloff, span, sps, 'sqrt');

% 2.6 Pulse shaping
fprintf('RRC filtering in progress...\n');
tx_filtered = upfirdn(tx_symbols, h_rrc, sps);

% 2.7 Extraction of stable central part
% We avoid transients at the beginning and end of the filtered signal
fprintf('Extracting stable central part...\n');
center_idx = floor(length(tx_filtered) / 2);
half_len = floor(target_length / 2);

start_idx = max(1, center_idx - half_len + 1);
end_idx = min(length(tx_filtered), center_idx + half_len);

PAin = tx_filtered(start_idx:end_idx);

% 2.8 Check final length
if length(PAin) ~= target_length
    fprintf('WARNING: Length adjusted from %d to %d samples\n', target_length, length(PAin));
    target_length = length(PAin);
end

%% ==========================================
%% 3. SIGNAL PREPARATION AND NORMALIZATION
%% ==========================================

% 3.1 Calculate signal statistics
signal_power = mean(abs(PAin).^2);
signal_peak = max(abs(PAin));
PAPR_linear = signal_peak^2 / signal_power;
PAPR_dB = 10 * log10(PAPR_linear);

fprintf('--- Raw signal statistics ---\n');
fprintf('Average power : %.4f\n', signal_power);
fprintf('Peak amplitude: %.4f\n', signal_peak);
fprintf('PAPR          : %.2f dB\n', PAPR_dB);

% 3.2 CRITICAL normalization for WebLab
% WebLab expects a signal with normalized peak amplitude
% We use a safety factor of 0.9 to avoid saturation
safety_factor = 0.9;
PAin_normalized = PAin / signal_peak * safety_factor;

% 3.3 Final value verification
if any(isnan(PAin_normalized(:)))
    error('ERROR: PAin contains NaN values after normalization!');
end

if any(isinf(PAin_normalized(:)))
    error('ERROR: PAin contains Inf values after normalization!');
end

% 3.4 Format to column vector
PAin = double(PAin_normalized(:));

% Convert to complex if necessary
if isreal(PAin)
    PAin = complex(PAin, zeros(size(PAin)));
end

% 3.5 Calculate final statistics
final_power = mean(abs(PAin).^2);
final_peak = max(abs(PAin));
final_PAPR_dB = 10 * log10(final_peak^2 / final_power);

fprintf('\n--- Statistics after normalization ---\n');
fprintf('Peak amplitude  : %.4f (safety limit = %.2f)\n', final_peak, safety_factor);
fprintf('Average power   : %.6f\n', final_power);
fprintf('Final PAPR      : %.2f dB\n', final_PAPR_dB);

%% ==========================================
%% 4. CONFIGURATION OF VARIABLES FOR MAIN.M
%% ==========================================

% 4.1 Sampling frequency (required by main.m)
Fs = fs_weblab;

% 4.2 Signal bandwidth
BW = fs_weblab / sps;

% 4.3 ACPR parameters for linearity calculation
% These values correspond to the LTE 20MHz standard
ACPR.BW = 18e6;       % Useful channel bandwidth
ACPR.Offset = 20e6;   % Adjacent channel offset

% 4.4 Additional information for debugging
signal_info.modulation = modulation_type;
signal_info.original_length = length(tx_filtered);
signal_info.final_length = length(PAin);
signal_info.sps = sps;
signal_info.rolloff = rolloff;
signal_info.span = span;
signal_info.PAPR_dB = final_PAPR_dB;
signal_info.peak_amplitude = final_peak;
signal_info.rms_amplitude = sqrt(final_power);

% 4.5 Compatibility check with main.m
fprintf('\n--- Compatibility check with main.m ---\n');

% Check that PAin is complex
if isreal(PAin)
    fprintf('WARNING: PAin is a real signal, converting to complex...\n');
    PAin = complex(PAin, zeros(size(PAin)));
end

% Calculate PAPR as in main.m
PAPRin = 10*log10(max(abs(PAin).^2) / mean(abs(PAin).^2));
fprintf('PAPR calculated (main.m method): %.2f dB\n', PAPRin);

% Calculate RMSin as in main.m (for reference)
RMSin_calculated = -8.5 - PAPRin - 2;
fprintf('Estimated RMSin (for main.m): %.2f dBm\n', RMSin_calculated);

%% ==========================================
%% 5. DATA SAVING
%% ==========================================

% 5.1 File name (MUST be 'PAinLTE20MHz.mat' for main.m)
filename = 'PAinLTE20MHz.mat';

% 5.2 List of variables to save
variables_to_save = {'PAin', 'Fs', 'BW', 'ACPR', 'tx_bits', 'signal_info'};

% 5.3 Save with verification
fprintf('\n--- Saving data ---\n');
fprintf('File: %s\n', filename);

try
    % Save essential variables
    save(filename, variables_to_save{:});
    
    % Verify that the file can be loaded
    loaded_data = load(filename);
    
    % Check each variable
    for i = 1:length(variables_to_save)
        var_name = variables_to_save{i};
        if ~isfield(loaded_data, var_name)
            error('Variable %s missing in saved file', var_name);
        end
    end
    
    fprintf('✅ File saved successfully\n');
    
catch ME
    fprintf('❌ ERROR during save: %s\n', ME.message);
    rethrow(ME);
end

%% ==========================================
%% 6. FINAL SUMMARY AND INSTRUCTIONS
%% ==========================================

fprintf('\n==========================================\n');
fprintf('GENERATION COMPLETED - FINAL SUMMARY\n');
fprintf('==========================================\n');

fprintf('Generated signal: %s\n', modulation_type);
fprintf('Samples        : %d (at %d MHz)\n', length(PAin), Fs/1e6);
fprintf('Duration       : %.3f ms\n', length(PAin)/Fs*1000);
fprintf('Bandwidth      : %.1f MHz\n', BW/1e6);
fprintf('Peak amplitude : %.4f\n', max(abs(PAin)));
fprintf('PAPR           : %.2f dB\n', PAPRin);

fprintf('\nVariables saved in %s:\n', filename);
for i = 1:length(variables_to_save)
    fprintf('  • %s\n', variables_to_save{i});
end

fprintf('\n==========================================\n');
fprintf('INSTRUCTIONS FOR USE WITH WEBLAB\n');
fprintf('==========================================\n');

fprintf('1. Ensure that %s is in the current directory\n', filename);
fprintf('2. Run main.m to launch PA characterization\n');
fprintf('3. The following steps will be executed:\n');
fprintf('   a) Measurement of AM-AM/AM-PM characteristics\n');
fprintf('   b) Calculation of ACPR and EVM\n');
fprintf('   c) Identification of DPD coefficients\n');
fprintf('   d) Application of digital predistortion\n');

% Final verification
fprintf('\n=== FINAL VERIFICATION ===\n');

% Check for absence of NaN/Inf
if any(isnan(PAin(:))) || any(isinf(PAin(:)))
    fprintf('❌ ERROR: PAin contains non-finite values!\n');
else
    fprintf('✅ PAin contains no NaN or Inf\n');
end

% Power verification
if final_power < 1e-6
    fprintf('⚠️  WARNING: Very low power (%.6f)\n', final_power);
else
    fprintf('✅ Power in a reasonable range\n');
end

% Peak amplitude verification
if final_peak > 1.0
    fprintf('❌ ERROR: Peak amplitude > 1.0 (risk of saturation)\n');
elseif final_peak < 0.1
    fprintf('⚠️  WARNING: Very low peak amplitude (%.4f)\n', final_peak);
else
    fprintf('✅ Correct peak amplitude (%.4f)\n', final_peak);
end

fprintf('\n=== READY FOR EXECUTION WITH WEBLAB ===\n');