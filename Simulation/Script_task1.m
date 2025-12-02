% This code is for simulating the telecommunition process in a script
%%   
clc;
clear;
close all;

%% Configuration
% Modulation type
MODULATION_TYPE = '16QAM';    % This script provide 'BPSK', 'QPSK' and '16QAM'

% Automatically set k and M according to the modulation type
% M = 2^k => k = log2(M)
switch MODULATION_TYPE
    case 'BPSK'
        M = 2; k = 1; phase_offset = 0;
    case 'QPSK'
        M = 4; k = 2; phase_offset = pi/4;
    case '16QAM'
        M = 16;k = 4; phase_offset = 0;
    otherwise
        error('Unknown modulation type');
end

% Basic parameters
N = 1e6;        % Numbre of bits
N = floor(N / k) * k;   % Ensure the number of bits is an integer multiple of k
fc = 5e3;       % carrier frequency
Rb = 1e3;       % Symbol rate
Rs = Rb*k;      % Bit rate
sps = 32;       % Samples Per Symbol
fs = Rb*sps;    % Sample frequency = symbol rate * the numbre of sample point for each symbol
n_total = N*sps;        % Sample point total
n_each = N/2*sps;       % Sample point for each bit
ts = 1/fs;      % Sample time
Tb = ts*sps*2;  % Code period
t1 = (0:N-1)*Tb/2;      % Code duration
t2 = (0:N-1)*ts;        % Total duration of codde    

%% Transmitter - Tx
% Generate binary data source
tx_bits = randi([0 1], N, 1);

% Map binary bits to constellation points on the complex plane
if (M == 2 || M ==4)
    tx_symbols = pskmod(tx_bits, M, phase_offset, 'InputType','bit');
elseif (M == 16)
    tx_symbols = qammod(tx_bits, M, 'InputType', 'bit', 'UnitAveragePower', true);
end

% filter
rolloff = 0.35; % rolloff factor for the Root Raised Cosine filter
span = 10;       % Filter length
h_rrc = rcosdesign(rolloff, span, sps, "sqrt");
tx_filtered = upfirdn(tx_symbols, h_rrc, sps);

%% Mixer - Up convension
len_sigal = length(tx_filtered);
t = (0 : len_sigal - 1) / fs;

% I et Q signal
I_signal = real(tx_filtered);
Q_signal = imag(tx_filtered);

% IQ Modulation
tx_rf = I_signal .* cos(2*pi*fc*t') - Q_signal .* sin(2*pi*fc*t');

%% Channel
rx_rf = tx_rf;
% rx_rf = awgn(tx_rf, snr, 'measured');

%% Mixer - Down convension
% Down convention
rx_baseband = rx_rf .* (cos(2*pi*fc*t') - 1j*sin(2*pi*fc*t'))* 2;

% Matched Filter
rx_filtered = upfirdn(rx_baseband, h_rrc, 1, 1);

% Downsampling and Delay Compensation
total_delay = span * sps; % Total delay (TX filter + RX filter)
rx_symbols = rx_filtered(total_delay + 1 : sps : end);

% Upfirdn
% Truncate redundant tail data
num_tx_symbols = length(tx_symbols);
rx_symbols = rx_symbols(1:num_tx_symbols);


%% Receiver - Rx
% Demodulation
if (M == 2 || M ==4)
    rx_bits = pskdemod(rx_symbols, M, phase_offset, 'OutputType', 'bit');
elseif (M == 16)
    rx_bits = qamdemod(rx_symbols, M, 'OutputType', 'bit', 'UnitAveragePower', true);
end

% BER Calculation
[numErrors, BER] = biterr(tx_bits, rx_bits);

%% Result
disp(['Modulation: ', MODULATION_TYPE]);
disp(['Bit Error Rate: ', num2str(BER)]);

scatterplot(tx_symbols); 
title('No noise transmitted signal constellation diagram'); 
grid on;

% Power/EVM
rx_power = mean(abs(rx_symbols).^2);
disp(['Received Signal Average Power: ', num2str(rx_power)]);
% In theory, after normalization, it should be close to 1.

%% Visualization
% Remove the transients from the head and the attenuation from the tail
% We only take the stable segment in the middle.
clean_signal = rx_filtered(total_delay + 1 : end - total_delay);

% Eye diagram
% 'period', 2*sps indicates that each chart displays 2 symbol periods.
eyediagram(clean_signal(1:2000), 2*sps); 
title('Optimized Eye Diagram (Transients Removed)');
% eyediagram(rx_filtered(1:2000), 2*sps);

scatterplot(rx_symbols); 
title('No noise received signal constellation diagram'); 
grid on;

L = length(tx_filtered);    % Signal length
f = (-L/2 : L/2-1)*(fs/L);  % Frequency axis (Hz)

% Frequency Domain (Spectrum)
% 1. Calculate the baseband signal spectrum (before transmission)v
spec_baseband = fftshift(fft(tx_filtered));
power_baseband = abs(spec_baseband).^2/L;

% 2. Calculate the radio frequency signal spectrum (after mixing)
spec_rf = fftshift(fft(tx_rf));
power_rf = abs(spec_rf).^2/L;

figure('Name', 'Spectrum Analysis');
max_p = max(power_baseband);
subplot(2,1,1);
plot(f, 10*log10(power_baseband));
% plot(f, 10*log10(power_baseband/max_p)); % Normalized to 0dB peak
grid on;
title(['Baseband Spectrum (' MODULATION_TYPE ')']);
xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([-fc*2, fc*2]); % Limit the display area for observation

subplot(2,1,2);
plot(f, 10*log10(power_rf));
grid on;
title(['RF Passband Spectrum (Carrier = ' num2str(fc) 'Hz)']);
xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([-fc*2, fc*2]);

% Time Domain
figure('Name', 'Time Domain Waveform');
if length(t) > 300
    plot_range = 100:300;
else
    plot_range = 1:length(t);
end
t_plot = t(plot_range);

subplot(2,1,1);
plot(t_plot, real(tx_filtered(plot_range)), 'LineWidth', 1.5);
hold on;
plot(t_plot, imag(tx_filtered(plot_range)), '--');
grid on;
legend('I (In-phase)', 'Q (Quadrature)');
title('Baseband Signal (Time Domain)');
xlabel('Time (s)');

subplot(2,1,2);
plot(t_plot, tx_rf(plot_range));
grid on;
title('RF Modulated Signal (After Mixer)');
xlabel('Time (s)');
ylabel('Amplitude');
