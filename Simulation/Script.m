% This code is for 
%%   
clc;
clear;
close all;

%% Configuration
% Modulation type
MODULATION_TYPE = 'QPSK';    % This script provide 'BPSK', 'QPSK' and '16QAM'

% Automatically set k and M according to the modulation type
% M = 2^k => k = log2(M)
switch MODULATION_TYPE
    case 'BPSK'
        M = 2;
        k = 1;
    case 'QPSK'
        M = 4;
        k = 2;
    case '16QAM'
        M = 16;
        k = 4;
    otherwise
        error('Unknown modulation type');
end

% Basic parameters
N = 1e5;        % Numbre of bits
fc = 5e3;       % carrier frequency
Rb = 1e3;       % Symbol rate
Rs = Rb*k;      % Bit rate
sps = 32;        % Samples Per Symbol
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
    tx_symbols = pskmod(tx_bits, M, pi/4, 'InputType','bit');
elseif (M == 16)
    tx_symbols = qammod(tx_bits, M, 'InputType','bit');
    %tx_symbols = qammod(tx_bits, M, 'InputType', 'bit', 'UnitAveragePower', true);
end

% filter
rolloff = 0.35; % rolloff factor for the Root Raised Cosine filter
span = 6;       % Filter length
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

%% Mixer - Down convension
% Down convention
rx_baseband = rx_rf .* (cos(2*pi*fc*t') - 1j*sin(2*pi*fc*t'));

% Matched Filter
rx_filtered = upfirdn(rx_baseband, h_rrc, 1, 1);

% Downsampling and Delay Compensation
delay = span * sps; 
rx_symbols = rx_filtered(delay + 1 : sps : end);

%% Receiver - Rx
% Demodulation
if (M == 2 || M ==4)
    rx_bits = pskdemod(rx_symbols, M, pi/4, 'OutputType', 'bit');
elseif (M == 16)
    rx_bits = qamdemod(rx_symbols, M, 'OutputType', 'bit');
end

% BER Calculation
[numErrors, BER] = biterr(tx_bits, rx_bits);

%% Resulte
disp(['Bit Error Rate: ', num2str(BER)]);
if BER == 0
    disp('OK');
else
    disp('Bit errors');
end
scatterplot(rx_symbols); title('No noise received signal constellation diagram'); %



