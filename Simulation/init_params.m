% Daniel FERREIRA LARA
% 18/11/2025
% Parameters for Task I & II
clear; clc;

M = 4;                  % Modulation Order (4 for QPSK) 
k = log2(M);            % Bits per symbol
Fs = 200e3;             % Sample rate (Hz) - Arbitrary for sim, but 200kHz is typical
Ts = 1/Fs;              % Sample time
Sps = 8;                % Samples per symbol (Oversampling factor)
SymbolRate = Fs/Sps;    % Symbol rate
EbNo_val = 10;          %
% Reporting
disp('Parameters loaded. Ready to run Simulink.');