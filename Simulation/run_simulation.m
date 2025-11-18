clear; clc;

% 1. Define Simulation Parameters
M = 4;                  % QPSK
k = log2(M);            % Bits per symbol
SymbolRate = 1000;      % 1 kSps
numBits = 10000;        % Bits per frame

% 2. Define SNR Range (Eb/No)
EbNo_Range = 0:1:10;    % Test from 0dB to 10dB
BER_Results = zeros(size(EbNo_Range));

% 3. Loop through SNRs
disp('Starting Simulation Loop...');
for i = 1:length(EbNo_Range)
    
    % Update the variable used inside the Simulink AWGN block
    EbNo_val = EbNo_Range(i); 
    
    % Run the simulation (Replace 'Task1_Baseband_Chain' with your filename)
    simOut = sim('Task1'); 
    
    % Extract the BER (Assumes "To Workspace" block is named "BER_Result")
    % The Error Rate block outputs a 1x3 vector: [BER, TotalErrors, TotalBits]
    % We take the last value recorded
    BER_Results(i) = simOut.BER_Result(end, 1);
    
    fprintf('Eb/No: %d dB -> BER: %e\n', EbNo_val, BER_Results(i));
end

% 4. Plot Results (Required by Task II)
figure;
semilogy(EbNo_Range, BER_Results, '-o');
grid on;
xlabel('Eb/No (dB)');
ylabel('Bit Error Rate (BER)');
title('QPSK Performance: BER vs Eb/No');