%% Digital Predistortion Linearization of PA in WebLab
% 
% 20/11/2019

clear all;close all;clc;

%% Frequency Object
Window='Blackman-Harris';%'blackman';
SegLength=2^10;
overlap=25;
Hs = spectrum.welch(Window,SegLength,overlap);

FuncMode='Testbench'; %DataAcquis Testbench ILC
%% Display control
EnAMbL='on';
EnSpec='on';
% EnTime='off';

%% GetStimulus
load('PAinLTE20MHz') % load your own baseband data

%% Get PA model

PAPRin=papr(PAin);
RMSin=-8.5-PAPRin-2;
[PAout, RMSout, Idc, Vdc] = RFWebLab_PA_meas_v1_2(PAin, RMSin); 
PAout=timealign(PAin,PAout);
        DisplayOptions.InterCorrArg.Enable='off';
        DisplayOptions.IQ.Scope=1e5;
        DisplayOptions.IQ.Enable='off';
        DisplayOptions.NMSE.Enable='on';
    Data.In=PAin;
    Data.Out=PAout;
    Data.ALimLinIn=0.2;
    
    PA=amampm(Data,EnAMbL);
%% Predistortion identification: Applying IBO

% Backoff=20*log10(PA.LimitPD/max(abs(PAin)));
% 
% % Set your input average power
% RMSin=RMSin+Backoff;
fprintf('\n=== Backoff Correction Calculation ===\n');

% Check and fix PA.LimitPD
if isnan(PA.LimitPD)
    fprintf('Fix: PA.LimitPD = NaN, using estimated InSat\n');
    if ~isnan(PA.InSat) && PA.InSat > 0
        PA.LimitPD = PA.InSat * 0.85;  % 85% saturation point
    else
        PA.LimitPD = 0.8;  % default value
    end
    fprintf('  Set PA.LimitPD = %.4f\n', PA.LimitPD);
end

% Calculate Backoff
current_peak = max(abs(PAin));
saturation_margin = PA.InSat / current_peak;  % Saturation margin

fprintf('\n=== Calculate Safe Backoff ===\n');
fprintf('Current peak: %.4f\n', current_peak);
fprintf('Saturation point (InSat): %.4f\n', PA.InSat);
fprintf('Saturation margin: %.4f (%.2f dB)\n', saturation_margin, 20*log10(saturation_margin));

% Calculate backoff to ensure sufficient margin
required_margin = 1.3;  % Increase margin to 30%
if saturation_margin < required_margin
    % Calculate required backoff: ensure 30% margin
    required_backoff = 20 * log10(required_margin / saturation_margin);
    fprintf('Signal close to saturation, backoff needed: %.2f dB\n', required_backoff);
    
    % Ensure backoff is positive (reduce power)
    Backoff = -abs(required_backoff);
else
    % If sufficient margin, use zero backoff
    Backoff = 0;
    fprintf('Signal has sufficient margin, using zero backoff\n');
end

% Add additional power limit
max_backoff = -10;  % Maximum backoff limit
if Backoff < max_backoff
    fprintf('Warning: Backoff too large (%.2f dB), limited to %.2f dB\n', Backoff, max_backoff);
    Backoff = max_backoff;
end

fprintf('Apply backoff: %.2f dB\n', Backoff);
fprintf('Update RMSin: %.2f dBm + %.2f dB = ', RMSin, Backoff);
RMSin = RMSin + Backoff;
fprintf('%.2f dBm\n', RMSin);

[PAout, RMSout, Idc, Vdc] = RFWebLab_PA_meas_v1_2(PAin, RMSin);

PAout=timealign(PAin,PAout);
%% Input/Output Pwr spectral density 
[ACPRin, PSDin]=acpr(PAin,Fs,ACPR);
[ACPRout, PSDout]=acpr(PAout,Fs,ACPR);
X=PSDin.Data;
PXin=10*log10(X);
Y=PSDout.Data;
PYout=10*log10(Y);
f=PSDout.Frequencies;

fprintf(['\n\t The output power before linearization is',...
    ' equal to %.2f dBm\n'],...
    RMSout)


%% Model PD first stage
clear ModelPD

PD.In=PAin;
PD.BW=BW;
%% Iterations
section=10e3; % number of samples for each iteration in indirect learning
N=length(PAin)-rem(length(PAin),section);

clear Data
Data.Fs=Fs;
Data.ACPR=ACPR;

Algorithm.NbIterPerStage=4;
Algorithm.DampingNewtonFactor=0.7;
Algorithm.NbSampPerIter=section;

Stimulus.wf=PAin;
Stimulus.RMSin=RMSin;
Stimulus.Fs=Fs;
Stimulus.ACPR=ACPR;

%% Launch DPD Test
tic;

SystIter=DPD(Stimulus,ACPRout,Algorithm,PD); % change in file "DPD.m" to implement your DPD algorithm

timerVal=toc;
fprintf(['\n\t Execution time: ', num2str(timerVal), ' seconds'])

% Define num_iter before using it
num_iter = length(SystIter);

fprintf('\nACPR Change Trend:\n');
fprintf('Iter | ACPR-1 (dB) | ACPR+1 (dB) | ACPR-2 (dB) | ACPR+2 (dB) | Improvement (dB)\n');
fprintf('-----|-------------|-------------|-------------|-------------|-------------------\n');

for i = 1:num_iter
    if isfield(SystIter(i).ACPR, 'L1')
        % Use L1/U1/L2/U2 fields
        if i == 1
            fprintf('%4d | %11.2f | %11.2f | %11.2f | %11.2f | %10s\n', ...
                i, SystIter(i).ACPR.L1, SystIter(i).ACPR.U1, ...
                SystIter(i).ACPR.L2, SystIter(i).ACPR.U2, 'Baseline');
        else
            improv_1 = mean([SystIter(i).ACPR.L1, SystIter(i).ACPR.U1]) - ...
                      mean([SystIter(1).ACPR.L1, SystIter(1).ACPR.U1]);
            improv_2 = mean([SystIter(i).ACPR.L2, SystIter(i).ACPR.U2]) - ...
                      mean([SystIter(1).ACPR.L2, SystIter(1).ACPR.U2]);
            fprintf('%4d | %11.2f | %11.2f | %11.2f | %11.2f | %+6.2f (Avg)\n', ...
                i, SystIter(i).ACPR.L1, SystIter(i).ACPR.U1, ...
                SystIter(i).ACPR.L2, SystIter(i).ACPR.U2, (improv_1+improv_2)/2);
        end
    else
        % Use L/R fields
        if i == 1
            fprintf('%4d | %11.2f | %11.2f | %11.2f | %11.2f | %10s\n', ...
                i, SystIter(i).ACPR.L(1), SystIter(i).ACPR.R(1), ...
                SystIter(i).ACPR.L(2), SystIter(i).ACPR.R(2), 'Baseline');
        else
            improv_1 = mean([SystIter(i).ACPR.L(1), SystIter(i).ACPR.R(1)]) - ...
                      mean([SystIter(1).ACPR.L(1), SystIter(1).ACPR.R(1)]);
            improv_2 = mean([SystIter(i).ACPR.L(2), SystIter(i).ACPR.R(2)]) - ...
                      mean([SystIter(1).ACPR.L(2), SystIter(1).ACPR.R(2)]);
            fprintf('%4d | %11.2f | %11.2f | %11.2f | %11.2f | %+6.2f (Avg)\n', ...
                i, SystIter(i).ACPR.L(1), SystIter(i).ACPR.R(1), ...
                SystIter(i).ACPR.L(2), SystIter(i).ACPR.R(2), (improv_1+improv_2)/2);
        end
    end
end

% Calculate final improvement
if isfield(SystIter(end).ACPR, 'L1')
    final_improv_1 = mean([SystIter(end).ACPR.L1, SystIter(end).ACPR.U1]) - ...
                     mean([SystIter(1).ACPR.L1, SystIter(1).ACPR.U1]);
    final_improv_2 = mean([SystIter(end).ACPR.L2, SystIter(end).ACPR.U2]) - ...
                     mean([SystIter(1).ACPR.L2, SystIter(1).ACPR.U2]);
else
    final_improv_1 = mean([SystIter(end).ACPR.L(1), SystIter(end).ACPR.R(1)]) - ...
                     mean([SystIter(1).ACPR.L(1), SystIter(1).ACPR.R(1)]);
    final_improv_2 = mean([SystIter(end).ACPR.L(2), SystIter(end).ACPR.R(2)]) - ...
                     mean([SystIter(1).ACPR.L(2), SystIter(1).ACPR.R(2)]);
end

fprintf('\nSummary:\n');
fprintf('  ACPR1 (±20MHz) Improvement: %.2f dB\n', final_improv_1);
fprintf('  ACPR2 (±40MHz) Improvement: %.2f dB\n', final_improv_2);

DisplayFigures;

%% Create comprehensive performance analysis plot

figure('Name', 'DPD Linearization Performance Analysis', 'Position', [100, 100, 1200, 800]);

% Subplot 1: ACPR vs Iteration
subplot(2, 3, 1);
hold on;

% === Data extraction and preprocessing ===
% Extract ACPR data from structure array
acpr_data = [SystIter.ACPR]; 

% Check which field structure is available and extract data accordingly
if isfield(acpr_data, 'L1')
    % Use L1/U1/L2/U2 fields
    L1 = [acpr_data.L1];
    L2 = [acpr_data.L2];
    U1 = [acpr_data.U1];
    U2 = [acpr_data.U2];
    
    % Combine into matrices
    L_mat = [L1; L2];
    R_mat = [U1; U2];
elseif isfield(acpr_data, 'L')
    % Use L/R fields
    L_mat = [acpr_data.L];
    R_mat = [acpr_data.R];
else
    error('ACPR data structure not recognized. Check available fields.');
end

% === Plotting ===
% Plot L(1) - Lower ACPR (1st adjacent)
plot(1:num_iter, L_mat(1, :), 'b-o', 'LineWidth', 2, 'MarkerSize', 8);

% Plot R(1) - Upper ACPR (1st adjacent)
plot(1:num_iter, R_mat(1, :), 'b--s', 'LineWidth', 1.5, 'MarkerSize', 6);

% Plot L(2) - Lower ACPR (2nd adjacent)
plot(1:num_iter, L_mat(2, :), 'r-o', 'LineWidth', 2, 'MarkerSize', 8);

% Plot R(2) - Upper ACPR (2nd adjacent)
plot(1:num_iter, R_mat(2, :), 'r--s', 'LineWidth', 1.5, 'MarkerSize', 6);

hold off; grid on; 
xlabel('Iteration Number'); 
ylabel('ACPR (dBc)');
title('ACPR Change with DPD Iterations');
legend('ACPR-1 (Left)', 'ACPR+1 (Right)', 'ACPR-2 (Left)', 'ACPR+2 (Right)', 'Location', 'best');

% Subplot 2: Power efficiency
subplot(2, 3, 2);
pout_dbm = [SystIter.RMSout];
pout_w = 10.^(pout_dbm/10 - 3);  % Convert dBm to watts
dc_power = [SystIter.Vdc] .* [SystIter.Idc];
efficiency = pout_w ./ dc_power * 100;
plot(1:num_iter, efficiency, 'g-^', 'LineWidth', 2, 'MarkerSize', 8);
grid on; 
xlabel('Iteration Number'); 
ylabel('Efficiency (%)');
title('PA Efficiency Change');
ylim([0, max(efficiency)*1.2]);

% Subplot 3: DC power consumption
subplot(2, 3, 3);
plot(1:num_iter, [SystIter.Idc], 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
plot(1:num_iter, [SystIter.Vdc], 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
grid on; 
xlabel('Iteration Number'); 
ylabel('Value');
title('DC Current and Voltage');
legend('Idc (A)', 'Vdc (V)', 'Location', 'best');

% Subplot 4: Spectrum comparison (last iteration)
subplot(2, 3, [4, 5, 6]);
hold on;
% Plot spectrum before and after DPD
plot(PSDout.Frequencies/1e6, 10*log10(PSDout.Data), 'b-', 'LineWidth', 2);
plot(PSDout.Frequencies/1e6, SystIter(end).PY, 'r-', 'LineWidth', 2);
grid on; 
xlabel('Frequency (MHz)'); 
ylabel('Power Spectral Density (dB)');
title('Spectrum Comparison Before and After DPD');
legend('Before DPD', 'After DPD', 'Location', 'best');

% Mark ACPR measurement points
xline(-20, 'k--', 'LineWidth', 1, 'Alpha', 0.5);
xline(20, 'k--', 'LineWidth', 1, 'Alpha', 0.5);
xline(-40, 'k--', 'LineWidth', 1, 'Alpha', 0.5);
xline(40, 'k--', 'LineWidth', 1, 'Alpha', 0.5);