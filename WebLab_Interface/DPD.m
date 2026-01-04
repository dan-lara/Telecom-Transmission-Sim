function SystIter=DPD(Stimulus,ACPRout,Algorithm,PD)
% SystIter=DPD(Stimulus,ACPRout,Algorithm,PD)
% This function implements DPD algorithms.
% It takes the following input arguments:
% Stimulus		: Structure with  three fields
%					wf	: for waveform, an Nx1 vector, input samples
%						additional code is needed to take the case 1xN)
%					Fs	: sampling frequency
%					ACPR: structure with 2 fields (BW: for bandwidth and Offset)
% ACPRout		: Structure with ACPR measurements of PA output; To avoid
%				computing such measurement inside this function, and used
%				as a refernces to compute ACPR improvement after predistortion
%				(Can be implemented in a better way / to be modified later)
% Algorithm		: Structure with three fields (in this function)
%                   Method : String to choose the method to use:
%                           'Backward' for the first algorithm.
%                           'Alternate' for the second one.
%					NbIter	or NbIterPerStage : Number of iteration
%                           (positive integer, verification
%							should be added later)
%                           NbIterPerStage when using 'Backward' method
%							NbIter when using 'Alternate' method
%					NbSampPerIter	: Number of samples to be used for each system
%										level iteration (positive integer, verification needed)
% PD			: Structure of a general parallel structure nonlinear model,
%				Each stage has its own number of iterations
%
% 	(c) Siqi WANG, 2019

RMSin=Stimulus.RMSin;
PAin=fixpwr(PD.In,RMSin);
Fs=Stimulus.Fs;
ACPR=Stimulus.ACPR;

clear PD.Model

% configure your DPD model parameters here
        PD.Model.NL(1).Kv=[0:6]; 
        PD.Model.Mem(1).Lv=[0:3];
        PD.Model.NL(2).Kv=[1:4];
        PD.Model.Mem(2).Lv=[0:3];
        PD.Model.Mem(2).Mv=[1:2];
        PD.Model.NL(3).Kv=[1:4];
        PD.Model.Mem(3).Lv=[0:3];
        PD.Model.Mem(3).Mv=[1:2];
    PD.Model.bias='off';
    PD.Model.Symmetric='off';
    PD.Model.Type='gmp';
    NbCoeff=nbcoefgmp(PD.Model);
    PD.Model.Coeff=zeros(NbCoeff,1);
    PD.Model.Coeff(1)=1;

section=Algorithm.NbSampPerIter;
mu=Algorithm.DampingNewtonFactor;
Data.Fs=Fs;
Data.ACPR=ACPR;

% Add regularization parameters to avoid matrix singularity
Options.Evaluation.Enable='on';
Options.Evaluation.Regularization = 1e-6;  % Add regularization
Options.Evaluation.MaxCondition = 1e6;     % Maximum condition number limit

SystIter(1).ACPR0=ACPRout;
NbIter=Algorithm.NbIterPerStage;
for n=1:NbIter
    fprintf('\n=== DPD Iteration %d/%d ===\n', n, NbIter);
    
    u=PAin(1:section);
    SystIter(n).u=u;
    
    if n>1
        PD.x=gmp(u,PD.Model);
        
        % Check PD.x amplitude to prevent overshoot
        pd_peak = max(abs(PD.x));
        fprintf('PD output peak: %.4f\n', pd_peak);
        
        if pd_peak > 1.2  % If PD output is too large
            fprintf('Warning: PD output too large (%.4f), applying clipping\n', pd_peak);
            PD.x = PD.x / pd_peak * 0.9;  % Clip to 0.9
        end
    else
        PD.x=u;
        figure('DefaultAxesFontSize',20,'WindowStyle','docked','NumberTitle','off');
    end
    
    % Check input power
    input_power_dbm = 10*log10(mean(abs(PD.x).^2)*1000);
    fprintf('Input power: %.2f dBm\n', input_power_dbm);
    
    % Set reasonable power range: -20 to -10 dBm
    target_power_dbm = -15;  % Target power
    
    if input_power_dbm > target_power_dbm
        power_reduction_db = input_power_dbm - target_power_dbm;
        power_reduction_linear = 10^(power_reduction_db/20);
        fprintf('Warning: Input power too high, reducing by %.1f dB\n', power_reduction_db);
        PD.x = PD.x / power_reduction_linear;
        
        % Recalculate power
        input_power_dbm = 10*log10(mean(abs(PD.x).^2)*1000);
        fprintf('Adjusted power: %.2f dBm\n', input_power_dbm);
    end
    
    % Output PA    
    [y, RMSout, Idc, Vdc] = RFWebLab_PA_meas_v1_2(PD.x, RMSin);
    
    % Check if y is a vector
    if ~isvector(y)
        fprintf('Warning: y is not a vector, dimensions: %s\n', mat2str(size(y)));
        y = y(:);  % Convert to column vector
    end
    
    % Align signals
    try
        y=safe_timealign(PD.x,y);
    catch ME
        fprintf('Signal alignment failed: %s\n', ME.message);
        fprintf('Attempting manual alignment...\n');
        
        % Manual simple alignment
        max_lag = 100;
        correlation = zeros(2*max_lag+1, 1);
        for lag = -max_lag:max_lag
            if lag > 0
                corr_segment = y(1+lag:end);
                ref_segment = PD.x(1:end-lag);
            else
                corr_segment = y(1:end+lag);
                ref_segment = PD.x(1-lag:end);
            end
            if length(corr_segment) > 100
                correlation(lag+max_lag+1) = abs(corr_segment' * ref_segment);
            end
        end
        
        [~, max_idx] = max(correlation);
        optimal_lag = max_idx - max_lag - 1;
        
        if optimal_lag > 0
            y = y(1+optimal_lag:end);
            PD.x = PD.x(1:end-optimal_lag);
        else
            PD.x = PD.x(1-optimal_lag:end);
            y = y(1:end+optimal_lag);
        end
        
        % Ensure length matches
        min_len = min(length(PD.x), length(y));
        PD.x = PD.x(1:min_len);
        y = y(1:min_len);
        
        fprintf('Manual alignment completed, lag: %d samples\n', optimal_lag);
    end

    SystIter(n).y=y;
    [SystIter(n).ACPR, PSDy]=acpr(y,Fs,ACPR);
    SystIter(n).ACPRimpr=acprdiff(ACPRout,SystIter(n).ACPR);
    Y=PSDy.Data;
    SystIter(n).PY=10*log10(Y);
    
    hold on; plot(PSDy.Frequencies*1e-6,10*log10(Y))
    
    y=fixpwr(y,RMSin);
    
    Data.In=y; 
    Data.Out=PD.x;  
    
    % Add numerical stability check
    if n > 1
        % Check data range
        if max(abs(y)) > 10 || max(abs(PD.x)) > 10
            fprintf('Warning: Data range too large, applying normalization\n');
            scale_factor = max([max(abs(y)), max(abs(PD.x))]);
            y = y / scale_factor;
            PD.x = PD.x / scale_factor;
        end
    end
    
        % Preprocess data to improve numerical stability
    try
        [Data.In, Data.Out] = preprocess_dpd_data(Data.In, Data.Out);
    catch ME
        fprintf('Data preprocessing failed: %s\n', ME.message);
        % Continue with original data
    end
    
    % Add stronger regularization options
    Options.Evaluation.Enable='on';
    Options.Evaluation.Regularization = 1e-3;  % Increase regularization
    Options.Evaluation.MaxCondition = 1e3;     % Further limit condition number
    Options.Evaluation.DampingFactor = 0.7;    % Damping factor
    
    % Use modified options for identification
    try
        [PD.Model,Eval]=gmpidentifier(Data,PD.Model,Options);
    catch ME
        fprintf('DPD identification failed: %s\n', ME.message);
        % Continue with identity model
        PD.Model.Coeff = zeros(size(PD.Model.Coeff));
        PD.Model.Coeff(1) = 1;
        Eval = struct();
    end

    % Use modified options for identification
    [PD.Model,Eval]=gmpidentifier(Data,PD.Model,Options);
    
    %% Test DPD effect    
    SystIter(n).Idc=Idc;
    SystIter(n).Vdc=Vdc;
    SystIter(n).RMSout=RMSout;
    SystIter(n).PD=PD;
    SystIter(n).Eval=Eval;
    SystIter(n).StageIterNb=n;
    SystIter(n).StageNb=1;
    
    % Display current iteration results
    if isfield(SystIter(n).ACPR, 'L1')
        % Use L1/U1/L2/U2 fields
        fprintf('Iteration %d completed: ACPR1 = %.2f dBc, ACPR2 = %.2f dBc\n', ...
            n, mean([SystIter(n).ACPR.L1, SystIter(n).ACPR.U1]), ...
            mean([SystIter(n).ACPR.L2, SystIter(n).ACPR.U2]));
    elseif isfield(SystIter(n).ACPR, 'L') && isfield(SystIter(n).ACPR, 'R')
        % Use L/R fields
        fprintf('Iteration %d completed: ACPR1 = %.2f dBc, ACPR2 = %.2f dBc\n', ...
            n, mean([SystIter(n).ACPR.L(1), SystIter(n).ACPR.R(1)]), ...
            mean([SystIter(n).ACPR.L(2), SystIter(n).ACPR.R(2)]));
    else
        fprintf('Iteration %d completed: ACPR structure format unknown\n', n);
    end
end