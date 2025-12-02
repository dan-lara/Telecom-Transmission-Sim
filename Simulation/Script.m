clc;
clear;
close all;

%% ==========================================
%% 1. 全局参数设置
%% ==========================================
N = 50000;              % 比特数 (用于BER计算)
N_plot = 2000;          % 绘图用的点数 (避免星座图太密集)
EbNo_dB_range = 0:2:12; % Eb/N0 扫描范围
mod_types = {'BPSK', 'QPSK', '16QAM'}; % 要测试的调制列表

% 滤波器与载波参数 (已修正为 fc=5e3)
fc = 5e3; Rb = 1e3; sps = 32; fs = Rb*sps;
rolloff = 0.35; span = 10;
h_rrc = rcosdesign(rolloff, span, sps, "sqrt");
total_delay = span * sps;

% 准备绘图
markers = {'-ob', '-sr', '-^g'}; 
fig_ber = figure('Name', 'Task 2: BER Curves'); hold on; grid on;

%% ==========================================
%% 2. 主循环：遍历调制方式 (BPSK -> QPSK -> 16QAM)
%% ==========================================
for m_idx = 1:length(mod_types)
    MODULATION_TYPE = mod_types{m_idx};
    
    % --- 参数配置 ---
    switch MODULATION_TYPE
        case 'BPSK', M = 2; k = 1; phase_offset = 0;
        case 'QPSK', M = 4; k = 2; phase_offset = pi/4;
        case '16QAM', M = 16; k = 4; phase_offset = 0;
    end
    
    % 确保 N 是 k 的整数倍
    N_aligned = floor(N / k) * k;
    
    disp(['Processing: ', MODULATION_TYPE, ' ...']);
    
    % --- [A] 信号生成 (Transmitter) ---
    tx_bits = randi([0 1], N_aligned, 1);
    
    if (M == 16)
        tx_symbols = qammod(tx_bits, M, 'InputType', 'bit', 'UnitAveragePower', true);
    else
        tx_symbols = pskmod(tx_bits, M, phase_offset, 'InputType', 'bit');
    end
    
    tx_filtered = upfirdn(tx_symbols, h_rrc, sps);
    
    % 上变频 (Mixer Up)
    len_sig = length(tx_filtered);
    t = (0 : len_sig - 1)' / fs; 
    I_sig = real(tx_filtered);
    Q_sig = imag(tx_filtered);
    tx_rf = I_sig .* cos(2*pi*fc*t) - Q_sig .* sin(2*pi*fc*t);
    
    %% ==========================================
    %% 3. 可视化部分 (新增功能)
    %% ==========================================
    % 仅使用前 N_plot 个符号进行绘图，避免卡顿
    num_plot_symbols = min(N_plot, length(tx_symbols));
    
    % --- 图1: 发射前理想星座图 ---
    h1 = scatterplot(tx_symbols(1:num_plot_symbols));
    title(['[Tx] Before Transmission: ', MODULATION_TYPE]);
    grid on;
    
    % --- 图2: 无噪声接收星座图 ---
    % 模拟无噪声信道
    rx_rf_nonoise = tx_rf; 
    
    % 接收处理链 (Rx Chain)
    rx_base_no = rx_rf_nonoise .* (cos(2*pi*fc*t) - 1j*sin(2*pi*fc*t)) * 2;
    rx_filt_no = upfirdn(rx_base_no, h_rrc, 1, 1);
    rx_sym_no  = rx_filt_no(total_delay + 1 : sps : end);
    rx_sym_no  = rx_sym_no(1:length(tx_symbols)); % 截断
    
    h2 = scatterplot(rx_sym_no(1:num_plot_symbols));
    title(['[Rx] After Transmission (No Noise): ', MODULATION_TYPE]);
    grid on;
    
    % --- 图3: 有噪声接收星座图 ---
    % 设定一个展示用的 SNR (例如 Eb/N0 = 10dB)
    EbNo_disp = 10; 
    snr_disp = EbNo_disp + 10*log10(k) - 10*log10(sps);
    rx_rf_noise = awgn(tx_rf, snr_disp, 'measured');
    
    % 接收处理链
    rx_base_n = rx_rf_noise .* (cos(2*pi*fc*t) - 1j*sin(2*pi*fc*t)) * 2;
    rx_filt_n = upfirdn(rx_base_n, h_rrc, 1, 1);
    rx_sym_n  = rx_filt_n(total_delay + 1 : sps : end);
    rx_sym_n  = rx_sym_n(1:length(tx_symbols));
    
    h3 = scatterplot(rx_sym_n(1:num_plot_symbols));
    title(['[Rx] With Noise (Eb/No=', num2str(EbNo_disp), 'dB): ', MODULATION_TYPE]);
    grid on;
    
    %% ==========================================
    %% 4. Task 2: BER 循环计算
    %% ==========================================
    ber_results = zeros(size(EbNo_dB_range));
    
    for i = 1:length(EbNo_dB_range)
        EbNo = EbNo_dB_range(i);
        
        % 计算当前 Waveform SNR
        snr_val = EbNo + 10*log10(k) - 10*log10(sps);
        
        % 加噪声
        rx_rf_ber = awgn(tx_rf, snr_val, 'measured');
        
        % 接收解调
        rx_base_ber = rx_rf_ber .* (cos(2*pi*fc*t) - 1j*sin(2*pi*fc*t)) * 2;
        rx_filt_ber = upfirdn(rx_base_ber, h_rrc, 1, 1);
        rx_sym_ber  = rx_filt_ber(total_delay + 1 : sps : end);
        rx_sym_ber  = rx_sym_ber(1:length(tx_symbols));
        
        % 判决
        if (M == 16)
            rx_bits_out = qamdemod(rx_sym_ber, M, 'OutputType', 'bit', 'UnitAveragePower', true);
        else
            rx_bits_out = pskdemod(rx_sym_ber, M, phase_offset, 'OutputType', 'bit');
        end
        
        % 误码率
        [~, ber_results(i)] = biterr(tx_bits, rx_bits_out);
    end
    
    % 绘制到总图上
    figure(fig_ber); % 切换回 BER 图窗口
    semilogy(EbNo_dB_range, ber_results, markers{m_idx}, 'LineWidth', 1.5, 'DisplayName', MODULATION_TYPE);
end

%% ==========================================
%% 5. 最终图表美化
%% ==========================================
figure(fig_ber);
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
title('Task 2: BER Performance Comparison');
legend('show');
ylim([1e-5 1]); 
grid on;
set(gca, 'YScale', 'log'); % 强制对数坐标
disp('All simulations completed!');