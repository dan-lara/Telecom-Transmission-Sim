clc;
clear;
close all;

%% 1. 全局参数设置
N = 50000;              % 比特数 (调试时可用5e4，正式跑图建议1e5或1e6)
EbNo_dB_range = 0:2:12; % Eb/N0 扫描范围 (例如 0, 2, ..., 12 dB)
mod_types = {'BPSK', 'QPSK', '16QAM'}; % 要测试的调制列表

% 滤波器与载波参数 (保持 Task 1 配置)
fc = 1e6; Rb = 1e3; sps = 32; fs = Rb*sps;
rolloff = 0.35; span = 10;
h_rrc = rcosdesign(rolloff, span, sps, "sqrt");
total_delay = span * sps;

% 准备颜色和标记，用于画图
markers = {'-ob', '-sr', '-^g'}; 
figure; hold on; grid on;

%% 2. 外层循环：遍历不同的调制方式 (BPSK -> QPSK -> 16QAM)
for m_idx = 1:length(mod_types)
    MODULATION_TYPE = mod_types{m_idx};
    
    % 根据调制类型设定参数
    switch MODULATION_TYPE
        case 'BPSK'
            M = 2; k = 1; phase_offset = 0;
        case 'QPSK'
            M = 4; k = 2; phase_offset = pi/4;
        case '16QAM'
            M = 16; k = 4; phase_offset = 0;
    end
    
    % 确保 N 是 k 的整数倍
    N_aligned = floor(N / k) * k;
    
    % 预分配 BER 结果数组
    ber_results = zeros(size(EbNo_dB_range));
    
    disp(['Running simulation for: ', MODULATION_TYPE, ' ...']);
    
    %% 3. 内层循环：遍历信噪比 (Eb/N0)
    for i = 1:length(EbNo_dB_range)
        EbNo = EbNo_dB_range(i);
        
        % --- [A] Transmitter ---
        tx_bits = randi([0 1], N_aligned, 1);
        
        if (M == 16)
            tx_symbols = qammod(tx_bits, M, 'InputType', 'bit', 'UnitAveragePower', true);
        else
            tx_symbols = pskmod(tx_bits, M, phase_offset, 'InputType', 'bit');
        end
        
        tx_filtered = upfirdn(tx_symbols, h_rrc, sps);
        
        % Up-conversion (Mixer)
        len_sig = length(tx_filtered);
        t = (0 : len_sig - 1)' / fs; % 列向量
        I_sig = real(tx_filtered);
        Q_sig = imag(tx_filtered);
        tx_rf = I_sig .* cos(2*pi*fc*t) - Q_sig .* sin(2*pi*fc*t);
        
        % --- [B] Channel (添加噪声) ---
        % 关键公式：Waveform Level SNR 计算
        snr_val = EbNo + 10*log10(k) - 10*log10(sps);
        
        % 添加高斯白噪声
        rx_rf = awgn(tx_rf, snr_val, 'measured');
        
        % --- [C] Receiver ---
        % Down-conversion (Mixer) * 2 补偿损耗
        rx_baseband = rx_rf .* (cos(2*pi*fc*t) - 1j*sin(2*pi*fc*t)) * 2;
        
        % Matched Filter
        rx_filtered = upfirdn(rx_baseband, h_rrc, 1, 1);
        
        % Downsampling
        rx_symbols = rx_filtered(total_delay + 1 : sps : end);
        
        % Truncation
        num_tx_symbols = length(tx_bits) / k;
        rx_symbols = rx_symbols(1:num_tx_symbols);
        
        % Demodulation
        if (M == 16)
            rx_bits = qamdemod(rx_symbols, M, 'OutputType', 'bit', 'UnitAveragePower', true);
        else
            rx_bits = pskdemod(rx_symbols, M, phase_offset, 'OutputType', 'bit');
        end
        
        % --- [D] BER Calc ---
        [~, ber] = biterr(tx_bits, rx_bits);
        ber_results(i) = ber;
    end
    
    %% 4. 绘制该调制的曲线
    semilogy(EbNo_dB_range, ber_results, markers{m_idx}, 'LineWidth', 1.5, 'DisplayName', MODULATION_TYPE);
end

%% 5. 图表美化
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
title('BER Performance Comparison (Waveform Level Simulation)');
legend('show');
ylim([1e-5 1]); % 限制 Y 轴范围
grid on;
set(gca, 'YScale', 'log');
disp('Task II Simulation Finished!');