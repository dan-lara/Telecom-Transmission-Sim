clc;
clear;
close all;

%% ==========================================
%% 1. Configuration des paramètres globaux
%% ==========================================
% Définir tous les types de modulation à tester
MODULATION_TYPES = {'BPSK', 'QPSK', '16QAM'};

% Paramètres de base
fc = 5e3;       % Fréquence porteuse
Rb = 1e3;       % Débit symbole
sps = 32;       % Facteur de suréchantillonnage
fs = Rb*sps;    % Fréquence d'échantillonnage
N = 10000;      % Nombre de bits initial

% Configuration du filtre (commune à toutes les modulations)
rolloff = 0.35; span = 10;
h_rrc = rcosdesign(rolloff, span, sps, "sqrt");
total_delay = span * sps;

% SNR pour test fixe
snr_target_fixed = 20; % dB

%% ==========================================
%% 2. Test avec SNR fixe - Vérification sans bruit et avec bruit
%% ==========================================
for mod_idx = 1:length(MODULATION_TYPES)
    MODULATION_TYPE = MODULATION_TYPES{mod_idx};
    disp('==========================================');
    disp(['Test de la modulation : ', MODULATION_TYPE]);
    disp('==========================================');
    
    % Configuration des paramètres de modulation
    switch MODULATION_TYPE
        case 'BPSK'
            M = 2; k = 1; phase_offset = 0;
        case 'QPSK'
            M = 4; k = 2; phase_offset = pi/4;
        case '16QAM'
            M = 16; k = 4; phase_offset = 0;
        otherwise
            error('Type de modulation inconnu');
    end

    % Assurer l'alignement du nombre de bits
    N_adjusted = floor(N / k) * k; 
    tx_bits = randi([0 1], N_adjusted, 1);

    %% ==========================================
    %% 2.1 Émetteur (Transmitter) - Deux modes indépendants
    %% ==========================================
    
    % --- Mode A : Automatique (MATLAB Toolbox) ---
    if M == 16
        tx_sym_auto = qammod(tx_bits, M, 'InputType', 'bit', 'UnitAveragePower', true);
    else
        tx_sym_auto = pskmod(tx_bits, M, phase_offset, 'InputType', 'bit');
    end

    % --- Mode B : Manuel (Calcul mathématique) ---
    tx_sym_manual = zeros(length(tx_bits)/k, 1);

    if M == 2 % BPSK manuel
        % Mapping manuel cohérent
        tx_sym_manual = 1 - 2*tx_bits;  % 0->1, 1->-1
        
    elseif M == 4 % QPSK manuel
        bits_r = reshape(tx_bits, 2, []);
        dec_idx = bits_r(1,:)*2 + bits_r(2,:);
        
        % Table de correspondance du code de Gray
        gray_map = [0, 1, 3, 2];
        sym_idx = gray_map(dec_idx + 1);
        
        % Générer les symboles
        tx_sym_manual = exp(1j * (phase_offset + sym_idx * 2*pi/4)).';
        
    elseif M == 16 % 16-QAM manuel
        % NOTE: Ce mapping est différent de celui de MATLAB qammod
        % C'est intentionnel pour démontrer la diversité des implémentations
        bits_r = reshape(tx_bits, 4, []);
        
        % Mapping personnalisé
        map_level = [-3, -1, 3, 1]; 
        
        idx_I = bits_r(1,:)*2 + bits_r(2,:);
        idx_Q = bits_r(3,:)*2 + bits_r(4,:);
        
        val_I = map_level(idx_I + 1);
        val_Q = map_level(idx_Q + 1);
        
        % Normalisation de la puissance
        tx_sym_manual = (val_I + 1j * val_Q).' / sqrt(10);
    end

    %% ==========================================
    %% 2.2 Test sans bruit (Vérification de la chaîne)
    %% ==========================================
    disp('--- Test sans bruit ---');
    
    % Utiliser la chaîne automatique pour la vérification
    tx_symbols_test = tx_sym_auto;
    
    % Filtrage RRC
    tx_filtered_test = upfirdn(tx_symbols_test, h_rrc, sps);
    
    % Modulation
    t_test = (0 : length(tx_filtered_test) - 1)' / fs;
    tx_rf_test = real(tx_filtered_test).*cos(2*pi*fc*t_test) - imag(tx_filtered_test).*sin(2*pi*fc*t_test);
    
    % Canal sans bruit
    rx_rf_test = tx_rf_test;  % Pas de bruit
    
    % Réception
    rx_base_test = rx_rf_test .* (cos(2*pi*fc*t_test) - 1j*sin(2*pi*fc*t_test)) * 2;
    rx_filt_test = upfirdn(rx_base_test, h_rrc, 1, 1);
    rx_sym_test = rx_filt_test(total_delay + 1 : sps : end);
    rx_sym_test = rx_sym_test(1:length(tx_symbols_test));
    
    % Démodulation automatique
    if M == 16
        rx_bits_test = qamdemod(rx_sym_test, M, 'OutputType', 'bit', 'UnitAveragePower', true);
    else
        rx_bits_test = pskdemod(rx_sym_test, M, phase_offset, 'OutputType', 'bit');
    end
    
    % Calcul du TEB sans bruit
    [~, ber_no_noise] = biterr(tx_bits, rx_bits_test);
    disp(['TEB sans bruit : ', num2str(ber_no_noise)]);

    %% ==========================================
    %% 2.3 Test avec bruit (SNR fixe)
    %% ==========================================
    disp('--- Test avec bruit (SNR fixe) ---');
    
    % Deux chaînes indépendantes
    ber_results = zeros(2, 1); % [Auto; Manuel]
    
    for chain_idx = 1:2
        if chain_idx == 1
            % Chaîne automatique
            tx_symbols = tx_sym_auto;
            chain_name = 'Automatique';
        else
            % Chaîne manuelle
            tx_symbols = tx_sym_manual;
            chain_name = 'Manuelle';
        end
        
        % Filtrage RRC
        tx_filtered = upfirdn(tx_symbols, h_rrc, sps);
        
        % Modulation
        t = (0 : length(tx_filtered) - 1)' / fs;
        tx_rf = real(tx_filtered).*cos(2*pi*fc*t) - imag(tx_filtered).*sin(2*pi*fc*t);
        
        % Canal avec bruit
        snr_val = snr_target_fixed + 10*log10(k) - 10*log10(sps);
        rx_rf = awgn(tx_rf, snr_val, 'measured');
        
        % Réception
        rx_base = rx_rf .* (cos(2*pi*fc*t) - 1j*sin(2*pi*fc*t)) * 2;
        rx_filt = upfirdn(rx_base, h_rrc, 1, 1);
        rx_sym = rx_filt(total_delay + 1 : sps : end);
        rx_sym = rx_sym(1:length(tx_symbols));
        
        % Démodulation
        if chain_idx == 1
            % Démodulation automatique
            if M == 16
                rx_bits = qamdemod(rx_sym, M, 'OutputType', 'bit', 'UnitAveragePower', true);
            else
                rx_bits = pskdemod(rx_sym, M, phase_offset, 'OutputType', 'bit');
            end
        else
            % Démodulation manuelle
            rx_bits = demapping_manual(rx_sym, M, phase_offset, k);
        end
        
        % Calcul du TEB
        [~, ber] = biterr(tx_bits, rx_bits);
        ber_results(chain_idx) = ber;
        
        % Stocker les résultats pour visualisation
        if chain_idx == 1
            rx_sym_auto = rx_sym;
            ber_auto = ber;
        else
            rx_sym_manual_chain = rx_sym;
            ber_manual = ber;
        end
    end
    
    disp(['TEB Chaîne Automatique (SNR=', num2str(snr_target_fixed), 'dB) : ', num2str(ber_results(1))]);
    disp(['TEB Chaîne Manuelle (SNR=', num2str(snr_target_fixed), 'dB) : ', num2str(ber_results(2))]);

    %% ==========================================
    %% 2.4 Visualisation des résultats (SNR fixe)
    %% ==========================================
    
    % Figure 1 : Constellations émises
    figure('Name', [MODULATION_TYPE, ' - Constellations Émises'], ...
           'Position', [100, 100, 1200, 400]);
    
    n_show = min(200, length(tx_sym_auto));
    
    subplot(1, 3, 1);
    scatter(real(tx_sym_auto(1:n_show)), imag(tx_sym_auto(1:n_show)), 40, 'b', 'filled');
    grid on; axis equal; xlim([-1.5 1.5]); ylim([-1.5 1.5]);
    xlabel('I'); ylabel('Q');
    title(['Auto - Émis (', MODULATION_TYPE, ')']);
    
    subplot(1, 3, 2);
    scatter(real(tx_sym_manual(1:n_show)), imag(tx_sym_manual(1:n_show)), 40, 'r', 'filled');
    grid on; axis equal;
    xlabel('I'); ylabel('Q');
    title(['Manuel - Émis (', MODULATION_TYPE, ')']);
    
    subplot(1, 3, 3);
    plot(real(tx_sym_auto(1:50)), 'b-', 'LineWidth', 1.5); hold on;
    plot(real(tx_sym_manual(1:50)), 'r--', 'LineWidth', 1.5);
    grid on; xlabel('Index du symbole'); ylabel('Amplitude I');
    title('Comparaison des symboles émis (Partie réelle)');
    legend('Auto', 'Manuel');
    
    % Figure 2 : Signaux reçus avec bruit
    figure('Name', [MODULATION_TYPE, ' - Signaux Reçus avec Bruit'], ...
           'Position', [100, 550, 1200, 400]);
    
    subplot(1, 3, 1);
    scatter(real(rx_sym_auto(1:n_show)), imag(rx_sym_auto(1:n_show)), 30, 'b', 'filled', 'MarkerFaceAlpha', 0.6);
    grid on; axis equal;
    xlabel('I'); ylabel('Q');
    title(['Auto - Reçu (SNR=', num2str(snr_target_fixed), 'dB)']);
    
    subplot(1, 3, 2);
    scatter(real(rx_sym_manual_chain(1:n_show)), imag(rx_sym_manual_chain(1:n_show)), 30, 'r', 'filled', 'MarkerFaceAlpha', 0.6);
    grid on; axis equal;
    xlabel('I'); ylabel('Q');
    title(['Manuel - Reçu (SNR=', num2str(snr_target_fixed), 'dB)']);
    
    subplot(1, 3, 3);
    plot(real(rx_sym_auto(1:100)), 'b-', 'LineWidth', 1); hold on;
    plot(imag(rx_sym_auto(1:100)), 'b--', 'LineWidth', 1);
    plot(real(rx_sym_manual_chain(1:100)), 'r-', 'LineWidth', 1);
    plot(imag(rx_sym_manual_chain(1:100)), 'r--', 'LineWidth', 1);
    grid on; xlabel('Échantillon'); ylabel('Amplitude');
    title('Signaux reçus (temporel)');
    legend('I Auto', 'Q Auto', 'I Manuel', 'Q Manuel');
    
    % Figure 3 : Densité spectrale de puissance
    figure('Name', [MODULATION_TYPE, ' - Densité Spectrale de Puissance'], ...
           'Position', [100, 1000, 1200, 400]);
    
    % Calcul de la DSP
    [psd_auto, f_auto] = pwelch(tx_rf, 1024, 512, 1024, fs);
    [psd_manual, f_manual] = pwelch(tx_rf, 1024, 512, 1024, fs);
    
    subplot(1, 2, 1);
    plot(f_auto/1000, 10*log10(psd_auto), 'b-', 'LineWidth', 1.5);
    grid on; xlabel('Fréquence (kHz)'); ylabel('Puissance (dB)');
    title(['DSP - Chaîne Auto (', MODULATION_TYPE, ')']);
    xlim([0, fs/2000]);
    
    subplot(1, 2, 2);
    plot(f_manual/1000, 10*log10(psd_manual), 'r-', 'LineWidth', 1.5);
    grid on; xlabel('Fréquence (kHz)'); ylabel('Puissance (dB)');
    title(['DSP - Chaîne Manuel (', MODULATION_TYPE, ')']);
    xlim([0, fs/2000]);
    
    % Figure 4 : Comparaison des performances
    figure('Name', [MODULATION_TYPE, ' - Comparaison des Performances'], ...
           'Position', [100, 1450, 600, 400]);
    
    bar([ber_auto, ber_manual]);
    set(gca, 'XTickLabel', {'Auto', 'Manuel'});
    ylabel('Taux d''Erreur Binaire (TEB)');
    title([MODULATION_TYPE, ' - Comparaison des TEB (SNR=', num2str(snr_target_fixed), 'dB)']);
    grid on;
    
    % Ajouter les valeurs sur les barres
    text(1, ber_auto + max([ber_auto, ber_manual])*0.1, sprintf('%.4f', ber_auto), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    text(2, ber_manual + max([ber_auto, ber_manual])*0.1, sprintf('%.4f', ber_manual), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    
    disp(['Test de modulation ', MODULATION_TYPE, ' terminé']);
    disp(' ');
end

%% ==========================================
%% 3. Analyse comparative : BER vs SNR pour différentes modulations
%% ==========================================
disp('==========================================');
disp('Analyse comparative : BER vs SNR');
disp('==========================================');

% Paramètres pour l'analyse BER vs SNR
SNR_dB_range = 0:2:20;  % Plage de SNR à tester
num_snr_points = length(SNR_dB_range);
num_modulations = length(MODULATION_TYPES);

% Initialisation des matrices de résultats
BER_auto_matrix = zeros(num_modulations, num_snr_points);
BER_manual_matrix = zeros(num_modulations, num_snr_points);

% Paramètres de simulation pour l'analyse BER
N_ber = 20000;  % Nombre de bits pour l'analyse BER

for mod_idx = 1:num_modulations
    MODULATION_TYPE = MODULATION_TYPES{mod_idx};
    disp(['Analyse BER vs SNR pour : ', MODULATION_TYPE]);
    
    % Configuration des paramètres de modulation
    switch MODULATION_TYPE
        case 'BPSK'
            M = 2; k = 1; phase_offset = 0;
        case 'QPSK'
            M = 4; k = 2; phase_offset = pi/4;
        case '16QAM'
            M = 16; k = 4; phase_offset = 0;
    end
    
    % Assurer l'alignement du nombre de bits
    N_adjusted_ber = floor(N_ber / k) * k;
    
    for snr_idx = 1:num_snr_points
        snr_target = SNR_dB_range(snr_idx);
        
        % Génération des bits aléatoires
        tx_bits_ber = randi([0 1], N_adjusted_ber, 1);
        
        % Chaîne automatique
        if M == 16
            tx_sym_auto_ber = qammod(tx_bits_ber, M, 'InputType', 'bit', 'UnitAveragePower', true);
        else
            tx_sym_auto_ber = pskmod(tx_bits_ber, M, phase_offset, 'InputType', 'bit');
        end
        
        % Transmission et réception automatique
        [ber_auto, ~] = simulate_chain(tx_sym_auto_ber, tx_bits_ber, M, phase_offset, ...
                                       snr_target, k, sps, fc, fs, h_rrc, total_delay, 'auto');
        BER_auto_matrix(mod_idx, snr_idx) = ber_auto;
        
        % Chaîne manuelle
        if M == 2
            tx_sym_manual_ber = 1 - 2*tx_bits_ber;
        elseif M == 4
            bits_r = reshape(tx_bits_ber, 2, []);
            dec_idx = bits_r(1,:)*2 + bits_r(2,:);
            gray_map = [0, 1, 3, 2];
            sym_idx = gray_map(dec_idx + 1);
            tx_sym_manual_ber = exp(1j * (phase_offset + sym_idx * 2*pi/4)).';
        elseif M == 16
            bits_r = reshape(tx_bits_ber, 4, []);
            map_level = [-3, -1, 3, 1];
            idx_I = bits_r(1,:)*2 + bits_r(2,:);
            idx_Q = bits_r(3,:)*2 + bits_r(4,:);
            val_I = map_level(idx_I + 1);
            val_Q = map_level(idx_Q + 1);
            tx_sym_manual_ber = (val_I + 1j * val_Q).' / sqrt(10);
        end
        
        % Transmission et réception manuelle
        [ber_manual, ~] = simulate_chain(tx_sym_manual_ber, tx_bits_ber, M, phase_offset, ...
                                        snr_target, k, sps, fc, fs, h_rrc, total_delay, 'manual');
        BER_manual_matrix(mod_idx, snr_idx) = ber_manual;
    end
    
    disp(['  Terminé pour ', MODULATION_TYPE]);
end

%% ==========================================
%% 4. Visualisation des résultats BER vs SNR
%% ==========================================

% Figure 5 : BER vs SNR pour toutes les modulations (Chaîne automatique)
figure('Name', 'BER vs SNR - Chaîne Automatique', ...
       'Position', [800, 100, 800, 600]);

colors = {'b-o', 'r-s', 'g-^'};
line_width = 1.5;
marker_size = 8;

for mod_idx = 1:num_modulations
    semilogy(SNR_dB_range, BER_auto_matrix(mod_idx, :), colors{mod_idx}, ...
             'LineWidth', line_width, 'MarkerSize', marker_size);
    hold on;
end

grid on;
xlabel('SNR (dB)', 'FontSize', 12);
ylabel('Bit Error Rate (BER)', 'FontSize', 12);
title('Performance des différentes modulations (Chaîne Automatique)', 'FontSize', 14);
legend('BPSK', 'QPSK', '16QAM', 'Location', 'best');
ylim([1e-5, 1]);

% Ajouter les courbes théoriques (approximation)
hold on;
SNR_lin = 10.^(SNR_dB_range/10);

% Courbe théorique BPSK
ber_bpsk_theo = 0.5 * erfc(sqrt(SNR_lin));
semilogy(SNR_dB_range, ber_bpsk_theo, 'b--', 'LineWidth', 1);

% Courbe théorique QPSK (identique à BPSK en terme de Eb/N0)
ber_qpsk_theo = 0.5 * erfc(sqrt(SNR_lin));
semilogy(SNR_dB_range, ber_qpsk_theo, 'r--', 'LineWidth', 1);

% Courbe théorique 16QAM (approximation)
ber_16qam_theo = (3/8) * erfc(sqrt(0.4 * SNR_lin));
semilogy(SNR_dB_range, ber_16qam_theo, 'g--', 'LineWidth', 1);

legend('BPSK Sim', 'QPSK Sim', '16QAM Sim', 'BPSK Theo', 'QPSK Theo', '16QAM Theo', 'Location', 'best');

% Figure 6 : Comparaison Auto vs Manuel pour chaque modulation
figure('Name', 'Comparaison Auto vs Manuel - BER vs SNR', ...
       'Position', [800, 750, 1200, 800]);

for mod_idx = 1:num_modulations
    subplot(2, 2, mod_idx);
    MODULATION_TYPE = MODULATION_TYPES{mod_idx};
    
    semilogy(SNR_dB_range, BER_auto_matrix(mod_idx, :), 'b-o', ...
             'LineWidth', line_width, 'MarkerSize', marker_size);
    hold on;
    semilogy(SNR_dB_range, BER_manual_matrix(mod_idx, :), 'r-s', ...
             'LineWidth', line_width, 'MarkerSize', marker_size);
    
    grid on;
    xlabel('SNR (dB)', 'FontSize', 10);
    ylabel('BER', 'FontSize', 10);
    title(['Comparaison Auto vs Manuel - ', MODULATION_TYPE], 'FontSize', 12);
    legend('Auto', 'Manuel', 'Location', 'best');
    ylim([1e-5, 1]);
end

% Sous-graphique supplémentaire : Efficacité spectrale
subplot(2, 2, 4);
spectral_efficiency = [1, 2, 4];  % BPSK:1, QPSK:2, 16QAM:4
bar(spectral_efficiency);
set(gca, 'XTickLabel', {'BPSK', 'QPSK', '16QAM'});
xlabel('Modulation');
ylabel('Efficacité spectrale (bits/s/Hz)');
title('Efficacité spectrale des différentes modulations');
grid on;

% Ajouter les valeurs sur les barres
for i = 1:length(spectral_efficiency)
    text(i, spectral_efficiency(i) + 0.1, num2str(spectral_efficiency(i)), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end

%% ==========================================
%% 5. Analyse des résultats
%% ==========================================
disp('==========================================');
disp('Résumé des performances :');
disp('==========================================');

for mod_idx = 1:num_modulations
    MODULATION_TYPE = MODULATION_TYPES{mod_idx};
    
    % SNR requis pour BER = 10^-3
    ber_target = 1e-3;
    ber_auto = BER_auto_matrix(mod_idx, :);
    ber_manual = BER_manual_matrix(mod_idx, :);
    
    % Trouver le SNR minimal pour atteindre la cible
    idx_auto = find(ber_auto <= ber_target, 1);
    idx_manual = find(ber_manual <= ber_target, 1);
    
    if ~isempty(idx_auto)
        snr_auto = SNR_dB_range(idx_auto);
    else
        snr_auto = NaN;
    end
    
    if ~isempty(idx_manual)
        snr_manual = SNR_dB_range(idx_manual);
    else
        snr_manual = NaN;
    end
    
    disp(['Modulation : ', MODULATION_TYPE]);
    disp(['  SNR requis (Auto) pour BER=10^-3 : ', num2str(snr_auto), ' dB']);
    disp(['  SNR requis (Manuel) pour BER=10^-3 : ', num2str(snr_manual), ' dB']);
    
    % Calcul du gain de codage (différence entre auto et manuel)
    if ~isnan(snr_auto) && ~isnan(snr_manual)
        coding_gain = snr_manual - snr_auto;
        disp(['  Gain de codage : ', num2str(coding_gain), ' dB']);
    end
    disp(' ');
end

%% ==========================================
%% 6. Fonctions auxiliaires
%% ==========================================

% Fonction de simulation de chaîne
function [ber, rx_sym] = simulate_chain(tx_symbols, tx_bits, M, phase_offset, ...
                                        snr_target, k, sps, fc, fs, h_rrc, total_delay, chain_type)
    
    % Filtrage RRC
    tx_filtered = upfirdn(tx_symbols, h_rrc, sps);
    
    % Modulation
    t = (0 : length(tx_filtered) - 1)' / fs;
    tx_rf = real(tx_filtered).*cos(2*pi*fc*t) - imag(tx_filtered).*sin(2*pi*fc*t);
    
    % Canal avec bruit
    snr_val = snr_target + 10*log10(k) - 10*log10(sps);
    rx_rf = awgn(tx_rf, snr_val, 'measured');
    
    % Réception
    rx_base = rx_rf .* (cos(2*pi*fc*t) - 1j*sin(2*pi*fc*t)) * 2;
    rx_filt = upfirdn(rx_base, h_rrc, 1, 1);
    rx_sym = rx_filt(total_delay + 1 : sps : end);
    rx_sym = rx_sym(1:length(tx_symbols));
    
    % Démodulation
    if strcmp(chain_type, 'auto')
        % Démodulation automatique
        if M == 16
            rx_bits = qamdemod(rx_sym, M, 'OutputType', 'bit', 'UnitAveragePower', true);
        else
            rx_bits = pskdemod(rx_sym, M, phase_offset, 'OutputType', 'bit');
        end
    else
        % Démodulation manuelle
        rx_bits = demapping_manual(rx_sym, M, phase_offset, k);
    end
    
    % Calcul du BER
    [~, ber] = biterr(tx_bits, rx_bits);
end

% Fonction de démodulation manuelle
function rx_bits = demapping_manual(rx_sym, M, phase_offset, k)
    num_symbols = length(rx_sym);
    rx_bits = zeros(num_symbols * k, 1);
    
    if M == 2 % BPSK
        % Décision basée sur la partie réelle
        decisions = real(rx_sym) < 0;
        rx_bits = double(decisions);
        
    elseif M == 4 % QPSK
        % Décision par quadrant (méthode corrigée)
        for i = 1:num_symbols
            real_part = real(rx_sym(i));
            imag_part = imag(rx_sym(i));
            start_idx = (i-1)*2 + 1;
            
            if real_part >= 0 && imag_part >= 0
                rx_bits(start_idx:start_idx+1) = [0; 0];
            elseif real_part < 0 && imag_part >= 0
                rx_bits(start_idx:start_idx+1) = [0; 1];
            elseif real_part < 0 && imag_part < 0
                rx_bits(start_idx:start_idx+1) = [1; 1];
            else % real_part >= 0 && imag_part < 0
                rx_bits(start_idx:start_idx+1) = [1; 0];
            end
        end
        
    elseif M == 16 % 16-QAM
        % Restauration de l'amplitude
        rx_scaled = rx_sym * sqrt(10);
        I_rec = real(rx_scaled);
        Q_rec = imag(rx_scaled);
        
        % Table de correspondance inverse
        bit_patterns = [
            0, 0;  % 0 -> 00 (correspond à -3)
            0, 1;  % 1 -> 01 (correspond à -1)
            1, 0;  % 2 -> 10 (correspond à +3)
            1, 1   % 3 -> 11 (correspond à +1)
        ];
        
        for i = 1:num_symbols
            bits_I_dec = demapping_logic_16qam(I_rec(i));
            bits_Q_dec = demapping_logic_16qam(Q_rec(i));
            
            b_I = bit_patterns(bits_I_dec + 1, :);
            b_Q = bit_patterns(bits_Q_dec + 1, :);
            
            start_idx = (i-1)*4 + 1;
            rx_bits(start_idx:start_idx+1) = b_I';
            rx_bits(start_idx+2:start_idx+3) = b_Q';
        end
    end
end

% Fonction de décision pour 16-QAM
function bits_dec = demapping_logic_16qam(level_val)
    if level_val < -2
        bits_dec = 0;      % -3 -> 00
    elseif level_val < 0
        bits_dec = 1;      % -1 -> 01
    elseif level_val < 2
        bits_dec = 3;      % +1 -> 11
    else
        bits_dec = 2;      % +3 -> 10
    end
end

disp('==========================================');
disp('Tous les tests sont terminés avec succès !');
disp('==========================================');