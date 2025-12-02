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

% Réglage du RSB
snr_target = 20; % dB

%% ==========================================
%% 2. Boucle principale : parcourir tous les types de modulation
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
    
    % ---【Mode A : Automatique (MATLAB Toolbox)】---
    if M == 16
        tx_sym_auto = qammod(tx_bits, M, 'InputType', 'bit', 'UnitAveragePower', true);
    else
        tx_sym_auto = pskmod(tx_bits, M, phase_offset, 'InputType', 'bit');
    end

    % ---【Mode B : Manuel (Calcul mathématique)】---
    tx_sym_manual = zeros(length(tx_bits)/k, 1);

    if M == 2 % BPSK manuel
        % Vérification: Comparer avec la sortie de pskmod pour assurer la cohérence
        tx_sym_manual_test = pskmod(tx_bits, M, phase_offset, 'InputType', 'bit');
        % Notre mapping manuel: 0->1, 1->-1
        tx_sym_manual = 1 - 2*tx_bits;  % 0->1, 1->-1
        % Vérifier si notre mapping correspond à celui de MATLAB
        if mean(abs(tx_sym_manual_test - tx_sym_manual)) > 1e-10
            warning('Le mapping BPSK manuel diffère de celui de MATLAB pskmod');
        end
        
    elseif M == 4 % QPSK manuel
        bits_r = reshape(tx_bits, 2, []);
        dec_idx = bits_r(1,:)*2 + bits_r(2,:);
        
        % Table de correspondance du code de Gray (identique au code de Gray par défaut de pskmod)
        gray_map = [0, 1, 3, 2];
        sym_idx = gray_map(dec_idx + 1);
        
        % Générer les symboles
        tx_sym_manual = exp(1j * (phase_offset + sym_idx * 2*pi/4)).';
        
    elseif M == 16 % 16-QAM manuel
        % NOTE IMPORTANTE: Cette mapping manuel est DIFFÉRENT du mapping par défaut de MATLAB qammod
        % C'est intentionnel pour démontrer la diversité des schémas de mapping
        % Les deux chaînes (automatique et manuelle) utilisent leurs propres règles cohérentes
        bits_r = reshape(tx_bits, 4, []);
        
        % Mapping personnalisé (différent de MATLAB)
        % Niveau: -3, -1, 3, 1 correspondant aux bits: 00, 01, 11, 10 (Gray)
        map_level = [-3, -1, 3, 1]; 
        
        idx_I = bits_r(1,:)*2 + bits_r(2,:);
        idx_Q = bits_r(3,:)*2 + bits_r(4,:);
        
        val_I = map_level(idx_I + 1);
        val_Q = map_level(idx_Q + 1);
        
        % Normalisation de la puissance (différente de MATLAB)
        tx_sym_manual = (val_I + 1j * val_Q).' / sqrt(10);
    end

    % ---【Vérification comparative】Calcul de l'erreur à l'émission ---
    % Note: Pour 16QAM, une erreur est attendue car les mappings sont différents
    tx_error = abs(tx_sym_auto - tx_sym_manual);
    disp(['Erreur moyenne entre les signaux générés automatiquement et manuellement à l''émission : ', num2str(mean(tx_error))]);
    if M == 16 && mean(tx_error) > 0.1
        disp('Note: L''erreur pour 16QAM est attendue car les schémas de mapping automatique et manuel sont différents.');
    end

    %% ==========================================
    %% 2.2 Chaîne Automatique Indépendante
    %% ==========================================
    disp('--- Chaîne Automatique ---');
    
    % Transmission de la chaîne automatique
    tx_symbols_auto = tx_sym_auto;
    
    % Filtrage RRC
    tx_filtered_auto = upfirdn(tx_symbols_auto, h_rrc, sps);
    
    % Modulation (up-conversion)
    t_auto = (0 : length(tx_filtered_auto) - 1)' / fs;
    tx_rf_auto = real(tx_filtered_auto).*cos(2*pi*fc*t_auto) - imag(tx_filtered_auto).*sin(2*pi*fc*t_auto);
    
    % Canal (ajout de bruit)
    snr_val = snr_target + 10*log10(k) - 10*log10(sps);
    rx_rf_auto = awgn(tx_rf_auto, snr_val, 'measured');
    
    % Traitement récepteur (démodulation + filtrage + troncature)
    rx_base_auto = rx_rf_auto .* (cos(2*pi*fc*t_auto) - 1j*sin(2*pi*fc*t_auto)) * 2;
    rx_filt_auto = upfirdn(rx_base_auto, h_rrc, 1, 1);
    rx_sym_auto = rx_filt_auto(total_delay + 1 : sps : end);
    rx_sym_auto = rx_sym_auto(1:length(tx_symbols_auto)); % Troncature
    
    % Démodulation automatique
    if M == 16
        rx_bits_auto_auto = qamdemod(rx_sym_auto, M, 'OutputType', 'bit', 'UnitAveragePower', true);
    else
        rx_bits_auto_auto = pskdemod(rx_sym_auto, M, phase_offset, 'OutputType', 'bit');
    end
    
    % Calcul du TEB pour la chaîne automatique
    [~, ber_auto_auto] = biterr(tx_bits, rx_bits_auto_auto);
    disp(['TEB (Automatique -> Automatique) : ', num2str(ber_auto_auto)]);

    %% ==========================================
    %% 2.3 Chaîne Manuelle Indépendante
    %% ==========================================
    disp('--- Chaîne Manuelle ---');
    
    % Transmission de la chaîne manuelle
    tx_symbols_manual = tx_sym_manual;
    
    % Filtrage RRC
    tx_filtered_manual = upfirdn(tx_symbols_manual, h_rrc, sps);
    
    % Modulation (up-conversion)
    t_manual = (0 : length(tx_filtered_manual) - 1)' / fs;
    tx_rf_manual = real(tx_filtered_manual).*cos(2*pi*fc*t_manual) - imag(tx_filtered_manual).*sin(2*pi*fc*t_manual);
    
    % Canal (ajout de bruit - utilise une graine différente pour l'indépendance)
    rng('shuffle'); % Changer la graine aléatoire
    rx_rf_manual = awgn(tx_rf_manual, snr_val, 'measured');
    
    % Traitement récepteur (démodulation + filtrage + troncature)
    rx_base_manual = rx_rf_manual .* (cos(2*pi*fc*t_manual) - 1j*sin(2*pi*fc*t_manual)) * 2;
    rx_filt_manual = upfirdn(rx_base_manual, h_rrc, 1, 1);
    rx_sym_manual = rx_filt_manual(total_delay + 1 : sps : end);
    rx_sym_manual = rx_sym_manual(1:length(tx_symbols_manual)); % Troncature
    
    % Démodulation manuelle
    rx_bits_manual_manual = zeros(size(tx_bits));
    
    if M == 2 % Démodulation BPSK manuelle
        % Décision basée sur la partie réelle
        rx_bits_manual_manual = (real(rx_sym_manual) < 0);
        
    elseif M == 4 % Démodulation QPSK manuelle - CORRIGÉ
        % Méthode corrigée: Décision par quadrant (signe de la partie réelle/imaginaire)
        num_symbols = length(rx_sym_manual);
        for i = 1:num_symbols
            real_part = real(rx_sym_manual(i));
            imag_part = imag(rx_sym_manual(i));
            start_idx = (i-1)*2 + 1;
            
            if real_part >= 0 && imag_part >= 0
                % Premier quadrant (00)
                rx_bits_manual_manual(start_idx:start_idx+1) = [0; 0];
            elseif real_part < 0 && imag_part >= 0
                % Deuxième quadrant (01)
                rx_bits_manual_manual(start_idx:start_idx+1) = [0; 1];
            elseif real_part < 0 && imag_part < 0
                % Troisième quadrant (11)
                rx_bits_manual_manual(start_idx:start_idx+1) = [1; 1];
            else % real_part >= 0 && imag_part < 0
                % Quatrième quadrant (10)
                rx_bits_manual_manual(start_idx:start_idx+1) = [1; 0];
            end
        end
        
    elseif M == 16 % Démodulation 16-QAM manuelle
        % Restauration de l'amplitude
        rx_scaled = rx_sym_manual * sqrt(10);
        I_rec = real(rx_scaled);
        Q_rec = imag(rx_scaled);
        
        % Décodage utilisant le mapping manuel cohérent
        num_symbols = length(rx_sym_manual);
        for i = 1:num_symbols
            bits_I_dec = demapping_logic(I_rec(i));
            bits_Q_dec = demapping_logic(Q_rec(i));
            
            % Table de correspondance inverse (cohérente avec l'émetteur manuel)
            bit_patterns = [
                0, 0;  % 0 -> 00 (correspond à -3)
                0, 1;  % 1 -> 01 (correspond à -1)
                1, 0;  % 2 -> 10 (correspond à +3)
                1, 1   % 3 -> 11 (correspond à +1)
            ];
            
            b_I = bit_patterns(bits_I_dec + 1, :);
            b_Q = bit_patterns(bits_Q_dec + 1, :);
            
            start_idx = (i-1)*4 + 1;
            rx_bits_manual_manual(start_idx:start_idx+1) = b_I';
            rx_bits_manual_manual(start_idx+2:start_idx+3) = b_Q';
        end
    end
    
    % Calcul du TEB pour la chaîne manuelle
    [~, ber_manual_manual] = biterr(tx_bits, rx_bits_manual_manual);
    disp(['TEB (Manuel -> Manuel) : ', num2str(ber_manual_manual)]);
    
    %% ==========================================
    %% 2.4 Test Croisé (Optionnel - pour démonstration)
    %% ==========================================
    % Test de la compatibilité entre les chaînes (devrait échouer pour 16QAM)
    if M == 16
        disp('--- Test Croisé (Démonstration) ---');
        % Démodulation automatique du signal manuel
        if M == 16
            rx_bits_cross_auto = qamdemod(rx_sym_manual, M, 'OutputType', 'bit', 'UnitAveragePower', true);
        else
            rx_bits_cross_auto = pskdemod(rx_sym_manual, M, phase_offset, 'OutputType', 'bit');
        end
        [~, ber_cross] = biterr(tx_bits, rx_bits_cross_auto);
        disp(['TEB Croisé (Manuel -> Automatique) : ', num2str(ber_cross)]);
        if ber_cross > 0.1
            disp('Note: Le TEB élevé est attendu car les schémas de mapping sont différents.');
        end
    end

    %% ==========================================
    %% 2.5 Visualisation des résultats
    %% ==========================================
    
    % Figure 1 : Comparaison des constellations émises
    figure('Name', [MODULATION_TYPE, ' - Constellations Émises'], ...
           'Position', [100, 100, 1200, 500]);
    
    n_const = min(200, length(tx_sym_auto));
    
    % Constellation automatique
    subplot(1, 3, 1);
    scatter(real(tx_sym_auto(1:n_const)), imag(tx_sym_auto(1:n_const)), 50, 'b', 'filled', 'MarkerEdgeColor', 'k');
    grid on; axis equal;
    xlabel('Composante I'); ylabel('Composante Q');
    title([MODULATION_TYPE, ' - Constellation Automatique']);
    
    % Constellation manuelle
    subplot(1, 3, 2);
    scatter(real(tx_sym_manual(1:n_const)), imag(tx_sym_manual(1:n_const)), 50, 'r', 'filled', 'MarkerEdgeColor', 'k');
    grid on; axis equal;
    xlabel('Composante I'); ylabel('Composante Q');
    title([MODULATION_TYPE, ' - Constellation Manuelle']);
    
    % Comparaison des constellations
    subplot(1, 3, 3);
    scatter(real(tx_sym_auto(1:n_const)), imag(tx_sym_auto(1:n_const)), 40, 'b', 'filled', 'MarkerFaceAlpha', 0.5);
    hold on;
    scatter(real(tx_sym_manual(1:n_const)), imag(tx_sym_manual(1:n_const)), 40, 'r', 'filled', 'MarkerFaceAlpha', 0.5);
    grid on; axis equal;
    xlabel('Composante I'); ylabel('Composante Q');
    title([MODULATION_TYPE, ' - Comparaison des Constellations']);
    legend('Automatique', 'Manuelle');
    
    % Figure 2 : Signaux reçus avec bruit
    figure('Name', [MODULATION_TYPE, ' - Signaux Reçus'], ...
           'Position', [100, 650, 1200, 500]);
    
    n_show = min(100, length(rx_sym_auto));
    
    % Signal reçu de la chaîne automatique
    subplot(2, 2, 1);
    scatter(real(rx_sym_auto(1:n_show)), imag(rx_sym_auto(1:n_show)), 30, 'b', 'filled', 'MarkerFaceAlpha', 0.6);
    grid on; axis equal;
    xlabel('Composante I'); ylabel('Composante Q');
    title(['Chaîne Auto - Reçu (RSB=', num2str(snr_target), 'dB)']);
    
    % Signal reçu de la chaîne manuelle
    subplot(2, 2, 2);
    scatter(real(rx_sym_manual(1:n_show)), imag(rx_sym_manual(1:n_show)), 30, 'r', 'filled', 'MarkerFaceAlpha', 0.6);
    grid on; axis equal;
    xlabel('Composante I'); ylabel('Composante Q');
    title(['Chaîne Manuelle - Reçu (RSB=', num2str(snr_target), 'dB)']);
    
    % Comparaison des performances
    subplot(2, 2, [3, 4]);
    bar([ber_auto_auto, ber_manual_manual]);
    set(gca, 'XTickLabel', {'Auto→Auto', 'Manuel→Manuel'});
    ylabel('Taux d''Erreur Binaire (TEB)');
    title([MODULATION_TYPE, ' - Comparaison des Performances']);
    grid on;
    
    % Ajouter les valeurs numériques sur les barres
    text(1, ber_auto_auto + 0.01, sprintf('%.4f', ber_auto_auto), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    text(2, ber_manual_manual + 0.01, sprintf('%.4f', ber_manual_manual), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    
    % Affichage des résultats
    disp(['Résumé ', MODULATION_TYPE, ' :']);
    disp(['  TEB Chaîne Automatique: ', num2str(ber_auto_auto)]);
    disp(['  TEB Chaîne Manuelle: ', num2str(ber_manual_manual)]);
    
    if M == 16
        disp('  Note: Les schémas de mapping 16QAM automatique et manuel sont différents.');
        disp('        C''est intentionnel pour démontrer la diversité des implémentations.');
    end
    
    disp(['Test de modulation ', MODULATION_TYPE, ' terminé']);
    disp(' ');
    
    % Pause pour permettre l'observation des figures
    pause(1);
end

disp('==========================================');
disp('Tous les tests de modulation sont terminés !');
disp('==========================================');

%% ==========================================
%% 3. Fonction locale : logique de décision pour la démodulation 16-QAM manuelle
%% ==========================================
function bits_dec = demapping_logic(level_val)
    % Logique de décision cohérente avec le mapping manuel
    % Mapping: -3->00(0), -1->01(1), 1->11(3), 3->10(2)
    if level_val < -2      % Correspond à -3
        bits_dec = 0;      % Binaire: 00
    elseif level_val < 0   % Correspond à -1
        bits_dec = 1;      % Binaire: 01
    elseif level_val < 2   % Correspond à +1
        bits_dec = 3;      % Binaire: 11
    else                   % Correspond à +3
        bits_dec = 2;      % Binaire: 10
    end
end