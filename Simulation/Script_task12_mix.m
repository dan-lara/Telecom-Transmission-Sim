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
N = 30000;      % Nombre de bits initial

% Configuration du filtre (commune à toutes les modulations)
rolloff = 0.35; span = 10;
h_rrc = rcosdesign(rolloff, span, sps, "sqrt");
h_rrc = h_rrc / norm(h_rrc);    % Énergie de filtre normalisée
total_delay = span * sps;

% Paramètres des canaux dégradés
CHANNEL_IMPAIRMENTS = {'AWGN', 'Impulsive', 'Rayleigh', 'All'};
USE_IMPAIRMENTS = true;  % Activer/désactiver les canaux dégradés
SELECTED_IMPAIRMENT = 'Rayleigh';  % Choisir 'AWGN', 'Impulsive', 'Rayleigh', ou 'All'

% Paramètres du bruit impulsionnel
IMPULSE_PROBABILITY = 0.01;     % Probabilité d'impulsion
IMPULSE_AMPLITUDE = 5;          % Amplitude du bruit impulsionnel (multiplicateur)

% Paramètres du canal de Rayleigh (fading)
NUM_RAYS = 3;                   % Nombre de trajets multiples
MAX_DELAY = 5;                  % Délai maximum en échantillons (trajets multiples)
DOPPLER_SHIFT = 10;             % Décalage Doppler en Hz (pour fading temporel)

% SNR pour test fixe
snr_target_fixed = 30; % dB

%% ==========================================
%% 2. Fonctions auxiliaires pour les canaux dégradés
%% ==========================================

% Fonction pour générer du bruit impulsionnel
function noise = impulsive_noise(signal_length, probability, amplitude)
    noise = zeros(signal_length, 1);
    impulse_indices = find(rand(signal_length, 1) < probability);
    if ~isempty(impulse_indices)
        % Générer des impulsions de signe aléatoire
        noise(impulse_indices) = amplitude * sign(randn(length(impulse_indices), 1)) .* ...
                                 (1 + 0.5 * randn(length(impulse_indices), 1));
    end
end

% % Fonction pour simuler un canal de Rayleigh (fading multi-trajets)
% function [output, channel_coeffs] = rayleigh_fading_channel(input, num_rays, max_delay, doppler_shift, fs)
%     % Générer les coefficients du canal
%     channel_coeffs = zeros(num_rays, 1);
%     delays = randi([0, max_delay], num_rays, 1);
%     
%     % Trier les retards
%     [delays, sort_idx] = sort(delays);
%     
%     % Générer des gains complexes (distribution de Rayleigh)
%     for i = 1:num_rays
%         % Gain complexe avec amplitude Rayleigh et phase uniforme
%         amplitude = sqrt(0.5) * randn();
%         phase = 2*pi*rand();
%         channel_coeffs(i) = amplitude * exp(1j*phase);
%     end
%     
%     % Normaliser l'énergie des coefficients
%     channel_coeffs = channel_coeffs / sqrt(sum(abs(channel_coeffs).^2));
%     
%     % Appliquer le fading Doppler (simplifié)
%     t = (0:length(input)-1)' / fs;
%     doppler_effect = exp(1j*2*pi*doppler_shift*t);
%     
%     % Convolution avec le canal multi-trajets
%     output = zeros(size(input));
%     for i = 1:num_rays
%         delay = delays(i);
%         coeff = channel_coeffs(i);
%         if delay == 0
%             output = output + coeff * input .* doppler_effect;
%         else
%             output(delay+1:end) = output(delay+1:end) + coeff * input(1:end-delay) .* doppler_effect(1:end-delay);
%         end
%     end
%     
%     % Normaliser la puissance de sortie
%     output = output * sqrt(mean(abs(input).^2) / mean(abs(output).^2));
% end

%% ==========================================
%% 3. Fonction de simulation de chaîne avec canaux dégradés
%% ==========================================

function [ber, rx_sym, rx_filt] = simulate_chain_with_impairments(tx_symbols, tx_bits_with_pilot, M, phase_offset, ...
                                        snr_target, k, sps, fc, fs, h_rrc, total_delay, ...
                                        chain_type, impairment_type, ...
                                        impulse_prob, impulse_amp, num_rays, max_delay, doppler_shift, SELECTED_IMPAIRMENT)
    % Paramètres supplémentaires pour les canaux dégradés
    
    % Filtrage RRC
    tx_filtered = upfirdn(tx_symbols, h_rrc, sps);
    
    % Modulation
    t = (0 : length(tx_filtered) - 1)' / fs;
    tx_rf = real(tx_filtered).*cos(2*pi*fc*t) - imag(tx_filtered).*sin(2*pi*fc*t);
    
    % Canal avec différentes dégradations
    snr_val = snr_target + 10*log10(k) - 10*log10(sps);

    %rx_rf = tx_rf;
    
    switch impairment_type
        case 'AWGN'
            % Seulement AWGN
            rx_rf = awgn(tx_rf, snr_val, 'measured');
            
        case 'Impulsive'
            % AWGN + bruit impulsionnel
            rx_rf = awgn(tx_rf, snr_val, 'measured');
            impulse_noise = impulsive_noise(length(rx_rf), impulse_prob, impulse_amp);
            rx_rf = rx_rf + impulse_noise;
            
        case 'Rayleigh'
            % Canal de Rayleigh (fading) + AWGN
            rayChan = comm.RayleighChannel(...
                'SampleRate', fs, ...
                'PathDelays', linspace(0, max_delay/fs, num_rays), ...  % 均匀分布时延
                'AveragePathGains', -linspace(0, 10, num_rays), ...     % 指数衰减功率（dB）
                'MaximumDopplerShift', doppler_shift, ...
                'DopplerSpectrum', doppler('Jakes'), ...
                'PathGainsOutputPort', true);
            
            [rx_rf_faded, channel_coeffs] = rayChan(tx_rf);
            
            % Normaliser la puissance de sortie pour maintenir l'énergie du signal
            rx_rf_faded = rx_rf_faded * sqrt(mean(abs(tx_rf).^2) / mean(abs(rx_rf_faded).^2));
            
            % Ajouter AWGN
            rx_rf = awgn(rx_rf_faded, snr_val, 'measured');
            
        case 'All'
            % Combinaison de toutes les dégradations
            % 1. Appliquer le fading avec comm.RayleighChannel
            rayChan = comm.RayleighChannel(...
                'SampleRate', fs, ...
                'PathDelays', linspace(0, max_delay/fs, num_rays), ...  % 均匀分布时延
                'AveragePathGains', -linspace(0, 10, num_rays), ...     % 指数衰减功率（dB）
                'MaximumDopplerShift', doppler_shift, ...
                'DopplerSpectrum', doppler('Jakes'), ...
                'PathGainsOutputPort', true);
            
            [rx_rf_faded, channel_coeffs] = rayChan(tx_rf);
            
            % Normaliser la puissance de sortie
            rx_rf_faded = rx_rf_faded * sqrt(mean(abs(tx_rf).^2) / mean(abs(rx_rf_faded).^2));
            
            % 2. Ajouter AWGN
            rx_rf = awgn(rx_rf_faded, snr_val, 'measured');
            % 3. Ajouter le bruit impulsionnel
            impulse_noise = impulsive_noise(length(rx_rf), impulse_prob, impulse_amp);
            rx_rf = rx_rf + impulse_noise;
            
        otherwise
            % Par défaut, seulement AWGN
            rx_rf = awgn(tx_rf, snr_val, 'measured');
    end
    
    % Réception
    rx_base = rx_rf .* (cos(2*pi*fc*t) - 1j*sin(2*pi*fc*t)) * 2;

    % === Traitement du bruit impulsionnel (Clipping) ===
    
    % 1. Calculer l'amplitude moyenne (RMS) du signal courant.
    current_rms = sqrt(mean(abs(rx_base).^2));
    
    % 2. Réglage du seuil (Threshold)
    threshold = 4 * current_rms; 
    
    % 3. Limitation (Clipping)
    outlier_indices = abs(rx_base) > threshold;
    
    % Méthode B: Mise à zéro (Blanking)
    rx_base(outlier_indices) = 0; 

    rx_filt = upfirdn(rx_base, h_rrc, 1, 1);
    
    % Trouver le point de départ des symboles
    start_idx = total_delay + 1;
    rx_sym = rx_filt(start_idx : sps : end);
    
    % === Égaliseur LMS (pour contrer les effets de Rayleigh et multi-trajets) ===
    if ismember(impairment_type, {'Rayleigh', 'All'})
        disp('Application de l''égaliseur LMS pour corriger le fading Rayleigh...');
        
        % 1. Génération dynamique de la constellation
        if M == 16
            eq_const = qammod(0:M-1, M, 'UnitAveragePower', true);
        else
            eq_const = pskmod(0:M-1, M, phase_offset);
        end

        % 2. Calcul dynamique de la longueur de comparaison
        safe_len = min([length(rx_sym), length(tx_symbols), 1000]);
        
        if safe_len > 10
            % Alignement au niveau des symboles
            d_sym = finddelay(tx_symbols(1:safe_len), rx_sym(1:safe_len));
            
            % Compensation du délai
            if d_sym > 0 && d_sym < length(rx_sym)
                rx_in = rx_sym(d_sym + 1 : end);
                tx_ref = tx_symbols(1 : length(rx_in));
            else
                rx_in = rx_sym;
                tx_ref = tx_symbols;
            end
        else
            rx_in = rx_sym;
            tx_ref = tx_symbols;
        end

        % 3. Configuration de l'égaliseur avec paramètres optimisés
        lms_eq = comm.LinearEqualizer(...
            'Algorithm', 'LMS', ...
            'NumTaps', 61, ... % Augmenté pour meilleure performance
            'StepSize', 0.0005, ... % Réduit pour convergence plus stable
            'ReferenceTap', 31, ...
            'Constellation', complex(eq_const), ...
            'AdaptAfterTraining', true);
        
        % 4. Exécution de l'égalisation
        train_len = min(length(rx_in), length(tx_ref));
        if train_len > 100  
            [rx_eq, ~] = lms_eq(rx_in, tx_ref);
        else
            rx_eq = rx_in; 
        end

        % 5. Correction de phase pour éviter les inversions
        if length(rx_eq) > 50 && length(tx_ref) > 50
            align_len = min([1000, length(rx_eq), length(tx_ref)]);
            phase_rot = angle(mean(conj(rx_eq(1:align_len)) .* tx_ref(1:align_len)));
            rx_eq = rx_eq * exp(-1j * phase_rot);
        else
            disp('    Symboles insuffisants pour une correction de phase précise.');
        end
        
        % 6. Troncature adaptative
        converge_skip = 500;
        if length(rx_eq) > converge_skip + 1000
            rx_sym = rx_eq(converge_skip + 1 : end);
        else
            rx_sym = rx_eq;
        end
    end

    % Démodulation (automatique ou manuelle)
    if isempty(rx_sym)
        disp('Avertissement : vecteur de symboles reçus vide, saut de la démodulation.');
        rx_bits = [];
    else
        if strcmp(chain_type, 'auto')
            if M == 16
                rx_bits = qamdemod(rx_sym, M, 'OutputType', 'bit', 'UnitAveragePower', true);
            else
                rx_bits = pskdemod(rx_sym, M, phase_offset, 'OutputType', 'bit');
            end
        else
            % Démodulation manuelle
            rx_bits = demapping_manual(rx_sym, M, phase_offset, k);
        end
    end

    % Calcul du BER avec alignement robuste
    total_rx_bits = length(rx_bits);
    total_tx_bits = length(tx_bits_with_pilot);
    
    if total_rx_bits < 10
        ber = 0.5;
        num_errors = floor(total_tx_bits / 2);
    else
        % 1. Recherche du délai au niveau des bits
        bit_comp_len = min([total_tx_bits, total_rx_bits, 2000]);
        d_bit = finddelay(tx_bits_with_pilot(1:bit_comp_len), rx_bits(1:bit_comp_len));
        
        % 2. Alignement élastique
        if d_bit > 0 && d_bit < total_rx_bits
            rx_temp = rx_bits(d_bit + 1 : end);
            tx_temp = tx_bits_with_pilot;
        elseif d_bit < 0 && abs(d_bit) < total_tx_bits
            tx_temp = tx_bits_with_pilot(-d_bit + 1 : end);
            rx_temp = rx_bits;
        else
            tx_temp = tx_bits_with_pilot;
            rx_temp = rx_bits;
        end
        
        % 3. Alignement final avec vérification des limites
        L = min(length(tx_temp), length(rx_temp));
        if L > 0
            tx_aligned = tx_temp(1:L);
            rx_aligned = rx_temp(1:L);
            [num_errors, ber] = biterr(tx_aligned, rx_aligned);
        else
            ber = 0.5;
        end
    end
    
    disp(['BER: ', num2str(ber)]);
end

%% ==========================================
%% 4. Fonction pour tracer l'oeil (eye diagram)
%% ==========================================

function plot_eye_diagram(signal, sps, modulation_name, fs, signal_type)
    % signal: signal en bande de base (surechantillonné)
    % sps: échantillons par symbole
    % signal_type: 'baseband' ou 'filtered'
    
    if nargin < 5
        signal_type = 'baseband';
    end
    
    % S'assurer que le signal est suffisamment long
    if length(signal) < 3 * sps
        fprintf('Signal trop court pour le diagramme de l''oeil. Longueur: %d, besoin au moins: %d\n', ...
                length(signal), 3*sps);
        return;
    end
    
    % Tronquer le signal pour avoir un nombre entier de périodes d'oeil
    num_symbols = floor(length(signal) / sps);
    num_eyes = num_symbols - 1;  % Pour avoir des tracés de 2 symboles
    
    if num_eyes < 2
        fprintf('Pas assez de symboles pour tracer un diagramme de l''oeil.\n');
        return;
    end
    
    % Créer une figure pour l'oeil
    figure('Name', [modulation_name, ' - Diagramme de l''Oeil (', signal_type, ')'], ...
           'Position', [100, 100, 900, 700]);
    
    % Déterminer le nombre de tracés à afficher (limité pour la clarté)
    num_traces = min(50, num_eyes);
    
    % Préparer les données pour le tracé
    eye_data_I = zeros(num_traces, 2*sps);
    eye_data_Q = zeros(num_traces, 2*sps);
    
    for i = 1:num_traces
        start_idx = (i-1)*sps + 1;
        end_idx = start_idx + 2*sps - 1;
        if end_idx <= length(signal)
            eye_data_I(i, :) = real(signal(start_idx:end_idx));
            eye_data_Q(i, :) = imag(signal(start_idx:end_idx));
        end
    end
    
    % Tracer l'oeil pour la partie réelle (I)
    subplot(2, 2, 1);
    hold on;
    time_eye = (0:2*sps-1) / fs * 1000;  % en ms
    for i = 1:num_traces
        plot(time_eye, eye_data_I(i, :), 'b-', 'LineWidth', 0.5, 'Color', [0, 0, 1, 0.3]);
    end
    
    % Ajouter la moyenne
    mean_I = mean(eye_data_I, 1, 'omitnan');
    plot(time_eye, mean_I, 'k-', 'LineWidth', 2);
    
    grid on;
    xlabel('Temps (ms)', 'FontSize', 10);
    ylabel('Amplitude I', 'FontSize', 10);
    title([modulation_name, ' - Oeil (Partie Réelle I)'], 'FontSize', 12);
    
    % Tracer l'oeil pour la partie imaginaire (Q)
    subplot(2, 2, 2);
    hold on;
    for i = 1:num_traces
        plot(time_eye, eye_data_Q(i, :), 'r-', 'LineWidth', 0.5, 'Color', [1, 0, 0, 0.3]);
    end
    
    % Ajouter la moyenne
    mean_Q = mean(eye_data_Q, 1, 'omitnan');
    plot(time_eye, mean_Q, 'k-', 'LineWidth', 2);
    
    grid on;
    xlabel('Temps (ms)', 'FontSize', 10);
    ylabel('Amplitude Q', 'FontSize', 10);
    title([modulation_name, ' - Oeil (Partie Imaginaire Q)'], 'FontSize', 12);
    
    % Tracer l'oeil en 2D (constellation temporelle)
    subplot(2, 2, [3, 4]);
    hold on;
    
    for i = 1:min(30, num_traces)
        plot(eye_data_I(i, :), eye_data_Q(i, :), ...
             '-', 'LineWidth', 0.5, 'Color', [0.5, 0.5, 0.5, 0.3]);
    end
    
    % Ajouter les points de décision si applicable
    if contains(modulation_name, 'BPSK')
        scatter([-1, 1], [0, 0], 100, 'r', 'x', 'LineWidth', 2);
    elseif contains(modulation_name, 'QPSK')
        scatter([-1, 1, -1, 1], [-1, -1, 1, 1], 100, 'r', 'x', 'LineWidth', 2);
    elseif contains(modulation_name, '16QAM')
        % Points pour 16QAM
        levels = [-3, -1, 1, 3] / sqrt(10);
        [X, Y] = meshgrid(levels, levels);
        scatter(X(:), Y(:), 80, 'r', 'x', 'LineWidth', 1.5);
    end
    
    grid on;
    axis equal;
    xlabel('Amplitude I', 'FontSize', 10);
    ylabel('Amplitude Q', 'FontSize', 10);
    title([modulation_name, ' - Diagramme de l''Oeil 2D'], 'FontSize', 12);
    
    % Calculer et afficher les métriques de l'oeil
    eye_samples = round(sps/2);  % Point central de l'oeil
    
    if eye_samples <= size(eye_data_I, 2)
        % Extraire les valeurs au point central
        center_I = eye_data_I(:, eye_samples);
        center_Q = eye_data_Q(:, eye_samples);
        
        % Calcul de l'ouverture verticale
        if ~isempty(center_I) && length(center_I) > 1
            % Pour BPSK et QPSK
            if contains(modulation_name, 'BPSK')
                opening_I = min(abs(diff(sort(center_I))));
                opening_Q = 0;  % Pour BPSK, Q est nul
            else
                opening_I = min(abs(diff(sort(center_I))));
                opening_Q = min(abs(diff(sort(center_Q))));
            end
            
            % Calcul du bruit d'ouverture
            noise_I = std(center_I, 'omitnan');
            noise_Q = std(center_Q, 'omitnan');
            
            % Afficher les métriques
            text(0.05, 0.95, sprintf('Ouverture I: %.4f', opening_I), ...
                 'Units', 'normalized', 'FontSize', 10, 'BackgroundColor', 'white');
            text(0.05, 0.90, sprintf('Ouverture Q: %.4f', opening_Q), ...
                 'Units', 'normalized', 'FontSize', 10, 'BackgroundColor', 'white');
            text(0.05, 0.85, sprintf('Bruit I: %.4f', noise_I), ...
                 'Units', 'normalized', 'FontSize', 10, 'BackgroundColor', 'white');
            text(0.05, 0.80, sprintf('Bruit Q: %.4f', noise_Q), ...
                 'Units', 'normalized', 'FontSize', 10, 'BackgroundColor', 'white');
            
            % Calcul du rapport d'ouverture (Eye Opening Factor)
            if opening_I > 0
                eof_I = opening_I / (2 * noise_I);
                text(0.05, 0.75, sprintf('EOF I: %.2f', eof_I), ...
                     'Units', 'normalized', 'FontSize', 10, 'BackgroundColor', 'white');
            end
        end
    end
    
    % Ajouter une grille de temps sur les axes x et y
    grid on;
    box on;
    
    % Ajuster la vue pour mieux voir l'oeil
    if contains(modulation_name, 'BPSK')
        xlim([-1.5, 1.5]);
        ylim([-1.5, 1.5]);
    end
    
    fprintf('Diagramme de l''oeil généré pour %s avec %d tracés\n', modulation_name, num_traces);
end

%% ==========================================
%% 5. Fonction pour tracer l'oeil sur le signal filtré en émission
%% ==========================================

function plot_eye_diagram_tx(tx_filtered, sps, modulation_name, fs)
    % Cette fonction trace l'oeil sur le signal filtré en émission
    % tx_filtered: signal après filtrage RRC en émission
    
    if length(tx_filtered) < 3 * sps
        fprintf('Signal TX trop court pour le diagramme de l''oeil.\n');
        return;
    end
    
    % Créer une figure séparée pour l'oeil en émission
    figure('Name', [modulation_name, ' - Oeil en Émission'], ...
           'Position', [100, 100, 900, 400]);
    
    % Tracer l'oeil pour la partie réelle
    subplot(1, 2, 1);
    hold on;
    
    num_traces = min(30, floor(length(tx_filtered) / sps) - 1);
    
    for i = 1:num_traces
        start_idx = (i-1)*sps + 1;
        end_idx = start_idx + 2*sps - 1;
        if end_idx <= length(tx_filtered)
            time_eye = (0:2*sps-1) / fs * 1000;
            plot(time_eye, real(tx_filtered(start_idx:end_idx)), 'b-', 'LineWidth', 0.5);
        end
    end
    
    grid on;
    xlabel('Temps (ms)');
    ylabel('Amplitude');
    title([modulation_name, ' - Oeil en Émission (Partie Réelle)']);
    
    % Tracer l'oeil pour la partie imaginaire
    subplot(1, 2, 2);
    hold on;
    
    for i = 1:num_traces
        start_idx = (i-1)*sps + 1;
        end_idx = start_idx + 2*sps - 1;
        if end_idx <= length(tx_filtered)
            time_eye = (0:2*sps-1) / fs * 1000;
            plot(time_eye, imag(tx_filtered(start_idx:end_idx)), 'r-', 'LineWidth', 0.5);
        end
    end
    
    grid on;
    xlabel('Temps (ms)');
    ylabel('Amplitude');
    title([modulation_name, ' - Oeil en Émission (Partie Imaginaire)']);
    
    sgtitle([modulation_name, ' - Diagramme de l''Oeil en Émission'], 'FontSize', 14);
end

%% ==========================================
%% 6. Test avec différents canaux dégradés
%% ==========================================
if USE_IMPAIRMENTS
    disp('==========================================');
    disp('Test avec canaux dégradés');
    disp('==========================================');
    
    % Définir les types de dégradation à tester
    if strcmp(SELECTED_IMPAIRMENT, 'All')
        impairments_to_test = {'AWGN', 'Impulsive', 'Rayleigh', 'All'};
    else
        impairments_to_test = {SELECTED_IMPAIRMENT};
    end
    
    % Boucle sur les modulations
    for mod_idx = 1:length(MODULATION_TYPES)
        MODULATION_TYPE = MODULATION_TYPES{mod_idx};
        disp(['Modulation : ', MODULATION_TYPE]);
        
        % Configuration des paramètres de modulation
        switch MODULATION_TYPE
            case 'BPSK'
                M = 2; k = 1; phase_offset = 0;
            case 'QPSK'
                M = 4; k = 2; phase_offset = pi/4;
            case '16QAM'
                M = 16; k = 4; phase_offset = 0;
        end
        
        N_test = N;  % Réduire pour accélérer le test
        N_adjusted = floor(N_test / k) * k;
        tx_bits = randi([0 1], N_adjusted, 1);
        
        pilot_len = 3000; % Longueur du préambule de guidage
        if M == 16
            pilot_indices = randi([0 M-1], pilot_len, 1);
            pilot_sym = qammod(pilot_indices, M, 'UnitAveragePower', true);
            pilot_bits = qamdemod(pilot_sym, M, 'OutputType', 'bit', 'UnitAveragePower', true);
        else
            pilot_indices = randi([0 M-1], pilot_len, 1);
            pilot_sym = pskmod(pilot_indices, M, phase_offset);
            pilot_bits = pskdemod(pilot_sym, M, phase_offset, 'OutputType', 'bit');
        end
        
        
        if M == 16
            tx_sym_auto = qammod(tx_bits, M, 'InputType', 'bit', 'UnitAveragePower', true);
        else
            tx_sym_auto = pskmod(tx_bits, M, phase_offset, 'InputType', 'bit');
        end
        tx_sym_with_pilot = [pilot_sym; tx_sym_auto];
        tx_bits_with_pilot = [pilot_bits; tx_bits];
        
        disp(['    Longueur totale: ', num2str(length(tx_sym_with_pilot)), ' symboles (', ...
              num2str(pilot_len), ' guide频 + ', num2str(length(tx_sym_auto)), ' données)']);
        
        % Tester chaque type de dégradation
        ber_results = zeros(length(impairments_to_test), 1);
        
        for impair_idx = 1:length(impairments_to_test)
            impairment = impairments_to_test{impair_idx};
            disp(['  Dégradation : ', impairment]);
            
            % Simulation avec la dégradation
            [ber, rx_sym, rx_filt] = simulate_chain_with_impairments(...
                tx_sym_with_pilot, tx_bits_with_pilot, M, phase_offset, ...
                snr_target_fixed, k, sps, fc, fs, h_rrc, total_delay, ...
                'auto', impairment, ...
                IMPULSE_PROBABILITY, IMPULSE_AMPLITUDE, NUM_RAYS, MAX_DELAY, DOPPLER_SHIFT, SELECTED_IMPAIRMENT);
            
            ber_results(impair_idx) = ber;
            
            % Visualisation pour chaque type de dégradation
            figure('Name', [MODULATION_TYPE, ' - Canal ', impairment], ...
                   'Position', [100, 100, 1400, 800]);
            
            % Sous-graphique 1: Constellation
            subplot(2, 3, 1);
            scatter(real(rx_sym), imag(rx_sym), 30, 'b', 'filled', 'MarkerFaceAlpha', 0.6);
            grid on; axis equal;
            xlabel('I'); ylabel('Q');
            title(['Constellation - ', impairment]);
            
            % Sous-graphique 2: Signal temporel (réel et imaginaire)
            subplot(2, 3, 2);
            plot(real(rx_sym(1:min(200, length(rx_sym)))), 'b-', 'LineWidth', 1); hold on;
            plot(imag(rx_sym(1:min(200, length(rx_sym)))), 'r--', 'LineWidth', 1);
            grid on; xlabel('Symbole'); ylabel('Amplitude');
            title(['Signal Temporel - ', impairment]);
            legend('I', 'Q');
            
            % Sous-graphique 3: BER
            subplot(2, 3, 3);
            bar(ber);
            ylabel('BER');
            title(['BER = ', num2str(ber, '%.4f')]);
            grid on;
            
            eye_pilot_len = min(100, pilot_len);
            eye_tx_sym = tx_sym_with_pilot(1:eye_pilot_len);
            eye_tx_bits = tx_bits_with_pilot(1:eye_pilot_len * k);
            subplot(2, 3, [4, 5, 6]);
            
            % Recalculer le signal reçu pour l'oeil
            [~, ~, rx_filt_complete] = simulate_chain_with_impairments(...
                eye_tx_sym, eye_tx_bits, M, phase_offset, snr_target_fixed, k, sps, fc, fs, h_rrc, total_delay, ...
                'auto', impairment, ...
                IMPULSE_PROBABILITY, IMPULSE_AMPLITUDE, NUM_RAYS, MAX_DELAY, DOPPLER_SHIFT, SELECTED_IMPAIRMENT);
            
            % Tracer l'oeil
            if ~isempty(rx_filt_complete) && length(rx_filt_complete) > 3*sps
                plot_eye_diagram(rx_filt_complete(total_delay+1:end), sps, ...
                                [MODULATION_TYPE, ' - ', impairment], fs, 'filtre réception');
            else
                text(0.5, 0.5, 'Signal trop court pour l''oeil', ...
                     'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', 12);
            end
            
            % Ajuster le titre principal
            sgtitle([MODULATION_TYPE, ' - Canal ', impairment, ' (SNR=', num2str(snr_target_fixed), 'dB)'], 'FontSize', 14);
            
            disp(['    BER = ', num2str(ber)]);
            
            % Figure supplémentaire: Oeil en émission
            figure('Name', [MODULATION_TYPE, ' - Oeil en Émission']);
            tx_filtered_eye = upfirdn(eye_tx_sym, h_rrc, sps);
            plot_eye_diagram_tx(tx_filtered_eye, sps, MODULATION_TYPE, fs);
        end
        
        % Figure comparative des BER pour différentes dégradations
        figure('Name', [MODULATION_TYPE, ' - Comparaison des Canaux'], ...
               'Position', [100, 100, 800, 600]);
        
        bar(ber_results);
        set(gca, 'XTickLabel', impairments_to_test);
        ylabel('Bit Error Rate (BER)');
        title([MODULATION_TYPE, ' - Impact des Différents Canaux (SNR=', num2str(snr_target_fixed), 'dB)']);
        grid on;
        
        % Ajouter les valeurs sur les barres
        for i = 1:length(ber_results)
            if ber_results(i) > 0
                text(i, ber_results(i) + max(ber_results)*0.1, sprintf('%.4f', ber_results(i)), ...
                    'HorizontalAlignment', 'center', 'FontWeight', 'bold');
            end
        end
        
        disp(' ');
    end
end

%% ==========================================
%% 7. Fonctions auxiliaires
%% ==========================================

% Fonction de démodulation manuelle
function rx_bits = demapping_manual(rx_sym, M, phase_offset, k)
    num_symbols = length(rx_sym);
    rx_bits = zeros(num_symbols * k, 1);
    
    if M == 2 % BPSK
        % Décision basée sur la partie réelle
        decisions = real(rx_sym) < 0;
        rx_bits = double(decisions);
        
    elseif M == 4 % QPSK
        % Décision par quadrant
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
        
        % === 修改开始 ===
        % 1. 动态获取当前接收到的符号数量 (因为均衡器可能截断了前几百个点)
        current_num_symbols = length(rx_sym); 

        % 2. 重新初始化 rx_bits 容器 (长度也要变短)
        rx_bits = zeros(current_num_symbols * 4, 1); % 16QAM 是 *4，QPSK 是 *2

        % 3. 循环次数改为 current_num_symbols
        for i = 1:current_num_symbols 
            bits_I_dec = demapping_logic_16qam(I_rec(i));
            bits_Q_dec = demapping_logic_16qam(Q_rec(i));

            b_I = bit_patterns(bits_I_dec + 1, :);
            b_Q = bit_patterns(bits_Q_dec + 1, :);

            start_idx = (i-1)*4 + 1;
            rx_bits(start_idx:start_idx+1) = b_I';
            rx_bits(start_idx+2:start_idx+3) = b_Q';
        end
        % === 修改结束 ===
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