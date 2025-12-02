clc; clear; close all;

%% ==========================================
%% 1. Configuration des paramètres matériels WebLab (respect strict des spécifications)
%% ==========================================
fs_weblab = 200e6;           % Fréquence d'échantillonnage fixe WebLab : 200 MHz (immuable)
target_length = 40000;       % Longueur cible du signal : 40k échantillons (équilibre efficacité entraînement DPD et limites du buffer WebLab)
num_symbols = 4000;          % Nombre de symboles (la longueur finale est déterminée par le filtrage)
modulation_type = '16QAM';   % Sélection du type de modulation (BPSK/QPSK/16QAM)

% Paramètres du filtre RRC (assurer la propriété de Nyquist)
rolloff = 0.35;
span = 8;                    % Réduire span pour contrôler la longueur de sortie
sps = floor(fs_weblab / (fs_weblab/8));  % Calculer le facteur de suréchantillonnage (~8)

% Paramètres de back-off de puissance
backoff_db = 6;              % Back-off de 6dB pour éviter la saturation de l'AP
safety_margin = 0.95;        % Marge de sécurité de 5% pour éviter l'écrêtage matériel

%% ==========================================
%% 2. Génération du signal conforme aux spécifications matérielles WebLab
%% ==========================================
disp('==========================================');
disp('Génération du signal compatible matériel WebLab');
disp(['Type de modulation : ', modulation_type]);
disp(['Fréquence d''échantillonnage : ', num2str(fs_weblab/1e6), ' MHz']);
disp('==========================================');

% 2.1 Déterminer les paramètres selon le type de modulation
switch upper(modulation_type)
    case 'BPSK'
        M = 2; k = 1; phase_offset = 0;
    case 'QPSK'
        M = 4; k = 2; phase_offset = pi/4;
    case '16QAM'
        M = 16; k = 4; phase_offset = 0;
    otherwise
        error('Type de modulation inconnu');
end

% 2.2 Générer le flux binaire
num_bits_needed = num_symbols * k;
tx_bits = randi([0 1], num_bits_needed, 1);

% 2.3 Modulation (utiliser la toolbox pour garantir l'exactitude)
if M == 16
    tx_symbols = qammod(tx_bits, M, 'InputType', 'bit', 'UnitAveragePower', true);
else
    tx_symbols = pskmod(tx_bits, M, phase_offset, 'InputType', 'bit');
end

% 2.4 Conception du filtre RRC (contrôle de la longueur)
h_rrc = rcosdesign(rolloff, span, sps, 'sqrt');
filt_delay = span * sps;  % Délai introduit par le filtre

% 2.5 Filtrage de mise en forme d'impulsion
tx_filtered = upfirdn(tx_symbols, h_rrc, sps);

% 2.6 Tronquer à la longueur cible (supprimer la transition du filtre)
if length(tx_filtered) > target_length
    % Commencer la troncature après la stabilisation du filtre
    start_idx = filt_delay + 1;
    end_idx = start_idx + target_length - 1;
    tx_truncated = tx_filtered(start_idx:end_idx);
else
    % Si trop court, compléter avec des zéros (modifie les propriétés statistiques, non recommandé)
    warning('La longueur du signal est inférieure à la longueur cible, envisagez d''augmenter le nombre de symboles');
    tx_truncated = tx_filtered;
    tx_truncated(target_length) = 0;  % Compléter avec des zéros jusqu'à la longueur cible
end

% 2.7 Normalisation matérielle sécurisée (étape critique)
peak_amplitude = max(abs(tx_truncated));
if peak_amplitude == 0
    error('L''amplitude du signal est nulle');
end

% Double normalisation : marge de sécurité + back-off de puissance
scale_factor = safety_margin * 10^(-backoff_db/20) / peak_amplitude;
PAin = tx_truncated * scale_factor;

% Forcer l'amplitude maximale à ne pas dépasser 1 (protection matérielle)
PAin = PAin / max(abs(PAin)) * safety_margin;

% 2.8 S'assurer que la longueur est paire (exigence WebLab)
if mod(length(PAin), 2) ~= 0
    PAin = PAin(1:end-1);
end

%% ==========================================
%% 3. Vérification de la qualité du signal et de la compatibilité matérielle
%% ==========================================
disp('--- Vérification de la compatibilité matérielle ---');

% 3.1 Calcul des indicateurs clés
peak_val = max(abs(PAin));
rms_val = sqrt(mean(abs(PAin).^2));
PAPR_db = 20*log10(peak_val / rms_val);
avg_power_w = mean(abs(PAin).^2);
avg_power_dbm = 10*log10(avg_power_w/1e-3);

% 3.2 Calculer théoriquement le RMSin pouvant être défini
% RMSin de WebLab est la puissance d'entrée attendue (dBm), à calculer en fonction du signal réel
RMSin_estimate = 10*log10(rms_val^2/50) + 30;  % Conversion en valeur dBm sur une charge de 50 ohms

% 3.3 Vérification des limites matérielles WebLab
is_compliant = true;
verification_messages = {};

if length(PAin) < 1000
    verification_messages{end+1} = '❌ Signal trop court (<1000 échantillons)';
    is_compliant = false;
elseif length(PAin) > 1e6
    verification_messages{end+1} = '❌ Signal trop long (>1e6 échantillons)';
    is_compliant = false;
else
    verification_messages{end+1} = '✅ Longueur du signal conforme';
end

if PAPR_db > 20
    verification_messages{end+1} = '❌ PAPR trop élevé (>20dB)';
    is_compliant = false;
else
    verification_messages{end+1} = '✅ PAPR conforme';
end

if peak_val > 1
    verification_messages{end+1} = '❌ Pic dépassant 1 (risque d''écrêtage matériel)';
    is_compliant = false;
else
    verification_messages{end+1} = '✅ Amplitude crête conforme';
end

% Afficher les résultats de la vérification
for i = 1:length(verification_messages)
    disp(verification_messages{i});
end

% 3.4 Afficher les indicateurs techniques
disp('--- Indicateurs techniques du signal ---');
fprintf('Longueur finale : %d échantillons (%.2f µs)\n', length(PAin), length(PAin)/fs_weblab*1e6);
fprintf('Amplitude crête : %.4f\n', peak_val);
fprintf('Amplitude RMS : %.4f\n', rms_val);
fprintf('PAPR : %.2f dB\n', PAPR_db);
fprintf('Puissance moyenne : %.2f dBm (estimation)\n', avg_power_dbm);
fprintf('Réglage RMSin utilisable pour WebLab : ~%.2f dBm\n', RMSin_estimate);

if ~is_compliant
    warning('⚠️ Le signal n''est pas entièrement conforme aux spécifications WebLab, veuillez ajuster les paramètres');
end

%% ==========================================
%% 4. Vérification visuelle (pour le rapport technique)
%% ==========================================
figure('Position', [100, 100, 1400, 900]);

% 4.1 Diagramme de constellation (après normalisation)
subplot(2, 3, 1);
scatter(real(PAin(1:min(2000, end))), imag(PAin(1:min(2000, end))), ...
        20, 'b', 'filled', 'MarkerFaceAlpha', 0.6);
grid on; axis equal; xlim([-1.1 1.1]); ylim([-1.1 1.1]);
xlabel('En phase'); ylabel('En quadrature');
title(['Diagramme de constellation normalisé (', modulation_type, ')']);
annotation('textbox', [0.15, 0.85, 0.3, 0.05], 'String', ...
    sprintf('Pic: %.3f\nRMS: %.3f', peak_val, rms_val), ...
    'BackgroundColor', 'white', 'EdgeColor', 'black');

% 4.2 Enveloppe temporelle
subplot(2, 3, 2);
plot(1:min(1000, length(PAin)), abs(PAin(1:min(1000, length(PAin)))), ...
     'b-', 'LineWidth', 1.5);
hold on;
plot([1, min(1000, length(PAin))], [1, 1], 'r--', 'LineWidth', 1, ...
     'DisplayName', 'Limite matérielle');
plot([1, min(1000, length(PAin))], [rms_val, rms_val], 'g--', ...
     'LineWidth', 1, 'DisplayName', 'Valeur RMS');
grid on; xlabel('Indice d''échantillon'); ylabel('Amplitude');
title('Enveloppe temporelle (avec ligne de limite matérielle)');
legend('Enveloppe du signal', 'Limite matérielle (1.0)', 'Valeur RMS', 'Location', 'best');
ylim([0 1.1]);

% 4.3 Densité spectrale de puissance
subplot(2, 3, 3);
[pxx, f] = pwelch(PAin, 1024, 512, 1024, fs_weblab, 'centered');
plot(f/1e6, 10*log10(pxx), 'b-', 'LineWidth', 1.5);
grid on; xlabel('Fréquence (MHz)'); ylabel('Densité spectrale de puissance (dB/Hz)');
title('Spectre de puissance normalisé');
xlim([-fs_weblab/2e6, fs_weblab/2e6]);
annotation('textbox', [0.72, 0.85, 0.2, 0.05], 'String', ...
    sprintf('Fréq. échant.: %.0f MHz\nBande passante: ~%.1f MHz', fs_weblab/1e6, fs_weblab/sps/1e6), ...
    'BackgroundColor', 'white', 'EdgeColor', 'black');

% 4.4 Statistiques de distribution d'amplitude
subplot(2, 3, 4);
[counts, edges] = histcounts(abs(PAin), 50);
bar(edges(1:end-1), counts, 'FaceColor', 'b', 'EdgeColor', 'none');
grid on; xlabel('Valeur d''amplitude'); ylabel('Fréquence');
title('Histogramme de distribution d''amplitude');
xlim([0 1]);
annotation('textbox', [0.15, 0.38, 0.3, 0.05], 'String', ...
    sprintf('PAPR: %.1f dB\nBack-off: %.0f dB', PAPR_db, backoff_db), ...
    'BackgroundColor', 'white', 'EdgeColor', 'black');

% 4.5 Formes d'onde temporelles des composantes I/Q
subplot(2, 3, 5);
t_show = 1:min(200, length(PAin));
plot(t_show, real(PAin(t_show)), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Composante I');
hold on;
plot(t_show, imag(PAin(t_show)), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Composante Q');
grid on; xlabel('Indice d''échantillon'); ylabel('Amplitude');
title('Formes d''onde temporelles des composantes I/Q');
legend('Location', 'best');
ylim([-1.1 1.1]);

% 4.6 Tableau récapitulatif des indicateurs
subplot(2, 3, 6);
axis off;
text(0.1, 0.9, '🚀 Récapitulatif des indicateurs techniques du signal WebLab', 'FontSize', 12, 'FontWeight', 'bold');
text(0.1, 0.8, sprintf('Type de modulation : %s', modulation_type), 'FontSize', 10);
text(0.1, 0.75, sprintf('Fréquence d''échantillonnage : %.0f MHz', fs_weblab/1e6), 'FontSize', 10);
text(0.1, 0.7, sprintf('Longueur du signal : %d échantillons', length(PAin)), 'FontSize', 10);
text(0.1, 0.65, sprintf('Durée : %.2f µs', length(PAin)/fs_weblab*1e6), 'FontSize', 10);
text(0.1, 0.6, sprintf('Amplitude crête : %.4f', peak_val), 'FontSize', 10);
text(0.1, 0.55, sprintf('PAPR : %.2f dB', PAPR_db), 'FontSize', 10);
text(0.1, 0.5, sprintf('Back-off de puissance : %.0f dB', backoff_db), 'FontSize', 10);
text(0.1, 0.45, sprintf('Marge de sécurité : %.0f%%', (1-safety_margin)*100), 'FontSize', 10);
text(0.1, 0.4, 'Conformité matérielle :', 'FontSize', 10, 'FontWeight', 'bold');
if is_compliant
    text(0.1, 0.35, '✅ Entièrement conforme aux spécifications WebLab', 'FontSize', 10, 'Color', 'green');
else
    text(0.1, 0.35, '⚠️ Paramètres à ajuster', 'FontSize', 10, 'Color', 'red');
end

%% ==========================================
%% 5. Sauvegarde au format compatible WebLab
%% ==========================================
% 5.1 Créer la structure d'information du signal
signal_info = struct();
signal_info.modulation = modulation_type;
signal_info.fs = fs_weblab;
signal_info.num_samples = length(PAin);
signal_info.peak_amplitude = peak_val;
signal_info.rms_amplitude = rms_val;
signal_info.PAPR_db = PAPR_db;
signal_info.backoff_db = backoff_db;
signal_info.safety_margin = safety_margin;
signal_info.generation_date = datestr(now, 'yyyy-mm-dd HH:MM:SS');
signal_info.rmsin_estimate = RMSin_estimate;

% 5.2 Sauvegarder le fichier (peut être chargé directement par le script principal WebLab)
save('WebLab_PAin_16QAM.mat', 'PAin', 'signal_info', '-v7.3');
disp(' ');
disp('📁 Signal sauvegardé sous : WebLab_PAin_16QAM.mat');

% 5.3 Générer les instructions d'utilisation
disp(' ');
disp('🔧 Guide d''utilisation WebLab :');
disp('   1. Dans le script principal WebLab main.m :');
disp('      load(''WebLab_PAin_16QAM.mat'');');
disp('   2. Définir le paramètre RMSin (suggestion : commencer à -15 dBm) :');
disp(sprintf('      RMSin = %.1f; %% dBm (à ajuster selon l''AP réel)', RMSin_estimate));
disp('   3. Appeler la fonction de mesure :');
disp('      [PAout, RMSout, Idc, Vdc] = RFWebLab_PA_meas_v1_2(PAin, RMSin);');
disp('   4. Aligner les signaux (étape clé) :');
disp('      PAout_aligned = timealign(PAin, PAout);');
disp(' ');
disp('✅ Tâche 3 - Phase 1 terminée !');
disp('==========================================');