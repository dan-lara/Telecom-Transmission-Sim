clc; clear; close all;

%% ==========================================
%% 1. CONFIGURATION DES PARAMÈTRES POUR WEBLAB
%% ==========================================

% Fréquence d'échantillonnage fixe du WebLab (200 MHz)
fs_weblab = 200e6;

% Longueur du signal : 10,000 échantillons (suffisant pour DPD)
target_length = 10000;

% Choix de la modulation - 16QAM recommandée pour la caractérisation DPD
modulation_type = '16QAM';

% Paramètres du filtre RRC (Root Raised Cosine)
rolloff = 0.25;      % Facteur de retombée (standard LTE)
span = 10;           % Durée du filtre en symboles
sps = 8;             % Échantillons par symbole (200 MHz / 8 = 25 MHz)

%% ==========================================
%% 2. GÉNÉRATION DU SIGNAL NUMÉRIQUE
%% ==========================================

% 2.1 Configuration des paramètres de modulation
switch upper(modulation_type)
    case 'BPSK'
        M = 2;  k = 1;   % 2 états, 1 bit par symbole
    case 'QPSK'
        M = 4;  k = 2;   % 4 états, 2 bits par symbole
    case '16QAM'
        M = 16; k = 4;   % 16 états, 4 bits par symbole
    otherwise
        error('Modulation non supportée. Utilisez BPSK, QPSK ou 16QAM.');
end

% 2.2 Calcul du nombre de bits nécessaire
% On ajoute 2*span pour compenser les effets de bord du filtrage
num_symbols_raw = ceil(target_length / sps) + 2*span;
num_bits_needed = num_symbols_raw * k;

% 2.3 Génération de bits aléatoires
fprintf('Génération de %d bits aléatoires...\n', num_bits_needed);
tx_bits = randi([0 1], num_bits_needed, 1);

% 2.4 Modulation en bande de base
fprintf('Modulation %s en cours...\n', modulation_type);
tx_symbols = qammod(tx_bits, M, 'InputType', 'bit', 'UnitAveragePower', true);

% Vérification de la présence de NaN ou Inf
if any(isnan(tx_symbols(:))) || any(isinf(tx_symbols(:)))
    error('ERREUR : La modulation a produit des valeurs NaN ou Inf !');
end

% 2.5 Conception du filtre RRC
fprintf('Conception du filtre RRC (rolloff=%.2f, span=%d, sps=%d)...\n', rolloff, span, sps);
h_rrc = rcosdesign(rolloff, span, sps, 'sqrt');

% 2.6 Filtrage de mise en forme
fprintf('Filtrage RRC en cours...\n');
tx_filtered = upfirdn(tx_symbols, h_rrc, sps);

% 2.7 Extraction de la partie centrale stable
% On évite les transitoires du début et de la fin du signal filtré
fprintf('Extraction de la partie centrale stable...\n');
center_idx = floor(length(tx_filtered) / 2);
half_len = floor(target_length / 2);

start_idx = max(1, center_idx - half_len + 1);
end_idx = min(length(tx_filtered), center_idx + half_len);

PAin = tx_filtered(start_idx:end_idx);

% 2.8 Vérification de la longueur finale
if length(PAin) ~= target_length
    fprintf('ATTENTION : Longueur ajustée de %d à %d échantillons\n', target_length, length(PAin));
    target_length = length(PAin);
end

%% ==========================================
%% 3. PRÉPARATION ET NORMALISATION DU SIGNAL
%% ==========================================

% 3.1 Calcul des statistiques du signal
signal_power = mean(abs(PAin).^2);
signal_peak = max(abs(PAin));
PAPR_linear = signal_peak^2 / signal_power;
PAPR_dB = 10 * log10(PAPR_linear);

fprintf('--- Statistiques du signal brut ---\n');
fprintf('Puissance moyenne : %.4f\n', signal_power);
fprintf('Amplitude crête   : %.4f\n', signal_peak);
fprintf('PAPR              : %.2f dB\n', PAPR_dB);

% 3.2 Normalisation CRITIQUE pour WebLab
% Le WebLab attend un signal avec amplitude crête normalisée
% Nous utilisons un facteur de sécurité de 0.9 pour éviter la saturation
safety_factor = 0.9;
PAin_normalized = PAin / signal_peak * safety_factor;

% 3.3 Vérification finale des valeurs
if any(isnan(PAin_normalized(:)))
    error('ERREUR : PAin contient des valeurs NaN après normalisation !');
end

if any(isinf(PAin_normalized(:)))
    error('ERREUR : PAin contient des valeurs Inf après normalisation !');
end

% 3.4 Formatage en vecteur colonne
PAin = double(PAin_normalized(:));

% Conversion en complexe si nécessaire
if isreal(PAin)
    PAin = complex(PAin, zeros(size(PAin)));
end

% 3.5 Calcul des statistiques finales
final_power = mean(abs(PAin).^2);
final_peak = max(abs(PAin));
final_PAPR_dB = 10 * log10(final_peak^2 / final_power);

fprintf('\n--- Statistiques après normalisation ---\n');
fprintf('Amplitude crête   : %.4f (limite sécurité = %.2f)\n', final_peak, safety_factor);
fprintf('Puissance moyenne : %.6f\n', final_power);
fprintf('PAPR final        : %.2f dB\n', final_PAPR_dB);

%% ==========================================
%% 4. CONFIGURATION DES VARIABLES POUR MAIN.M
%% ==========================================

% 4.1 Fréquence d'échantillonnage (requise par main.m)
Fs = fs_weblab;

% 4.2 Largeur de bande du signal
BW = fs_weblab / sps;

% 4.3 Paramètres ACPR pour le calcul de linéarité
% Ces valeurs correspondent au standard LTE 20MHz
ACPR.BW = 18e6;       % Largeur de bande du canal utile
ACPR.Offset = 20e6;   % Décalage du canal adjacent

% 4.4 Informations supplémentaires pour le débogage
signal_info.modulation = modulation_type;
signal_info.original_length = length(tx_filtered);
signal_info.final_length = length(PAin);
signal_info.sps = sps;
signal_info.rolloff = rolloff;
signal_info.span = span;
signal_info.PAPR_dB = final_PAPR_dB;
signal_info.peak_amplitude = final_peak;
signal_info.rms_amplitude = sqrt(final_power);

% 4.5 Vérification de la compatibilité avec main.m
fprintf('\n--- Vérification de compatibilité avec main.m ---\n');

% Vérification que PAin est complexe
if isreal(PAin)
    fprintf('ATTENTION : PAin est un signal réel, conversion en complexe...\n');
    PAin = complex(PAin, zeros(size(PAin)));
end

% Calcul du PAPR comme dans main.m
PAPRin = 10*log10(max(abs(PAin).^2) / mean(abs(PAin).^2));
fprintf('PAPR calculé (méthode main.m) : %.2f dB\n', PAPRin);

% Calcul de RMSin comme dans main.m (pour référence)
RMSin_calculated = -8.5 - PAPRin - 2;
fprintf('RMSin estimé (pour main.m) : %.2f dBm\n', RMSin_calculated);

%% ==========================================
%% 5. SAUVEGARDE DES DONNÉES
%% ==========================================

% 5.1 Nom du fichier (DOIT être 'PAinLTE20MHz.mat' pour main.m)
filename = 'PAinLTE20MHz.mat';

% 5.2 Liste des variables à sauvegarder
variables_to_save = {'PAin', 'Fs', 'BW', 'ACPR', 'tx_bits', 'signal_info'};

% 5.3 Sauvegarde avec vérification
fprintf('\n--- Sauvegarde des données ---\n');
fprintf('Fichier : %s\n', filename);

try
    % Sauvegarde des variables essentielles
    save(filename, variables_to_save{:});
    
    % Vérification que le fichier peut être chargé
    loaded_data = load(filename);
    
    % Vérification de chaque variable
    for i = 1:length(variables_to_save)
        var_name = variables_to_save{i};
        if ~isfield(loaded_data, var_name)
            error('Variable %s manquante dans le fichier sauvegardé', var_name);
        end
    end
    
    fprintf('✅ Fichier sauvegardé avec succès\n');
    
catch ME
    fprintf('❌ ERREUR lors de la sauvegarde : %s\n', ME.message);
    rethrow(ME);
end

%% ==========================================
%% 6. RÉSUMÉ FINAL ET INSTRUCTIONS
%% ==========================================

fprintf('\n==========================================\n');
fprintf('GÉNÉRATION TERMINÉE - RÉSUMÉ FINAL\n');
fprintf('==========================================\n');

fprintf('Signal généré : %s\n', modulation_type);
fprintf('Échantillons  : %d (à %d MHz)\n', length(PAin), Fs/1e6);
fprintf('Durée         : %.3f ms\n', length(PAin)/Fs*1000);
fprintf('Largeur bande : %.1f MHz\n', BW/1e6);
fprintf('Amplitude crête : %.4f\n', max(abs(PAin)));
fprintf('PAPR          : %.2f dB\n', PAPRin);

fprintf('\nVariables sauvegardées dans %s :\n', filename);
for i = 1:length(variables_to_save)
    fprintf('  • %s\n', variables_to_save{i});
end

fprintf('\n==========================================\n');
fprintf('INSTRUCTIONS POUR L''UTILISATION AVEC WEBLAB\n');
fprintf('==========================================\n');

fprintf('1. Assurez-vous que %s est dans le répertoire courant\n', filename);
fprintf('2. Exécutez main.m pour lancer la caractérisation PA\n');
fprintf('3. Les étapes suivantes seront :\n');
fprintf('   a) Mesure des caractéristiques AM-AM/AM-PM\n');
fprintf('   b) Calcul du ACPR et EVM\n');
fprintf('   c) Identification des coefficients DPD\n');
fprintf('   d) Application de la prédistorsion numérique\n');

% Vérification finale
fprintf('\n=== VÉRIFICATION FINALE ===\n');

% Vérification de l'absence de NaN/Inf
if any(isnan(PAin(:))) || any(isinf(PAin(:)))
    fprintf('❌ ERREUR : PAin contient des valeurs non finies !\n');
else
    fprintf('✅ PAin ne contient pas de NaN ou Inf\n');
end

% Vérification de la puissance
if final_power < 1e-6
    fprintf('⚠️  ATTENTION : Puissance très faible (%.6f)\n', final_power);
else
    fprintf('✅ Puissance dans une plage raisonnable\n');
end

% Vérification de l'amplitude crête
if final_peak > 1.0
    fprintf('❌ ERREUR : Amplitude crête > 1.0 (risque de saturation)\n');
elseif final_peak < 0.1
    fprintf('⚠️  ATTENTION : Amplitude crête très faible (%.4f)\n', final_peak);
else
    fprintf('✅ Amplitude crête correcte (%.4f)\n', final_peak);
end

fprintf('\n=== PRÊT POUR L''EXÉCUTION AVEC WEBLAB ===\n');