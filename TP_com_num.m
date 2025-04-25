%% Initialisation
clear all; close all; clc;

%% Paramètres généraux
%N = 16;           % Nombre de bits (Question 1)
N = 2048;
T = 1e-6;           % Durée symbole (1 µs) (Question 3)
F = 16;             % Facteur de surechantillonnage (Question 5)
alpha = 0.5;        % Roll-off pour SRRC (Question 6)
L = 8;              % Longueur du filtre (en symboles) (Question 6)
EbN0_dB = 2;       % Rapport signal/bruit (dB) (Question 11)


%% === Question 1 === (Avec rand et round)
bn = round(rand(1, N));          % Génération des bits
moyenne_bn = mean(bn);
variance_bn = var(bn);
disp('=== Question 1 ===');
disp(['Moyenne empirique : ', num2str(moyenne_bn), ' (théorie : 0.5)']);
disp(['Variance empirique : ', num2str(variance_bn), ' (théorie : 0.25)']);

%% === Question 2 ===
ak = 2 * bn - 1;                 % Mapping 2-PAM

moyenne_ak = mean(ak);           %moyenne
variance_ak = var(ak);           %variance
puissance_ak = mean(ak.^2);      %puissance moyenne
disp('=== Question 2 ===');
disp(['Moyenne empirique de ak : ', num2str(moyenne_ak), ' (théorie : 0)']);
disp(['Variance empirique de ak : ', num2str(variance_ak), ' (théorie : 1)']);
disp(['Puissance moyenne de ak : ', num2str(puissance_ak), ' (théorie : 1)']);

%% === Question 3 ===
t_a = (0:N-1) * T;
figure;
plot(t_a, ak, '*');
xlabel('Temps (s)'); ylabel('Amplitude (V)');
title('Question 3 : Symboles 2-PAM dans le temps');
grid on;

%% === Question 4 ===
figure;
plot(real(ak), imag(ak), 'x');
xlabel('Partie Réelle'); ylabel('Partie Imaginaire');
title('Question 4 : Constellation 2-PAM');
axis([-1.5 1.5 -0.1 0.1]);
grid on;

%% === Question 5 ===
st = zeros(1, N * F);
st(1:F:length(ak)*F) = F * ak;
t_s = (0:length(st)-1) * T / F;
figure;
plot(t_s, st);
xlabel('Temps (s)'); ylabel('Amplitude (V)');
title('Question 5 : Signal surechantillonné st');
grid on;

%% === Question 6 ===
t_filtre = 0:T/F:L*T - T/F;
h_e = gen_filters3('srrc', t_filtre, T, F, L, alpha);


% Réponse impulsionnelle
figure;
subplot(2,1,1);
plot(t_filtre, h_e);
title('Question 6 : Réponse impulsionnelle du filtre SRRC');
xlabel('Temps (s)'); ylabel('Amplitude');

% Réponse fréquentielle
subplot(2,1,2);
spectrum2(h_e, T/F, 'lin');
title('Réponse fréquentielle');

%% === Question 7 ===
% xt = conv(st, h_e);
% t_xt = (0:length(xt)-1) * T / F;
% figure;
% subplot(2,1,1);
% plot(t_s, st); title('Signal avant filtrage');
% xlabel('Temps (s)'); ylabel('Amplitude');
% subplot(2,1,2);
% plot(t_xt, xt); title('Signal après filtrage');
% xlabel('Temps (s)'); ylabel('Amplitude');

% Filtrage par convolution (xt)
xt = conv(st, h_e);      % 'full' pour garder toute la convolution
L_xt = length(xt);              % Longueur du signal filtré

% Vecteurs temps
t_st = (0:length(st)-1) * (T/F);               % Temps pour st
t_xt = (0:L_xt-1) * (T/F);                     % Temps pour xt

% Affichage (subplot)
figure;
subplot(2,1,1);
plot(t_st, st);
grid;
axis([0 L_xt*T/F -1.5 1.5]);    % Même échelle de temps pour les deux signaux
title('Signal avant filtre de mise en forme (st)');
xlabel('Temps [s]');
ylabel('Amplitude [V]');

subplot(2,1,2);
plot(t_xt, xt, 'r');
grid;
axis([0 L_xt*T/F -1.5 1.5]);
title('Signal après filtre de mise en forme SRRC (xt)');
xlabel('Temps [s]');
ylabel('Amplitude [V]');

%% === Question 8 ===
P_xt = mean(xt.^2);
disp('=== Question 8 ===');
disp(['Puissance de xt : ', num2str(P_xt), ' (théorie ≈ 1)']);

%% === Question 9 ===
figure;
spectrum2(xt, T/F, 'log');
title('Question 9 : DSP du signal émis (échelle log)');

%% === Question 10  ===
ak_non_centre = 2 * bn; % 0 → 0, 1 → 2
st_non_centre = zeros(1, N*F);
st_non_centre(1:F:end) = F * ak_non_centre;
xt_non_centre = conv(st_non_centre, h_e);
figure;
plot(t_xt, xt_non_centre);
title('Question 10 : Signal avec symboles non centrés');
%% === Question 11 === 
EbN0_lin = 10^(EbN0_dB/10);
%sigma_n = sqrt((F* P_xt) / (2 * EbN0_lin));
sigma_n = sqrt(F / (2 * EbN0_lin));

nt = sigma_n * randn(1, length(xt));
rt = xt + nt;
figure;
plot(t_xt, rt);
title('Question 11 : Signal reçu bruité');
%% === Question 13  ===

p_t= conv(h_e,fliplr(h_e));
%t2= T/F :2*(L*T)-T/F;
t2 = (0:length(h_e)-1) * T / F;
t3 = (0:length(p_t)-1) * T / F;

% affichage de h_e
figure;
subplot(2,1,1);
plot(t2, h_e);
title('Question 13 : Réponse impulsionnelle du filtre SRRC');
xlabel('Temps (s)'); ylabel('Amplitude');

%affiche de p_t
subplot(2,1,2);
plot(t3, p_t);
title('Réponse impulsionelle du filtre h_t convolué a h_e');
xlabel('Temps (s)'); ylabel('Amplitude');


%% === Question 14 ===

% 1. Vérification visuelle du critère de Nyquist
indices_nyquist = 1:F:length(p_t);
figure;
stem(indices_nyquist, p_t(indices_nyquist));
title('Échantillons de pt aux instants multiples de T');
xlabel('Index'); ylabel('Amplitude');
grid on;

% 2. Filtrage du signal reçu
y_t = conv(rt, fliplr(h_e));

% 3. Diagramme de l’œil
taille_visu = 2;
begin_offset = length(h_e);
end_offset = 2 * length(h_e);
eyepattern(y_t, T, F, taille_visu, begin_offset, end_offset);


%% === Question 18 ===
t0 = L * F;
y_k = zeros(1,N);
y_k(1:1:N) = y_t(t0:F:(N-1)*F + t0);


%% === Question 19 ===

figure;
plot(real(y_k), imag(y_k), 'x');  % Imaginaire nul ici
title('Constellation des symboles reçus y_k');
xlabel('Partie réelle'); ylabel('Partie imaginaire');
grid on;

b_decide = y_k>0;
erreurs = sum(xor(bn(1 :length(bn)),b_decide)) ;
teb=erreurs/N ;





