function [S_f,f] = spectrum2(s_t,Te,typ)
%SPECTRUM Affiche la réponse fréquentielle d'un signal.
%   s_t : signal temporel
%   Te : période d'échantillonnage (ici T/F)
%	Type : 'linear', 'log' : scale for display (linear by default)

% x_t=zeros(1,length(s_t)+1);
% x_t(2:length(x_t))=s_t;
% modification MArs 2015: 2 types d'affichage (Linéaire ou Log/Décibel)

if nargin == 2
	typ = 'linear';
end

%[S_f,f]=freqz(s_t,1,8192,1/Te);
[S_f,f]=freqz(s_t,1,65536,1/Te);

figure;
subplot(211);
if strcmp(typ,'log')
	plot(f,20*log10(abs(S_f)));
    grid;
	ylabel('Gain [dB]');
else
	plot(f,abs(S_f).^2);
	ylabel('Gain [Echelle Lineaire]');
    grid;
end
title('Module de la transformée de Fourier');
xlabel('Fréquence [Hz]');
subplot(212);
plot(f,angle(S_f));
grid;
title('Phase de la transformée de Fourier');
xlabel('Fréquence [Hz]');
ylabel('Phase [rad]');

