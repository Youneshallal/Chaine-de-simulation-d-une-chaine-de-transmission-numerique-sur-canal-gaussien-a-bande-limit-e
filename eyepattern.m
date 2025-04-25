function [] = eyepattern(y_t,T,F,taille_visu,begin_offset,end_offset)
%EYEPATTERN Affiche le diagramme de l'oeil.
%
% y_t : signal reçu (suréchantillonné).
% T : durée d'un symbole [s].
% F : facteur de suréchantillonnage.
% taille_visu : nombre de symboles à visualiser.
% begin_offset : nombre d'échantillons à écarter en début de signal.
% end_offset : nombre d'échantillons à écarter en fin de signal.

% Création du vecteur temps
tvisu=0:T/F:taille_visu*T-T/F;

% Troncature du signal à visualiser
y_t=y_t(begin_offset:end-end_offset);

% Visualisation sur un nombre entier de trajectoires
nb_trajectoires=floor(length(y_t)/length(tvisu));

% Remise en forme du signal pour diagramme de l'?il
y_t_mat=reshape(y_t(1:nb_trajectoires*taille_visu*F),taille_visu*F,nb_trajectoires);

% Affichage du diagramme de l'oeil
figure;
plot(tvisu,y_t_mat);
title('Diagramme de l''oeil du signal recu');
xlabel('Temps [s]');
ylabel('Amplitude [V]');
end

