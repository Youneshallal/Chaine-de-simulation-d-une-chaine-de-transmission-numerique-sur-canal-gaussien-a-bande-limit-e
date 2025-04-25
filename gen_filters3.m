function y = gen_filters3(type,t,T,F,L,alpha)
%GEN_FILTERS Retourne la réponse impulsionnelle de différents filtres de
%mise en forme.
%
% type : 'nrz' pour une porte de durée T
%        'rz' pour une porte de durée T/2
%        'srrc' pour un filtre en racine de cosinus surélevé
% t : vecteur temps [s]
% T : durée symbole [s]
% F : facteur de suréchantillonnage
% L : durée de la réponse impulsionnelle du filtre [T*s]
% alpha : coefficient d'excès de bande (pour filtre srrc)
%
% Par Damien ROQUE <damien.roque@grenoble-inp.fr>
% modification L. Ros Mars 2015: normalisation des filtres en 1/F.  
% modification L. Ros Mars 2017: sans utilisation de la fonction sinc(.)
% qui nécessite une bibliothèque spéciale

switch type
    case{'nrz'}
        % Création de l'impulsion rectangulaire NRZ (non-centrée)
        %y=1/sqrt(F)*[ones(1,F) zeros(1,length(t)-F)];
        y=(1/F)*[ones(1,F) zeros(1,length(t)-F)];
    case{'rz'}
        if mod(F,2)
            error('Le facteur de suréchantillonnage doit être un multiple de deux');
        end
        %y=1/sqrt(0.5*F)*[ones(1,F/2) zeros(1,length(t)-F/2)];
        y=(sqrt(2)/F)*[ones(1,F/2) zeros(1,length(t)-F/2)];
    case{'srrc'}
        % Décalage pour obtenir une réponse causale
        t=t-T*L/2;
        % Construction du vecteur de sortie dans le cas général (abs(t) != T/(4*alpha))
        %y=(1/sqrt(F))*(1./(1-(4*alpha.*t/T).^2)).*((1-alpha)*sinc((t/T)*(1-alpha))+(4*alpha/pi).*cos((pi*t/T)*(1+alpha)));
        y=(1/F)*(1./(1-(4*alpha.*t/T).^2)).*((1-alpha)*Sincard((t/T)*(1-alpha))+(4*alpha/pi).*cos((pi*t/T)*(1+alpha)));
       
        % Correction des éventuelles valeurs non-définies (t=0 et abs(t) = T/(4*alpha))
        for indice=1:length(t)
            %if t(indice)==0
            if abs(t(indice)) < 10^(-10)
               % y(indice)=(1/sqrt(F))*(1-alpha+4*(alpha/pi));
               y(indice)=(1/F)*(1-alpha+4*(alpha/pi));
            end
            %if abs(t(indice))==T/(4*alpha)
            if abs(abs(t(indice)) - T/(4*alpha)) < 10^(-10)
                %y(indice)=(alpha/sqrt(2*F))*((1+2/pi)*sin(pi/(4*alpha)) + (1-2/pi)*cos(pi/(4*alpha)));
             y(indice)=(alpha/(sqrt(2)*F))*((1+2/pi)*sin(pi/(4*alpha)) + (1-2/pi)*cos(pi/(4*alpha)));
            end
        end
end

end

% --- Define sincard function ------ 
function out = Sincard(in)
if(in==0)
    Sincard(in)=1;
else
	out = sin(pi*in)./(pi*in);
end
end