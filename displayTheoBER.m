% --- displayTheoBER.m 
% --- 
% Description 
%    Display Theoretical Bit Error Rate for classical modulations (PAM, QAM, QPSK) and different modulation size 
% --- 
% Syntax 
%    displayTheoBER 
% --- 
% v 1.0 - Robin Gerzaguet.

% --- Clear all stuff 
clear;

% --- modulations
modulations		= {'PAM','QAM','PSK'};
sizeConst		= [2 4 8 16 32];
ebNo			= 1:0.5:25;
ebLin			= 10.^(ebNo/10);

% --- Iterative display 
for iMod = 1 : 1 : size(modulations,2)
	modCurr			  = modulations{iMod};
	figure('name',['BER for ' modCurr]);
	strL			  = cell(1,length(sizeConst));
	col				  = hsv(length(sizeConst));
	% --- Iterative performance on constellation size
	for iConst = 1 : 1 : length(sizeConst)
		constCurr	  = sizeConst(iConst);
		ber			  = getBER(ebLin,modCurr,constCurr);
		semilogy(ebNo,ber,'lineWidth',3,'LineStyle','-','color',col(iConst,:));
    	hold on
		strL{iConst} = [modCurr ' - ' num2str(constCurr)];
	end 
	grid; 
    grid minor % minor;
	xlabel('Eb / No [dB]');
	ylabel('Bit Error Rate');
	legend(strL);
	title(['Bit Error Rate for Modulation : ' modCurr]);
	ylim([10^-10 1]);
end


