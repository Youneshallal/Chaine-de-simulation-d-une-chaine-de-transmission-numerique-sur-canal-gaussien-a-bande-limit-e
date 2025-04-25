% --- getBER.m 
% --- 
% Description 
%    Returns uncoded Bit Error Rate (BER) for several modulation scheme and different constellations of size M  for a given Eb / No vector.
% --- 
% Syntax 
%					  [ber]	  = getBER(ebNo,'mod',M)
%									  --- ebNo		: Vector of desired Energy ber bit [Eb] versus noise spectral density [No] (in linear case)
%									  --- 'mod'		: target modulation [PAM, PSK, QAM]
%									  --- M			: Constellation size : 2^(bit per Syntaxmbol)
%    
% --- 
% v 1.0 - Robin Gerzaguet


function [ber]   = getBER(ebNo,modulation,M)
	% --- Checking that modulation is a power of 2
	if floor(log2(M)) ~= log2(M)
		error('simRF:getBER','Modulation size must be a power of 2');
	end
	% --- Switching Modulation type 
	switch modulation 
		case 'PAM'
			% --- Bit Error rate for PAM modulations
			ber	  = 2*(M-1) / (M*log2(M))* Q(sqrt( 6 * log2(M) / (M^2-1) * ebNo));
		case 'PSK'
            if M >=4 
			ber	  = 2 / log2(M) * Q( sqrt(2 * ebNo * log2(M) * (sin(pi/M))^2));
            end
            if M==2 
                ber = Q( sqrt(2 * ebNo));
            end
		case 'QAM'
			% --- Bit Error Rate for QAM modulations
			ber	  = 4 * ( 1 - 1 / sqrt(M)) / log2(M) * Q(sqrt( 2*ebNo * 3 * log2(M) / (2*(M-1)  )));
		otherwise
			error('simRF:getBER','Unkown modulation type');
	end
end
% --- Define Q function 
function out = Q(in)
	out = 1/2*erfc(in / sqrt(2));
end
