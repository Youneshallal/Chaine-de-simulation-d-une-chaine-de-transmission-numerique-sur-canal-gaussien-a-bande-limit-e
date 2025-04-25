% --- Define sincard function ------ 
function out = Sincard(in)
if(in==0)
    Sincard(in)=1;
else
	out = sin(pi*in)./(pi*in);
end
end