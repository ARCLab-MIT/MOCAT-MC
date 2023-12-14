function IdB=dB20(A)
%     Compute 20*log10(A)  
%USAGE: IdB=dB20(A)
%INPUTS:
%    A = input amplitude vector  
%OUTPUTS:
%    dB if A>0 else NaN
%

% Chuck Rino
% Rino Consulting
% July 2010
%
IdB=NaN*ones(size(A));
iOK=find(abs(A)>0);
IdB(iOK)=20*log10(abs(A(iOK)));
return