function P=polyrec(Ps,MU)
%    P=polyrec(Ps,MU)
%    Recover original polynomial coefficients from call to [Ps,S,MU]=polyfit(x,y,N)
%    

N=length(Ps)-1;
%Coefficients are stored in decending order of n
for in=1:N+1;
    n=N-in+1;
    PsN(in)=Ps(in)/MU(2).^n;
end
P(1)=-PsN(1);     
SN=-1;
for in=2:N+1     %Current index
    m=N-in+1;    %Current power [N-1:-1:0]        
    P(in)=0;
    S=-1;
    for l=1:N-m+1
        S=S*SN;
        P(in)=P(in)+S*PsN(l)*C(N-l+1,m)*MU(1)^(N-l+1-m);
    end
    P=(-1)^(N-m)*P;
end  
return
function result=C(n,m)
if m==0
    result=1;
elseif m==1;
    result=n;
else
    result=factorial(n)/factorial(m)/factorial(n-m);
end
return