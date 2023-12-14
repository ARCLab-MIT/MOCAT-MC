function [falias,m]=alias(f,fNyq)
%
nf=length(f);
falias=zeros(1,nf);
for n=1:nf;
    err=1.e12; m=0;
    while err>fNyq/2
        err=abs(f(n)-m*fNyq);
        m=m+1;
    end
    falias(n)=f(n)-(m-1)*fNyq;
end
return