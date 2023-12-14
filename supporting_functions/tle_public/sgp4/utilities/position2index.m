function [i1,i2]=position2index(ipos,n1)
%
%   invert ipos=(i2-1)*n1+i1
%
nn=length(ipos);
i1=zeros(nn,1);
i2=zeros(nn,1);
for n=1:nn
  i1(n)=mod(ipos(n),n1);
  iZero= i1==0; 
  i1(iZero)=n1;
  i2(n)=(ipos(n)-i1(n))/n1+1;
end
return