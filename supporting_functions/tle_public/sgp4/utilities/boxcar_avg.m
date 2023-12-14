function   vbar=boxcar_avg(v,nspan);
%  M point running average
npts=length(v);
nmid=floor(nspan/2)+1; nupd=(nmid-1);
vbar=zeros(size(v));
for m=1:npts
    mspan=[max(1,m-nupd):min(m+nupd,npts)]; ntot=length(mspan);
    vbar(m)=sum(v(mspan))/ntot;
end
return