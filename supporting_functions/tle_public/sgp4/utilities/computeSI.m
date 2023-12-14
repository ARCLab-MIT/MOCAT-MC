function [SI,Ibar]=computeSI(I,nspan)
%USAGE:  [SI,Ibar]=computeSI(I,nspan)
%
%Symmetric boxcar average with reflected end corrections
%
npts=length(I);
nmid=floor(nspan/2)+1; nupd=(nmid-1);
Ibar=zeros(size(I)); SI=Ibar;
for m=1:npts
    mspan=[max(1,m-nupd):min(m+nupd,npts)]; ntot=length(mspan);
    Ibar(m)=sum(I(mspan))/ntot;
    SI(m)=sqrt(sum(I(mspan).^2)/(ntot*Ibar(m)^2)-1);
end
return