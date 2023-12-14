function simturb=turbsim2(rootSDF)
%        simturb=turbsim2(rootSDF)
%   
% Generate realization of field with SDF rootSDF.^2
%
[n1,n2]=size(rootSDF);
xi=zeros(n1,n2)+1i*zeros(n1,n2);
xi(1:n1/2+1,1:n2/2+1)  =((randn(n1/2+1,n2/2+1)+1i*randn(n1/2+1,n2/2+1)))/sqrt(2);
xi(n1/2+2:n1,2:n2/2)   =((randn(n1/2-1,n2/2-1)+1i*randn(n1/2-1,n2/2-1)))/sqrt(2);
xi(1,1)=0; 
xi(1,n2/2+1)     =real(xi(1,n2/2+1))/2;
xi(n1/2+1,1)     =real(xi(n1/2+1,1))/2;
xi(n1/2+1,n2/2+1)=real(xi(n1/2+1,n2/2+1))/2;

xi(n1/2+2:n1,1)     =conj(xi(n1/2:-1:2,1));
xi(n1/2+2:n1,n2/2+1)=conj(xi(n1/2:-1:2,n2/2+1));
xi(1,     n2/2+2:n2)=conj(xi(1,     n2/2:-1:2));
xi(n1/2+1,n2/2+2:n2)=conj(xi(n1/2+1,n2/2:-1:2));

xi(2:n1/2,n2/2+2:n2)   =conj(xi(n1:-1:n1/2+2,n2/2:-1:2));
xi(n1/2+2:n1,n2/2+2:n2)=conj(xi(n1/2:-1:2,   n2/2:-1:2));

simturb=real(fft2(fftshift(rootSDF.*xi)));  %NOTE: This scaling preserves variance
return