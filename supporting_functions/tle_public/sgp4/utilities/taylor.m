function weight = taylor(nx,slldb)
%
%INPUT:
%    nx = number of data samples
% slldb = sidelobe suppression in dB
%
%OUTPUT:
%     w = column vector of real weights
%
%Chris Wilson 20 August 2001
%
weight = ones(nx,1);                        %default weights
x = linspace(-(nx-1)/2,(nx-1)/2,nx)/nx;     %scaled positions of data samples
a = acosh(10^(abs(slldb)/20))/pi;           %scale parameter
nbar = round(2*a^2+0.5);                    %number of terms in sum
if nbar<2; return; end
sigmasq = nbar^2/(a^2+(nbar-0.5)^2);
%NOTE: vectorizing causes errors for extremely large nx!
f = zeros(1,nbar-1);
for m=1:(nbar-1)
    f(m) = (-1)^(m+1)/2;
    for n=1:(nbar-1)
        f(m)=f(m)*(1-m^2/sigmasq/(a^2+(n-0.5)^2));
        if m~=n
            f(m)=f(m)/(1-(m/n)^2);
        end
    end
    weight = weight +2*f(m)*cos(2*pi*m*x');
end
fsum = 1+2*sum(f);
weight = weight/fsum;
