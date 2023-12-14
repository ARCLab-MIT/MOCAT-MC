function [out_oe,error]=analytic_propagation(input_oe,param)

error=0;
% this code includes the solution for the mean elements
% as a function of time from Ref. 1.
%
% [Ref. 1] Martinusi, Vladimir, Lamberto Dell Elce, and Gaëtan Kerschen.
% "Analytic propagation of near-circular satellite orbits
% in the atmosphere of an oblate planet." Celestial Mechanics
% and Dynamical Astronomy 123, no. 1 (2015): 85-103.

re = param.req;
J2 = param.j2;
mu = param.mu;

if strcmpi(param.density_profile,'JB2008')

    % if statement ensure we get reasonable values for rho. Rho=0 break the
    % model below. Also altitude greater than 2300 km breaks the model.
    if (input_oe(1,1)-re)<2000 && (input_oe(1,1)-re)>200% changes from 200-1100
          rho_0 = lininterp2(param.alt(:,1),param.dens_times(1,:),param.dens_value,input_oe(1,1)-re, param.jd)*(1000)^3;
    else
        if (input_oe(1,1)-re)>1100
            rho_0 = lininterp1(param.dens_times(1,:),param.dens_value(end,:),param.jd)*(1000)^3;
        else
            rho_0 = lininterp1(param.dens_times(1,:),param.dens_value(1,:),param.jd)*(1000)^3;
        end
    end

elseif strcmpi(param.density_profile,'static')
    rho_0 = densityexp(input_oe(1,1)-re)*(1000)^3;
end
%constants
%C_0=param.C_0*rho_0;%1/2*Cd*s_ref/m*rho_0

%check these unit unit of 1/re and using 1/(6378.137 km) to covert to km
%units

% rho_reference for Bstar is 0.157 in units of kg/(m^2 * re).
% Using (1000^2 m^2/1 km^2) to convert A from m^2 to km^2, and dividing by 0.157 to scale
% according to newly computed density rho_0 (kg/km^3)
C_0 = max((abs(param.Bstar)/(1e6*0.157))*rho_0,1e-16); 

k2=mu*J2*re^2/2;

t=param.t;
t_0=param.t_0;

%initial conditions
a_0=input_oe(1,1);
e_0=input_oe(2,1);
inc_0=input_oe(3,1);
bigO_0=input_oe(4,1);
omega_0=input_oe(5,1);
Mo_0=input_oe(6,1);

c=cos(inc_0);
n_0=sqrt(mu)*a_0^(-3/2);

% p_0=a_0*(1-e_0^2);

% Cu=(3/4)*n_0*J2*(re/p_0)^2;

alpha_0=e_0/sqrt(a_0);

beta_0=(sqrt(3)/2)*e_0;

tan_coeff = max((tan(atan(beta_0)-beta_0*n_0*a_0*C_0*(t-t_0))),0);
a=(a_0/beta_0^2)*tan_coeff^2;

e=(1/(sqrt(3)/2))*tan_coeff;

Mo= (1/8)*(1/C_0)*(4/a+3*alpha_0^2*log(a/a_0))+...
    -(1/8)*(1/C_0)*(4/a_0+3*alpha_0^2*log(a_0/a_0))+...
    (3*k2*(3*c^2-1))/(16*mu)*(1/C_0)*((3*alpha_0^2/2)*1/a^2+4/(3*a^3))+...
    -(3*k2*(3*c^2-1))/(16*mu)*(1/C_0)*((3*alpha_0^2/2)*1/a_0^2+4/(3*a_0^3))+...
    Mo_0;

omega=(3*k2*(5*c^2-1))/(16*mu)*(1/C_0)*((5*alpha_0^2/2)*1/a^2+4/(3*a^3))+...
    -(3*k2*(5*c^2-1))/(16*mu)*(1/C_0)*((5*alpha_0^2/2)*1/a_0^2+4/(3*a_0^3))+omega_0;

bigO=-(3*k2*c)/(8*mu)*(1/C_0)*((5*alpha_0^2/2)*1/a^2+4/(3*a^3))+...
    (3*k2*c)/(8*mu)*(1/C_0)*((5*alpha_0^2/2)*1/a_0^2+4/(3*a_0^3))+bigO_0;

if isreal(inc_0) && isreal(bigO) && isreal(omega) && isreal(Mo)
    out_oe=[a;e;mod(inc_0, 2*pi);mod(bigO, 2*pi);mod(omega, 2*pi);mod(Mo, 2*pi)];
else
    out_oe=input_oe;
    error=1;
end
