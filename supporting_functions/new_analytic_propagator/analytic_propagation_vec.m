function [out_oe,errors]=analytic_propagation_vec(input_oe,param)

errors=zeros(size(input_oe,1),1);
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

% if statement ensure we get reasonable values for rho. Rho=0 break the
% model below. Also altitude greater than 2300 km breaks the model.
a_0 = input_oe(:,1);
a_minus_re = a_0-re;

if strcmpi(param.density_profile,'JB2008')
    
    rho_0 = zeros(length(a_0),1);
    check_above = a_minus_re>param.alt(end,1);
    check_below = a_minus_re<param.alt(1,1);
    check_in_range = ~check_above & ~check_below;
    rho_0(check_in_range) = lininterp2_vec_v2(param.alt(:,1),param.dens_times(1,:),param.dens_value,a_minus_re(check_in_range), param.jd)*1e9;
    rho_0(check_above) = lininterp1_vec(param.dens_times(1,:),param.dens_value(end,:),param.jd)*1e9;
    rho_0(check_below) = lininterp1_vec(param.dens_times(1,:),param.dens_value(1,:),param.jd)*1e9;
 
elseif strcmpi(param.density_profile,'static')
    rho_0 = densityexp_vec(a_minus_re)*1e9;
end

% rho_reference for Bstar is 0.157 in units of kg/(m^2 * re).
% Using (1000^2 m^2/1 km^2) to convert A from m^2 to km^2, and dividing by 0.157 to scale
% according to newly computed density rho_0 (kg/km^3)% 1/(6378.137 km)*(1000^2 m/1 km) to convert to kg/km^3

C_0 = max((param.Bstar/(1e6*0.157)).*rho_0,1e-20); %temporary lower limit on C_0 to avoid singularity

k2_over_mu = J2*re^2/2; %k2 = mu*J2*re^2/2

t = param.t;
t_0 = param.t_0;

%Initial conditions
e_0 = input_oe(:,2);
inc_0 = input_oe(:,3);
bigO_0 = input_oe(:,4);
omega_0 = input_oe(:,5);
Mo_0 = input_oe(:,6);

c = cos(inc_0);
c_sq = c.^2;

n_0 = sqrt(mu)*a_0.^(-3/2);

alpha0_sq = (e_0./sqrt(a_0)).^2;

beta_0 = (sqrt(3)/2)*e_0;

tan_atan_beta0 = max(tan(atan(beta_0)-beta_0.*n_0.*a_0.*C_0*(t-t_0)),0); %place lower limit on eccentricity and semi-major axis reduction
a = (a_0./beta_0.^2).*tan_atan_beta0.^2;
e = (2/sqrt(3))*tan_atan_beta0; 

check_beta = beta_0==0;
if any(check_beta) %avoid a=NaN when beta_0=0; no numerical issue when beta_0<<1
    a0_beta = a_0(check_beta);
    a(check_beta) = a0_beta.*(1-C_0(check_beta).*n_0(check_beta).*a0_beta*(t-t_0)); 
end

%Compute some variables to avoid repetition of operations
a_sq = a.^2;
four_thirds_over_a_cb = 4/3./(a_sq.*a);
a0_sq = a_0.^2;
four_thirds_over_a0_cb = 4/3./(a0_sq.*a_0);
alpha0sq_over_asq = alpha0_sq./a_sq;
alpha0sq_over_a0sq = alpha0_sq./a0_sq;

Mo = (0.5./a - 0.5./a_0 + 3/8*alpha0_sq.*log(a./a_0))./C_0 + 3*k2_over_mu/16*(3*c_sq-1).*(1.5*(alpha0sq_over_asq - alpha0sq_over_a0sq) + four_thirds_over_a_cb - four_thirds_over_a0_cb)./C_0 + Mo_0;

five_a0sq_over2_tau2_plus_4thirds_over_tau3_over_C0 = (2.5*(alpha0sq_over_asq - alpha0sq_over_a0sq) + four_thirds_over_a_cb - four_thirds_over_a0_cb)./C_0;

omega = 3*k2_over_mu/16*(5*c_sq-1).*five_a0sq_over2_tau2_plus_4thirds_over_tau3_over_C0 + omega_0;

bigO = -3*k2_over_mu/8*c.*five_a0sq_over2_tau2_plus_4thirds_over_tau3_over_C0 + bigO_0;

out_oe = [a,e,mod(inc_0, 2*pi),mod(bigO, 2*pi),mod(omega, 2*pi),mod(Mo, 2*pi)];

not_real = ~isreal(inc_0) | ~isreal(bigO) | ~isreal(omega) | ~isreal(Mo);
errors(not_real) = 1;

out_oe(not_real,:) = input_oe(not_real,:);