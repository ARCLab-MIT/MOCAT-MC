function [mat_sat_out]=prop_mit_vec(mat_sat_in,t,param)
%mat_sat_in = [a,ecco,inclo,nodeo,argpo,mo,Bstar,controlled]
%mat_sat_out = [a,ecco,inclo,nodeo,argpo,mo,errors,r_eci,v_eci];

req = param.req;

%set parameters
% param.F107=sat.F107;%current F107
% param.Ap=sat.Ap;%current Ap
param.t = t;% in units of seconds
param.t_0 = 0;% in units of seconds

n_sat = size(mat_sat_in,1);

in_mean_oe = [req*mat_sat_in(:,1),mat_sat_in(:,2:6)];

Bstar = abs(mat_sat_in(:,7));
Bstar(Bstar<1e-12) = 9.7071e-05;

%Preallocate output
out_mean_oe = zeros(n_sat,6);
errors = zeros(n_sat,1);

% check if a de-orbit has already happened or if controlled
idx_notdecay = in_mean_oe(:,1).*(1-in_mean_oe(:,2))>req+150;
idx_controlled = mat_sat_in(:,8)==1;
idx_propagate = idx_notdecay & ~idx_controlled; %objects that haven't decayed and are not controlled
if any(idx_propagate)
    param.Bstar = Bstar(idx_propagate);
    [out_mean_oe(idx_propagate,:),errors(idx_propagate)] = analytic_propagation_vec(in_mean_oe(idx_propagate,:),param);
end
out_mean_oe(~idx_propagate,:) = in_mean_oe(~idx_propagate,:); %decayed or controlled, assign input value
out_mean_oe(idx_controlled,6) = out_mean_oe(idx_controlled,6) + sqrt(param.mu./out_mean_oe(idx_controlled,1).^3)*t; %update mean anomaly of controlled objects, all other orbital elements remain the same

check_alt_ecc = out_mean_oe(:,1).*(1-out_mean_oe(:,2))>req+150 & out_mean_oe(:,2)<1; %check if decayed or hyperbolic
errors(~check_alt_ecc) = 1;

%%%%Issue with controlled objects decaying too quickly
% % check if a de-orbit has already happened and indicate error
% idx_propagate = in_mean_oe(:,1).*(1-in_mean_oe(:,2))>req+150;
% if any(idx_propagate)
%     param.Bstar = Bstar(idx_propagate);
%     [out_mean_oe(idx_propagate,:),errors(idx_propagate)] = analytic_propagation_vec(in_mean_oe(idx_propagate,:),param);
% end
% out_mean_oe(~idx_propagate,:) = in_mean_oe(~idx_propagate,:);
% errors(~idx_propagate) = 1;
% 
% check_alt_ecc = out_mean_oe(:,1).*(1-out_mean_oe(:,2))>req+150 & out_mean_oe(:,2)<1;
% errors(~check_alt_ecc) = 1;

%Mean to osculating orbital elements
osc_oe = zeros(n_sat,6);
[osc_oe(check_alt_ecc,:),~,E_osc] = mean2osc_m_vec(out_mean_oe(check_alt_ecc,:),param); %from Matlabfile exchange
out_mean_oe(:,1) = out_mean_oe(:,1)/req;
osc_oe(~check_alt_ecc,:) = out_mean_oe(~check_alt_ecc,:);

%Osculating orbital elements to state vector
r_eci = zeros(n_sat,3);
v_eci = zeros(n_sat,3);
[r_eci(check_alt_ecc,:),v_eci(check_alt_ecc,:)] = oe2rv_vec(osc_oe(check_alt_ecc,:),E_osc,param);

mat_sat_out = [out_mean_oe,errors,r_eci,v_eci];