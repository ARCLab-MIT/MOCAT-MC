function [popSSEM] = Fast_MC2SSEM_population(sats_info,param)

% Divides the output of the supercloud runs by species and altitude bin. 
% The structure param contains the parameters that need to be set.
%---
% Authors: Peng Mun Siew, MIT 11/4/2022
%---

% species: [S,D,N,Su,B,U];
ind = param.species~=0;

idx_active = sats_info{1}==1 & sats_info{3}==1;
idx_derelict = sats_info{1}==1 & sats_info{3}~=1;
idx_debris = sats_info{1}==3 | sats_info{1}==4 | sats_info{1}>=7;
idx_rbodies = sats_info{1}==5 | sats_info{1}==6;

obj_alt = sats_info{2}* param.re - param.re;

% Separate sat alt based on classes
alt_active = obj_alt(idx_active);
alt_derelict = obj_alt(idx_derelict);
alt_debris = obj_alt(idx_debris);
alt_rbodies = obj_alt(idx_rbodies);

[s_count,~] = histcounts(alt_active,param.R02);
[d_count,~] = histcounts(alt_derelict,param.R02);
[n_count,~] = histcounts(alt_debris,param.R02);
[b_count,~] = histcounts(alt_rbodies,param.R02);
su_count = zeros(1,param.N_shell);
u_count = zeros(1,param.N_shell);

popSSEM = [s_count' d_count' n_count' su_count' b_count' u_count' ];
popSSEM = popSSEM(:,ind);

end