function [mat_frag]=func_create_tlesv2_vec(ep, r_parent,v_parent,class_parent, fragments, max_frag,mu,req,maxID)

% Create new sat objects from fragmentation information

% ep:               epoch
% sat:              original satellite that is fragmented
% fragments:        Nx8 data on fragments: 
%           [diam,  Area,  AMR,   m,   total_dv, dv_X,  dv_Y, dv_Z];


% OUTPUT:
% mat_frag = [a,ecco,inclo,nodeo,argpo,mo,bstar,mass,radius,...
%             errors,controlled,a_desired,missionlife,constel,...
%             date_created,launch_date,r,v,frag_objectclass,ID_frag]; [num_sat x 24]

% PARAMETERS TO EXPELL FROM FUNCTION: 
% max_frag, mu, req


% %max_frag=300;
% max_frag=1200;

% mu = 398600.5;
% re = 6378.137;
% 
% slct=2;

% M=12;
% D=1;
% Y=2022;
% H=3;
% MI=4;
% S=30;

% req = 6378.14;

% needs to be current time 
% ep_date = datetime(Y,M,D,H,MI,S)+days(ep/(24*60));
% mtle=month(ep_date);
% dtle=day(ep_date);
% ytle=year(ep_date);
% tlehr=hour(ep_date);
% tlemin=minute(ep_date);
% tlesec=second(ep_date);

r0 = r_parent;
v0 = v_parent;

if size(fragments,1) > max_frag
    warning('number of fragments %i exceeds number of max_frag specified %i',...
        size(fragments,1), max_frag);
end
n_frag = min([size(fragments,1) max_frag]);

% sort by mass to minimize mass conservation issues when limited by n_frag
[~,si] = sort(fragments(:,4),'descend'); % sort by mass (large ones first)
fragments = fragments(si(1:n_frag),:);

v = [fragments(:,6)+v0(1),fragments(:,7)+v0(2),fragments(:,8)+v0(3)]; % add dV
r = [r0(1)*ones(n_frag,1),r0(2)*ones(n_frag,1),r0(3)*ones(n_frag,1)];

[~,a,ecc,incl,omega,argp,~,m,~,~,~] = rv2coe_vec(r,v,mu);

idx_a = find(a>0); %if in an elliptical orbit (reject hyperbolic fragments)
num_a = numel(idx_a);

a = a(idx_a)/req;
ecco = ecc(idx_a);
inclo = incl(idx_a);
nodeo = omega(idx_a);
argpo = argp(idx_a);
mo = m(idx_a);

% Bstar parameter, https://en.wikipedia.org/wiki/BSTAR
rho_0 = 0.157; %rho_0 in units of kg/(m^2*Re)
A_M = fragments(idx_a,3); % Area to mass ratio in units of m^2/kg
bstar = (0.5*2.2*rho_0)*A_M; % Bstar in units of 1/re

mass = fragments(idx_a,4);
radius = fragments(idx_a,1)/2; 

errors = zeros(num_a,1);

controlled = zeros(num_a,1);
a_desired = nan(num_a,1);
missionlife = nan(num_a,1);

constel = zeros(num_a,1);

date_created = ep*ones(num_a,1);
launch_date = nan(num_a,1);

frag_objectclass = ones(num_a,1)*filter_objclass_fragments_int(class_parent); %Assign object class to fragments according to parent particle

ID_frag = linspace(maxID+1,maxID+num_a,num_a)'; %generate ID for new fragments

mat_frag = [a,ecco,inclo,nodeo,argpo,mo,bstar,mass,radius,...
            errors,controlled,a_desired,missionlife,constel,...
            date_created,launch_date,r(idx_a,:),v(idx_a,:),frag_objectclass,ID_frag];