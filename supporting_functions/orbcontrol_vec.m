function [mat_sat_out,deorbit] = orbcontrol_vec(mat_sat_in,tsince,time0,orbtol,PMD,DAY2MIN,YEAR2DAY,param)
%mat_sat_in = [a,ecco,inclo,nodeo,argpo,mo,controlled,a_desired,missionlife,launched,r,v]
%mat_sat_out = [a,controlled,r,v];

current_time = time0+days(tsince/DAY2MIN);

controlled = mat_sat_in(:,7);
a_desired = mat_sat_in(:,8);
missionlife = mat_sat_in(:,9);
launched = mat_sat_in(:,10);

a_out = mat_sat_in(:,1);
r_out = mat_sat_in(:,11:13);
v_out = mat_sat_in(:,14:16);

is_controlled = find(controlled==1); %find controlled satellites
a_current = a_out(is_controlled); %semi-major axis of controlled satellites
find_control = is_controlled(abs(a_current-a_desired(is_controlled))>(orbtol/param.req)); %identify controlled satellites beyond tolerance
a_out(find_control) = a_desired(find_control); %reset semi-major axis of controlled satellites beyond tolerance

osc_oe = mean2osc_m_vec([a_out(find_control,1)*param.req,mat_sat_in(find_control,2:6)],param); %compute new osculating elements
[r_out(find_control,:),v_out(find_control,:)]=oe2rv_vec(osc_oe,[],param); %reset position and velocity

%Satellites past their mission life
find_life = is_controlled((juliandate(current_time) - launched(is_controlled))/YEAR2DAY > missionlife(is_controlled)); % compare in years
rand_life = rand(numel(find_life),1); %generate random number for each satellite beyond life
check_PMD = PMD<rand_life; %check if PMD is fulfilled
controlled(find_life(check_PMD)) = 0; % from active becomes inactive 
deorbit = find_life(~check_PMD); %deorbited satellites

mat_sat_out = [a_out,controlled,r_out,v_out];