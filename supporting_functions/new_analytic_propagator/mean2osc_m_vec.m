function [osc_orbital_elements,theta_osc,E_osc]=mean2osc_m_vec(x,param)

% input

%  x(1) = mean semimajor axis (kilometers)
%  x(2) = mean orbital eccentricity (non-dimensional)
%              (0 <= eccentricity < 1)
%  x(3) = mean orbital inclination (radians)
%              (0 <= inclination <= pi)
%  x(4) = mean right ascension of ascending node (radians)
%              (0 <= raan <= 2 pi)
%  x(5) = mean argument of perigee (radians)
%              (0 <= argument of perigee <= 2 pi)
%  x(6) = mean mean anomaly (radians)
%              (0 <= true anomaly <= 2 pi)

% convert from mean anomaly to true anomaly
e = x(:,2);
Epw = x(:,6);
Mo = Epw;
check_mag = 1:numel(Epw);
for l=1:1000
    Epw_check = Epw(check_mag);
    e_check = e(check_mag);
    DeltaEpw = - (Mo(check_mag) - Epw_check + e_check.*sin(Epw_check))./(-1+e_check.*cos(Epw_check));
    Epw(check_mag) = Epw_check+DeltaEpw;
    check_mag = check_mag(DeltaEpw > 10e-14);
    if isempty(check_mag) %if all elements have gone below threshold
        break
    end
end
sqrt_ep1_em1 = sqrt((1+e)./(1-e));
theta = 2*atan(sqrt_ep1_em1.*tan(Epw/2));

oeosc = mean2osc_vec([x(:,1:end-1),theta],param);
theta_osc = oeosc(:,6);
e_osc = oeosc(:,2);
%compute the osculating mean anomaly
E_osc = 2*atan(sqrt((1+e_osc)./(1-e_osc)).^(-1).*tan(theta_osc/2));
M_osc = E_osc-e_osc.*sin(E_osc);

osc_orbital_elements = [oeosc(:,1:5),M_osc];