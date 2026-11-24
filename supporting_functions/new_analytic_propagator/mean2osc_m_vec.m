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
%              (0 <= mean anomaly <= 2 pi)

% mean2osc_vec takes and returns the MEAN anomaly in column 6
% (delm solves Kepler's equation on its input and converts its output
% back to mean anomaly), so pass the mean anomaly straight through.
oeosc = mean2osc_vec(x,param);
M_osc = oeosc(:,6);
e_osc = oeosc(:,2);

% osculating eccentric anomaly from Kepler's equation
E_osc = M_osc;
check_mag = 1:numel(E_osc);
for l=1:1000
    E_check = E_osc(check_mag);
    e_check = e_osc(check_mag);
    DeltaE = - (M_osc(check_mag) - E_check + e_check.*sin(E_check))./(-1+e_check.*cos(E_check));
    E_osc(check_mag) = E_check+DeltaE;
    check_mag = check_mag(abs(DeltaE) > 1e-13);
    if isempty(check_mag) %if all elements have gone below threshold
        break
    end
end

% osculating true anomaly
theta_osc = 2*atan(sqrt((1+e_osc)./(1-e_osc)).*tan(E_osc/2));

osc_orbital_elements = [oeosc(:,1:5),M_osc];
