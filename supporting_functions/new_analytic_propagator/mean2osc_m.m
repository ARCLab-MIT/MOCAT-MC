function [osc_orbital_elements,theta_osc]=mean2osc_m(x,param)

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
e=x(2,1);
Epw = x(6,1)';
Mo=Epw;
for l =1:1000
    DeltaEpw = -(Mo-Epw+e*sin(Epw)) /(-1+e*cos(Epw));
    Epw = Epw + DeltaEpw;
    if DeltaEpw < 10e-14
        break
    end
end
E=Epw;
theta=2*atan(sqrt((1+e)/(1-e))*tan(E/2));
if e > 1
    warning('hyperbolic orbit in mean2osc_m');
    fprintf('e: %0.3e \t theta: %0.3e \n', e, theta);
end

oeosc = mean2osc ([x(1:end-1,1);theta],param);

%compute the osculating mean anomaly
theta_osc=oeosc(6);
e_osc=oeosc(2);
E_osc=2*atan(sqrt((1+e_osc)/(1-e_osc))^(-1)*tan(theta_osc/2));
M_osc=E_osc-e_osc*sin(E_osc);

osc_orbital_elements=[oeosc(1);oeosc(2);oeosc(3);oeosc(4);oeosc(5);M_osc];

