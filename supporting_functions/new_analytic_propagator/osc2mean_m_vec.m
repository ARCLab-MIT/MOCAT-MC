function [mean_orbital_elements,theta_mean,E_mean]=osc2mean_m_vec(x,param,theta_osc)

% input

%  x(1) = osc semimajor axis (kilometers)
%  x(2) = osc orbital eccentricity (non-dimensional)
%              (0 <= eccentricity < 1)
%  x(3) = osc orbital inclination (radians)
%              (0 <= inclination <= pi)
%  x(4) = osc right ascension of ascending node (radians)
%              (0 <= raan <= 2 pi)
%  x(5) = osc argument of perigee (radians)
%              (0 <= argument of perigee <= 2 pi)
%  x(6) = osc mean anomaly (radians)
%              (0 <= true anomaly <= 2 pi)

if isempty(theta_osc)
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
    theta_osc = 2*atan(sqrt_ep1_em1.*tan(Epw/2));
end

oemean = osc2mean_vec([x(:,1:end-1),theta_osc],param);
theta_mean = oemean(:,6);
e_mean = oemean(:,2);
%compute the osculating mean anomaly
E_mean = 2*atan(sqrt((1+e_mean)./(1-e_mean)).^(-1).*tan(theta_mean/2));
M_mean = E_mean-e_mean.*sin(E_mean);

mean_orbital_elements = [oemean(:,1:5),M_mean];