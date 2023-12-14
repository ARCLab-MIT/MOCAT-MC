function [eanom] = kepler1_vec (manom, ecc)

% solve Kepler's equation for circular,
% elliptic and hyperbolic orbits

% Danby's method

% input

%  manom = mean anomaly (radians)
%  ecc   = orbital eccentricity (non-dimensional)

% output

%  eanom = eccentric anomaly (radians)
%  tanom = true anomaly (radians)

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% global niter
% niter=0;
% define convergence criterion

ktol = 1.0e-10;

pi2 = 2.0 * pi;

xma = manom - pi2 * fix(manom/pi2);

% tanom = zeros(size(ecc));
eanom = zeros(size(ecc));
% initial guess
check_ecc_0 = ecc == 0; %circular orbit
check_ecc_above_1 = ecc >= 1.0; %hyperbolic orbit
check_ecc_hyp = find(check_ecc_above_1);
check_ecc_ellip = find(~check_ecc_above_1 & ~check_ecc_0); %elliptical, not circular

% if (ecc == 0.0)
    
% circular orbit
% tanom(check_ecc_0) = xma(check_ecc_0);
eanom(check_ecc_0) = xma(check_ecc_0);
    
% elseif (ecc < 1.0)
% elliptic orbit
eanom(check_ecc_ellip) = xma(check_ecc_ellip) + 0.85 * sign(sin(xma(check_ecc_ellip))) .* ecc(check_ecc_ellip);

% else
    % hyperbolic orbit
eanom(check_ecc_hyp) = log(2.0 * xma(check_ecc_hyp) ./ ecc(check_ecc_hyp) + 1.8);
% end

% perform iterations

niter = 0;

check_ecc_ellip_temp = check_ecc_ellip;
check_ecc_hyp_temp = check_ecc_hyp;
while ~isempty(check_ecc_ellip_temp) || ~isempty(check_ecc_hyp_temp) %while some elements in f are still above tolerance
    
%     if (ecc < 1)

% elliptic orbit

se = ecc(check_ecc_ellip_temp) .* sin(eanom(check_ecc_ellip_temp));
ce = ecc(check_ecc_ellip_temp) .* cos(eanom(check_ecc_ellip_temp));

fe = eanom(check_ecc_ellip_temp) - se - xma(check_ecc_ellip_temp);
fpe = 1 - ce;
fppe = se;
fpppe = ce;
        
%     else

% hyperbolic orbit

sh = ecc(check_ecc_hyp_temp) .* sinh(eanom(check_ecc_hyp_temp));
ch = ecc(check_ecc_hyp_temp) .* cosh(eanom(check_ecc_hyp_temp));

fh = sh - eanom(check_ecc_hyp_temp) - xma(check_ecc_hyp_temp);
fph = ch - 1;
fpph = sh;
fppph = ch;
        
%     end

    niter = niter + 1;

    % check for convergence
    idx_above_e = abs(fe) > ktol;
    idx_above_h = abs(fh) > ktol;
    check_ecc_ellip_temp = check_ecc_ellip_temp(idx_above_e);
    check_ecc_hyp_temp = check_ecc_hyp_temp(idx_above_h);
    if (niter > 20 || (isempty(check_ecc_ellip_temp) && isempty(check_ecc_hyp_temp)))
        
        break;
        
    end

    check_f = [check_ecc_ellip_temp;check_ecc_hyp_temp];
    f =  [fe(idx_above_e);fh(idx_above_h)]; 
    fp =  [fpe(idx_above_e);fph(idx_above_h)];
    fpp =  [fppe(idx_above_e);fpph(idx_above_h)];
    fppp =  [fpppe(idx_above_e);fppph(idx_above_h)];

    % update eccentric anomaly

    delta = -f ./ fp;

    deltastar = -f ./ (fp + 0.5 * delta .* fpp);

    deltak = -f ./ (fp + 0.5 * deltastar .* fpp ...
        + deltastar .* deltastar .* fppp / 6);

    eanom(check_f) = eanom(check_f) + deltak;
end

if (niter > 20)
    
    clc; home;
    
    fprintf('\n\n   more than 20 iterations in kepler1_vec \n\n');
    
    pause;
    
    return;
end

% % compute true anomaly
% 
% % if (ecc < 1)
% % elliptic orbit
% 
% stae = sqrt(1 - ecc(check_ecc_ellip) .* ecc(check_ecc_ellip)) .* sin(eanom(check_ecc_ellip));
% ctae = cos(eanom(check_ecc_ellip)) - ecc(check_ecc_ellip);
% % else
% % hyperbolic orbit
% 
% stah = sqrt(ecc(check_ecc_hyp) .* ecc(check_ecc_hyp) - 1) .* sinh(eanom(check_ecc_hyp));
% ctah = ecc(check_ecc_hyp) - cosh(eanom(check_ecc_hyp));
% % end
% 
% % tanom = atan3(sta, cta);
% 
% tanom([check_ecc_ellip;check_ecc_hyp]) = atan2([stae;stah], [ctae;ctah]);

