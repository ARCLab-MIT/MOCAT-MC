function [eanom, tanom] = kepler1 (manom, ecc)

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

global niter

% define convergence criterion

ktol = 1.0e-10;

pi2 = 2.0 * pi;

xma = manom - pi2 * fix(manom/pi2);

% initial guess

if (ecc == 0.0)
    
    % circular orbit

    tanom = xma;
    eanom = xma;
    return;
    
elseif (ecc < 1.0)
    % elliptic orbit

    eanom = xma + 0.85 * sign(sin(xma)) * ecc;

else
    % hyperbolic orbit

    eanom = log(2.0 * xma / ecc + 1.8);
end

% perform iterations

niter = 0;

while(1)
    
    if (ecc < 1)

        % elliptic orbit

        s = ecc * sin(eanom);
        c = ecc * cos(eanom);

        f = eanom - s - xma;
        fp = 1 - c;
        fpp = s;
        fppp = c;
        
    else

        % hyperbolic orbit

        s = ecc * sinh(eanom);
        c = ecc * cosh(eanom);

        f = s - eanom - xma;
        fp = c - 1;
        fpp = s;
        fppp = c;
        
    end

    niter = niter + 1;

    % check for convergence

    if (abs(f) <= ktol || niter > 20)
        
        break;
        
    end

    % update eccentric anomaly

    delta = -f / fp;

    deltastar = -f / (fp + 0.5 * delta * fpp);

    deltak = -f / (fp + 0.5 * deltastar * fpp ...
        + deltastar * deltastar * fppp / 6);

    eanom = eanom + deltak;
end

if (niter > 20)
    
    home;
    
    fprintf('\n\n   more than 20 iterations in kepler1 \n\n');
    
    pause;
    
    return;
end

% compute true anomaly

if (ecc < 1)
    % elliptic orbit

    sta = sqrt(1 - ecc * ecc) * sin(eanom);
    cta = cos(eanom) - ecc;
else
    % hyperbolic orbit

    sta = sqrt(ecc * ecc - 1) * sinh(eanom);
    cta = ecc - cosh(eanom);
end

tanom = atan3(sta, cta);


