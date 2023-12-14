function oeosc = mean2osc_vec (oemean,param)

% convert mean classical orbital elements to
% osculating classical orbital elements

% input

%  oemean(1) = mean semimajor axis (kilometers)
%  oemean(2) = mean orbital eccentricity (non-dimensional)
%              (0 <= eccentricity < 1)
%  oemean(3) = mean orbital inclination (radians)
%              (0 <= inclination <= pi)
%  oemean(4) = mean right ascension of ascending node (radians)
%              (0 <= raan <= 2 pi)
%  oemean(5) = mean argument of perigee (radians)
%              (0 <= argument of perigee <= 2 pi)
%  oemean(6) = mean true anomaly (radians)
%              (0 <= true anomaly <= 2 pi)

% output

%  oeosc(1) = osculating semimajor axis (kilometers)
%  oeosc(2) = osculating orbital eccentricity (non-dimensional)
%             (0 <= eccentricity < 1)
%  oeosc(3) = osculating orbital inclination (radians)
%             (0 <= inclination <= pi)
%  oeosc(4) = osculating right ascension of ascending node (radians)
%             (0 <= raan <= 2 pi)
%  oeosc(5) = osculating argument of perigee (radians)
%             (0 <= argument of perigee <= 2 pi)
%  oeosc(6) = osculating true anomaly (radians)
%             (0 <= true anomaly <= 2 pi)

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pi2 = 2.0 * pi;

% compute j2 effect on orbital elements

doe = delm(oemean,param);

oeosc = oemean + doe;
    
oeosc(:,3:6) = mod(oeosc(:,3:6) + pi2, pi2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function aout = anomly (nflg, ain, e)

% orbital anomaly function

switch nflg
    
    case {1 2}
        
        am = ain;
        
    case {3 4}
        
        ea = ain;
        
    case {5 6}
        
        f = ain;
        
end

if (nflg == 1 || nflg == 2)
    n_sats = numel(am);
    if n_sats<1e5 %C version runs faster for small (<1e5) number of objects
        ea = kepler1_C(am, e, pi, n_sats);
    else
        ea = kepler1_vec(am, e);
    end
end

if (nflg == 2 || nflg == 4)
    
    f = 2.0 * atan(sqrt((1.0 + e) ./ (1.0 - e)) .* tan(ea / 2.0));
    
end

if (nflg == 5 || nflg == 6)
    
    ea = 2.0 * atan(sqrt((1.0 - e) ./ (1.0 + e)) .* tan(f / 2.0));
    
end

if (nflg == 5 || nflg == 3)
    
    am = ea - e .* sin(ea);
end

switch nflg
    
    case {1 6}
        
        aout = ea;
        
    case {2 4}
        
        aout = f;
        
    case {3, 5}
        
        aout = am;
        
end

%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%

function dx = delm (x,param)

mu = param.mu;
req = param.req;
j2 = param.j2;

% j2 effect on orbital elements

% global mu req j2

pi2 = 2.0 * pi;

a = x(:,1);
e = x(:,2);
ai = x(:,3);
an = x(:,4);
w = x(:,5);
am = x(:,6);

f = anomly(2, am, e);

f = mod(f + pi2, pi2);

am = mod(am + pi2, pi2);

si = sin(ai);
ci = cos(ai);

ti = si ./ ci;

si2 = si .* si;
ci2 = ci .* ci;

sf = sin(f);
cf = cos(f);

s2f = sin(2.0 * f);

u = f + w;

e2 = e .* e;

esf = e .* sf;

d1 = 1.0 - e2;

d2 = sqrt(d1);

d3 = e .* cf;

d4 = 1.0 + d3;

d42 = d4 .* d4;

d5 = 1.0 + d2;

d6 = (3.0 * ci2 - 1.0) ./ d5;

p = a .* d1;

d7 = sqrt(mu ./ p);

r = p ./ d4;

rdot = d7 .* esf;

twou = 2.0 * u;
twow = 2.0 * w;

s2u = sin(twou);
c2u = cos(twou);

sf2w = sin(f + twow);
d8 = 3.0 * f + twow;
s3f2w = sin(d8);

cf2w = cos(f + twow);
c3f2w = cos(d8);

q1 = j2 * (req ./ p).^2;

di = 0.75 * q1 .* si .* ci .* (c2u + e .* cf2w + e / 3.0 .* c3f2w);

dp = 2.0 * p .* ti .* di;

dummy1 = f - am + esf - 0.5 * s2u - 0.5 * e .* sf2w - e .* s3f2w / 6.0;

dn = -1.5 * q1 .* ci .* dummy1;

dr = -0.25 * p .* q1 .* ((3.0 * ci2 - 1.0) ...
    .* (2.0 * d2 ./ d4 + d3 ./ d5 + 1.0) - si2 .* c2u);

drdot = 0.25 * d7 .* q1 .* (d6 .* esf .* (d2 .* d5 + d42) ...
    - 2.0 * si2 .* d42 .* s2u);

du = -0.125 * q1 .* (6.0 * (1.0 - 5.0 * ci2) ...
    .* (f - am) + 4.0 * esf .* ((1.0 - 6.0 * ci2) - d6) ...
    - d6 .* e2 .* s2f + 2.0 * (5.0 * ci2 - 2.0) .* e .* sf2w ...
    + (7.0 * ci2 - 1.0) .* s2u + 2.0 * ci2 .* e .* s3f2w);

pnw = p + dp;

ainw = ai + di;
annw = an + dn;

rnw = r + dr;

rdotnw = rdot + drdot;

unw = u + du;

aa = pnw ./ rnw - 1.0;
bb = sqrt(pnw ./ mu) .* rdotnw;

enw2 = aa .* aa + bb .* bb;
enw = sqrt(enw2);

%xfnw = atan3(bb, aa);
xfnw = atan2(bb, aa);

anw = pnw ./ (1.0 - enw2);
wnw = unw - xfnw;

amnw = anomly(5, xfnw, enw);

dx = zeros(numel(a),6);
dx(:,1) = anw - a;
dx(:,2) = enw - e;
dx(:,3) = ainw - ai;
dx(:,4) = annw - an;
dx(:,5) = wnw - w;
dx(:,6) = amnw - am;


