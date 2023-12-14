function oemean = osc2mean (oeosc,param)

% convert osculating classical orbital elements
% to mean classical orbital elements

% input

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

% output

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

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

req = param.req;
j2 = param.j2;

pi2 = 2.0 * pi;

oetmp = zeros(6, 1);

for i = 1:1:6
    
    oetmp(i) = oeosc(i);
    
end

a = sin(oetmp(6)) * sqrt(1.0 - oetmp(2) * oetmp(2));

b = oetmp(2) + cos(oetmp(6));

% eccentric anomaly

eanom = atan3(a, b);

% mean anomaly

oetmp(6) = mod(eanom - oetmp(2) * sin(eanom), pi2);

aos = oetmp(1);
eos = oetmp(2);
ios = oetmp(3);
ranos = oetmp(4);
apos = oetmp(5);
maos = oetmp(6);

aa = 1.0 / 3.0 - 0.5 * sin(ios)^2;

bb = 0.5 * sin(ios)^2;

if (eos < 0.01)
    
    lamos = mod(maos + apos, pi2);
    
    zos = eos * cos(apos);
    
    etaos = eos * sin(apos);
    
    sl = sin(lamos);
    cl = cos(lamos);
    
    s2l = sin(2 * lamos);
    c2l = cos(2 * lamos);
    
    s3l = sin(3 * lamos);
    c3l = cos(3 * lamos);
    
    s4l = sin(4 * lamos);
    c4l = cos(4 * lamos);
    
    s2i = sin(2 * ios);
    
    ci = cos(ios);
    
    am = aos;
    im = ios;
    ranm = ranos;
    mam = maos;
    lamm = lamos;
    zm = zos;
    etam = etaos;
    
    for n = 1:1:5
        
        asp = 3 * j2 * req * req / am * (bb * c2l + (1 - 3.5 * bb) ...
            * zm * cl + (1 - 2.5 * bb) * etam * sl + 3.5 * bb ...
            * (zm * c3l + etam * s3l));
        
        am = aos-asp;
        
        isp = 3 * j2 / 8 * req * req / am^2 * s2i * (c2l - zm * cl ...
            + etam * sl + 7/3 * (zm * c3l + etam * s3l));
        
        im = ios - isp;
        
        ci = cos(im);
        
        s2i = sin(2 * im);
        
        bb = 0.5 * sin(im)^2;
        
        ransp = 1.5 * j2 * req * req / am^2 * ci ...
            * (0.5 * s2l - 3.5 * zm * sl + 2.5 * etam * cl ...
            + 7/6 * (zm * s3l - etam * c3l));
        
        ranm = ranos - ransp;
        
        lamsp = 1.5 * j2 * req * req / (am * am) ...
            * (-0.5 * (1 - 5 * bb) * s2l + (7 - 77/4 * bb) ...
            * zm * sl - (6 - 55/4 * bb) * etam * cl - (7/6 - 77/12 * bb) ...
            * (zm * s3l - etam * c3l));
        
        lamm = lamos - lamsp;
        
        sl = sin(lamm);
        cl = cos(lamm);
        s2l = sin(2*lamm);
        c2l = cos(2*lamm);
        s3l = sin(3*lamm);
        c3l = cos(3*lamm);
        s4l = sin(4*lamm);
        c4l = cos(4*lamm);
        
        zsp = 1.5 * j2 * req * req / am^2 * ((1 - 2.5 * bb) * cl ...
            + 7/6 * bb * c3l + (1.5 - 5 * bb) * zm * c2l + (2 - 3 * bb) ...
            * etam * s2l + 17/4 * bb * (zm * c4l + etam * s4l));
        
        zm = zos - zsp;
        
        etasp = 1.5 * j2 * req * req / am^2 * ((1 - 3.5 * bb) * sl ...
            + 7/6 * bb * s3l + (1 - 6 * bb) * zm * s2l - (1.5 - 4 * bb) ...
            * etam * c2l + 17/4 * bb * (zm * s4l - etam * c4l));
        
        etam = etaos - etasp;
    end
    
    em = sqrt(etam^2 + zm^2);
    
    apm = 0;
    
    if (em > 1.0e-8)
        
        apm = atan2(etam, zm);
        
    end
    
    % "mean" mean anomaly
    
    mam = mod(lamm - apm, pi2);
    
else
    
    pm = aos*(1-eos^2);
    am = aos;
    em = eos;
    im = ios;
    apm = apos;
    ranm = ranos;
    mam = maos;
    
    [~, tam] = kepler1(mam, em);
    
    um = mod(apm+tam, pi2);
    
    hm = pm/(1+em*cos(tam));
    
    for n = 1:1:5
        
        asp = 3*j2*req*req/am*((am/hm)^3*(aa+bb*cos(2*um)) ...
            -aa*(1-em^2)^(-1.5));
        
        am = aos-asp;
        
        isp = 3/8*j2*req^2/pm^2*sin(2*im)*(cos(2*um) ...
            +em*cos(tam+2*apm)+1/3*em*cos(3*tam+2*apm));
        
        im = ios-isp;
        
        aa = 1/3-0.5*sin(im)^2;
        
        bb = 0.5*sin(im)^2;
        
        esp = 1.5*j2*req*req/am^2*((1-em^2)/em*((am/hm)^3 ...
            *(aa+bb*cos(2*um))-aa*(1-em^2)^(-1.5))-bb/(em*(1-em^2)) ...
            *(cos(2*um)+em*cos(tam+2*apm)+em*cos(3*tam+2.*apm)/3));
        
        em = eos-esp;
        
        pm = am * (1.0 - em^2);
        
        hm = pm / (1.0 + em * cos(tam));
        
        tam = mod(tam, 2.0 * pi);
        
        mam = mod(mam, 2.0 * pi);
        
        if ((abs(tam-pi) <= 1.0e-06) || (abs(mam-pi) <= 1.0e-06) || ...
                (abs(tam) <= 1.0e-06) || (abs(mam) <= 1.0e-06))
            
            eqoc = 0;
            
        else
            
            eqoc = tam-mam;
            
        end
        
        ransp = -1.5*j2*(req/pm)^2*cos(im)*(eqoc+em*sin(tam) ...
            -0.5*sin(2*um)-0.5*em*sin(tam+2*apm)-1/6*em*sin(3*tam+2*apm));
        
        ranm = ranos-ransp;
        
        apsp = 1.5*j2*(req/pm)^2*((2-5*bb)*(eqoc+em*sin(tam)) ...
            +(1-3*bb)*((1-0.25*em*em)*sin(tam)/em+0.5*sin(2*tam) ...
            +em*sin(3*tam)/12)-(0.5*bb+(0.5-15/8*bb)*em^2)/em*sin(tam+2*apm) ...
            +em/8*bb*sin(tam-2*apm)-0.5*(1-5*bb)*sin(2*um)+(7/6*bb ...
            -1/6*em^2*(1-19/4*bb))/em*sin(3*tam+2*apm) ...
            +0.75*bb*sin(4*tam+2*apm)+em/8*bb*sin(5*tam+2*apm));
        
        apm = apos-apsp;
        
        masp = 1.5*j2*(req/pm)^2*sqrt(1-em^2)/em*(-(1-3*bb)*((1-em^2/4) ...
            *sin(tam)+em/2*sin(2*tam)+em^2/12*sin(3*tam)) ...
            +bb*(0.5*(1+1.25*em^2)*sin(tam+2*apm)-em^2/8*sin(tam-2*apm) ...
            -7/6*(1-em^2/28)*sin(3*tam+2*apm)-0.75*em*sin(4*tam+2*apm) ...
            -em^2/8*sin(5*tam+2*apm)));
        
        mam = maos-masp;
        
        [~, tam] = kepler1(mam, em);
        
        um = mod(apm + tam, pi2);
    end
end

oemean(1) = am;
oemean(2) = em;
oemean(3) = im;
oemean(4) = mod(ranm, pi2);
oemean(5) = mod(apm, pi2);
oemean(6) = mod(mam, pi2);

[~, ta] = kepler1(oemean(6), em);

oemean(6) = ta;

