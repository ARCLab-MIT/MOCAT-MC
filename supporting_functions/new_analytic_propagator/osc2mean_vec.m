function oemean = osc2mean_vec(oeosc,param)

% convert osculating classical orbital elements
% to mean classical orbital elements

% input

%  oeosc(:,1) = osculating semimajor axis (kilometers)
%  oeosc(:,2) = osculating orbital eccentricity (non-dimensional)
%             (0 <= eccentricity < 1)
%  oeosc(:,3) = osculating orbital inclination (radians)
%             (0 <= inclination <= pi)
%  oeosc(:,4) = osculating right ascension of ascending node (radians)
%             (0 <= raan <= 2 pi)
%  oeosc(:,5) = osculating argument of perigee (radians)
%             (0 <= argument of perigee <= 2 pi)
%  oeosc(:,6) = osculating true anomaly (radians)
%             (0 <= true anomaly <= 2 pi)

% output

%  oemean(:,1) = mean semimajor axis (kilometers)
%  oemean(:,2) = mean orbital eccentricity (non-dimensional)
%              (0 <= eccentricity < 1)
%  oemean(:,3) = mean orbital inclination (radians)
%              (0 <= inclination <= pi)
%  oemean(:,4) = mean right ascension of ascending node (radians)
%              (0 <= raan <= 2 pi)
%  oemean(:,5) = mean argument of perigee (radians)
%              (0 <= argument of perigee <= 2 pi)
%  oemean(:,6) = mean true anomaly (radians)
%              (0 <= true anomaly <= 2 pi)

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

req = param.req;
j2 = param.j2;

pi2 = 2.0 * pi;

eos = oeosc(:,2);
M = oeosc(:,6);

a = sin(M) .* sqrt(1.0 - eos.^2);

b = eos + cos(M);

% eccentric anomaly
eanom = atan2(a, b);

aos = oeosc(:,1);
ios = oeosc(:,3);
ranos = oeosc(:,4);
apos = oeosc(:,5);

si2 = sin(ios).^2;

% mean anomaly
maos = mod(eanom - eos.*sin(eanom), pi2);

aa = 1.0/3.0 - 0.5*si2;

bb = 0.5 * si2;

% if (eos < 0.01)
check_e = eos<0.01;
find_e = find(check_e);

aos_e = aos(find_e);
eos_e = eos(find_e);
ios_e = ios(find_e);
ranos_e = ranos(find_e);
apos_e = apos(find_e);
maos_e = maos(find_e);
bb_e = bb(find_e);
    
lamos = mod(maos_e + apos_e, pi2);

zos = eos_e .* cos(apos_e);

etaos = eos_e .* sin(apos_e);

sl = sin(lamos);
cl = cos(lamos);

s2l = sin(2 * lamos);
c2l = cos(2 * lamos);

s3l = sin(3 * lamos);
c3l = cos(3 * lamos);

%s4l = sin(4*lamos);
%c4l = cos(4*lamos);

s2i = sin(2*ios_e);
%ci = cos(ios_e);

am_e = aos_e;
lamm = lamos;
zm = zos;
etam = etaos;

threej2req2 = 3*j2*req*req;
    
for n = 1:1:5

    asp_e = threej2req2./am_e.*(bb_e.*c2l + (1 - 3.5*bb_e) ...
        .*zm.*cl + (1 - 2.5*bb_e).*etam.*sl + 3.5*bb_e ...
        .*(zm.*c3l + etam.*s3l));

    am_e = aos_e-asp_e;
    am2_e = am_e.^2;
    
    threej2req2_over2_am2 = threej2req2/2./am2_e;

    isp_e = threej2req2_over2_am2./4.*s2i.*(c2l - zm.*cl ...
        + etam.*sl + 7/3*(zm.*c3l + etam.*s3l));

    im_e = ios_e - isp_e;
    ci = cos(im_e);
    s2i = sin(2*im_e);
    bb_e = 0.5*sin(im_e).^2;
    
    ransp_e = threej2req2_over2_am2.*ci ...
        .*(0.5*s2l - 3.5*zm.*sl + 2.5*etam.*cl ...
        + 7/6*(zm.*s3l - etam.*c3l));

    ranm_e = ranos_e - ransp_e;

    lamsp = threej2req2_over2_am2 ...
        .*(-0.5*(1 - 5*bb_e).*s2l + (7 - 77/4*bb_e) ...
        .*zm.*sl - (6 - 55/4*bb_e).*etam.*cl - (7/6 - 77/12*bb_e) ...
        .*(zm.*s3l - etam.*c3l));

    lamm = lamos - lamsp;

    sl = sin(lamm);
    cl = cos(lamm);
    s2l = sin(2*lamm);
    c2l = cos(2*lamm);
    s3l = sin(3*lamm);
    c3l = cos(3*lamm);
    s4l = sin(4*lamm);
    c4l = cos(4*lamm);

    zsp = threej2req2_over2_am2.*((1 - 2.5*bb_e).*cl ...
        + 7/6*bb_e.*c3l + (1.5 - 5*bb_e).*zm.*c2l + (2 - 3*bb_e) ...
        .*etam.*s2l + 17/4*bb_e.*(zm.*c4l + etam.*s4l));

    zm = zos - zsp;

    etasp = threej2req2_over2_am2.*((1 - 3.5*bb_e).*sl ...
        + 7/6*bb_e.*s3l + (1 - 6*bb_e).*zm.*s2l - (1.5 - 4*bb_e) ...
        .*etam.*c2l + 17/4*bb_e.*(zm.*s4l - etam.*c4l));

    etam = etaos - etasp;
end
    
em_e = sqrt(etam.^2 + zm.^2);

apm_e = zeros(size(em_e));
    
% if (em_note > 1.0e-8)
check_em = em_e > 1.0e-8;
apm_e(check_em) = atan2(etam(check_em), zm(check_em));
% end
    
% "mean" mean anomaly
mam_e = mod(lamm - apm_e, pi2);
    
% else (eos < 0.01)
find_note = find(~check_e);

aos_note = aos(find_note);
eos_note = eos(find_note);
ios_note = ios(find_note);
apos_note = apos(find_note);
ranos_note = ranos(find_note);
maos_note = maos(find_note);

am_note = aos_note;
em_note = eos_note;
em_note2 = em_note.^2;
one_minus_emnote = 1-em_note2;
pm_note = aos_note.*one_minus_emnote;
im_note = ios_note;
apm_note = apos_note;
mam_note = maos_note;

aa_note = aa(~check_e);
bb_note = bb(~check_e);

n_sats = numel(em_note);
if n_sats<1e5 %C version runs faster for small (<1e5) number of objects
    [~, tam_note] = kepler1_C_tanom(mam_note, em_note, pi, n_sats);
else
    [~, tam_note] = kepler1_vec_tanom(mam_note, em_note); %tanom
end

um = mod(apm_note+tam_note, pi2);

hm = pm_note./(1+em_note.*cos(tam_note));

for n = 1:1:5
    cos2um = cos(2*um);
    sin2um = sin(2*um);

    asp_note = threej2req2./am_note.*((am_note./hm).^3.*(aa_note+bb_note.*cos2um) ...
        -aa_note.*one_minus_emnote.^(-1.5));

    am_note = aos_note-asp_note;

    isp_note = threej2req2/8./pm_note.^2.*sin(2*im_note).*(cos2um ...
        +em_note.*cos(tam_note+2*apm_note)+1/3*em_note.*cos(3*tam_note+2*apm_note));

    im_note = ios_note-isp_note;
    
    si2_enote = sin(im_note).^2;
    
    aa_note = 1/3-0.5*si2_enote;

    bb_note = 0.5*si2_enote;

    esp_note = threej2req2/2./am_note.^2.*(one_minus_emnote./em_note.*((am_note./hm).^3 ...
        .*(aa_note+bb_note.*cos2um)-aa_note.*one_minus_emnote.^(-1.5))-bb_note./(em_note.*one_minus_emnote) ...
        .*(cos2um+em_note.*cos(tam_note+2*apm_note)+em_note.*cos(3*tam_note+2.*apm_note)/3));

    em_note = eos_note-esp_note;
    em_note2 = em_note.^2;
    one_minus_emnote = 1-em_note2;

    pm_note = am_note.*one_minus_emnote;

    hm = pm_note ./ (1.0 + em_note.*cos(tam_note));

    tam_note = mod(tam_note, pi2);

    mam_note = mod(mam_note, pi2);
    
    eqoc = zeros(size(mam_note));
    
    check_tam = (abs(tam_note-pi) <= 1.0e-06) | (abs(mam_note-pi) <= 1.0e-06) | (abs(tam_note) <= 1.0e-06) | (abs(mam_note) <= 1.0e-06);

%     if ((abs(tam_note-pi) <= 1.0e-06) || (abs(mam_note-pi) <= 1.0e-06) || ...
%             (abs(tam_note) <= 1.0e-06) || (abs(mam_note) <= 1.0e-06))

%         eqoc = 0;

%     else
    eqoc(~check_tam) = tam_note(~check_tam)-mam_note(~check_tam);

%         eqoc = tam_note-mam_note;

%     end
    
    threej2req2_over2pm2 = threej2req2/2./pm_note.^2;
    
    ransp_note = -threej2req2_over2pm2.*cos(im_note).*(eqoc+em_note.*sin(tam_note) ...
        -0.5*sin2um-0.5*em_note.*sin(tam_note+2*apm_note)-1/6*em_note.*sin(3*tam_note+2*apm_note));

    ranm_note = ranos_note-ransp_note;

    apsp = threej2req2_over2pm2.*((2-5*bb_note).*(eqoc+em_note.*sin(tam_note)) ...
        +(1-3*bb_note).*((1-0.25*em_note2).*sin(tam_note)./em_note+0.5*sin(2*tam_note) ...
        +em_note.*sin(3*tam_note)/12)-(0.5*bb_note+(0.5-15/8*bb_note).*em_note2)./em_note.*sin(tam_note+2*apm_note) ...
        +em_note/8.*bb_note.*sin(tam_note-2*apm_note)-0.5*(1-5*bb_note).*sin2um+(7/6*bb_note ...
        -1/6*em_note2.*(1-19/4*bb_note))./em_note.*sin(3*tam_note+2*apm_note) ...
        +0.75*bb_note.*sin(4*tam_note+2*apm_note)+em_note/8.*bb_note.*sin(5*tam_note+2*apm_note));

%     apm_note_old = apm_note;
    apm_note = apos_note-apsp;
%     diff_apm = sum(abs(apm_note-apm_note_old))

    masp = threej2req2_over2pm2.*sqrt(one_minus_emnote)./em_note.*(-(1-3*bb_note).*((1-em_note2/4) ...
        .*sin(tam_note)+em_note/2.*sin(2*tam_note)+em_note2/12.*sin(3*tam_note)) ...
        +bb_note.*(0.5*(1+1.25*em_note2).*sin(tam_note+2*apm_note)-em_note2/8.*sin(tam_note-2*apm_note) ...
        -7/6*(1-em_note2/28).*sin(3*tam_note+2*apm_note)-0.75*em_note.*sin(4*tam_note+2*apm_note) ...
        -em_note2/8.*sin(5*tam_note+2*apm_note)));
    
%     mam_note_old = mod(mam_note,pi2);
    mam_note = maos_note-masp;
%     diff_mam = sum(abs(mod(mam_note,pi2)-mam_note_old))

    n_sats = numel(em_note);
    if n_sats<1e5 %C version runs faster for small (<1e5) number of objects
        [~, tam_note] = kepler1_C_tanom(mam_note, em_note, pi, n_sats);
    else
        [~, tam_note] = kepler1_vec_tanom(mam_note, em_note); %tanom
    end
    
    um = mod(apm_note + tam_note, pi2);
end
% end

oemean = zeros(length(aos),6);
em = [em_e;em_note];
mod_mam = mod([mam_e;mam_note], pi2);

n_sats = numel(em);
if n_sats<1e5 %C version runs faster for small (<1e5) number of objects
    [~, tanom] = kepler1_C_tanom(mod_mam, em, pi, n_sats); %tanom
else
    [~, tanom] = kepler1_vec_tanom(mod_mam, em); %tanom
end

oemean([find_e;find_note],:) = [[am_e;am_note],em,[im_e;im_note],[ranm_e;ranm_note],[apm_e;apm_note],tanom];





