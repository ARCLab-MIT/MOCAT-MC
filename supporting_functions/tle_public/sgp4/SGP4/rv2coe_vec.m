%
% ------------------------------------------------------------------------------
%
%                           function rv2coe
%
%  this function finds the classical orbital elements given the geocentric
%    equatorial position and velocity vectors.
%
%  author        : david vallado                  719-573-2600   21 jun 2002
%
%  revisions
%    vallado     - fix special cases                              5 sep 2002
%    vallado     - delete extra check in inclination code        16 oct 2002
%    vallado     - add constant file use                         29 jun 2003
%    vallado     - add mu                                         2 apr 2007
%
%  inputs          description                    range / units
%    r           - ijk position vector            km
%    v           - ijk velocity vector            km / s
%    mu          - gravitational parameter        km3 / s2
%
%  outputs       :
%    p           - semilatus rectum               km
%    a           - semimajor axis                 km
%    ecc         - eccentricity
%    incl        - inclination                    0.0  to pi rad
%    omega       - longitude of ascending node    0.0  to 2pi rad
%    argp        - argument of perigee            0.0  to 2pi rad
%    nu          - true anomaly                   0.0  to 2pi rad
%    m           - mean anomaly                   0.0  to 2pi rad
%    arglat      - argument of latitude      (ci) 0.0  to 2pi rad
%    truelon     - true longitude            (ce) 0.0  to 2pi rad
%    lonper      - longitude of periapsis    (ee) 0.0  to 2pi rad
%
%  locals        :
%    hbar        - angular momentum h vector      km2 / s
%    ebar        - eccentricity     e vector
%    nbar        - line of nodes    n vector
%    c1          - v**2 - u/r
%    rdotv       - r dot v
%    hk          - hk unit vector
%    sme         - specfic mechanical energy      km2 / s2
%    i           - index
%    e           - eccentric, parabolic,
%                  hyperbolic anomaly             rad
%    temp        - temporary variable
%    typeorbit   - type of orbit                  ee, ei, ce, ci
%
%  coupling      :
%    mag         - magnitude of a vector
%    angl        - find the angl between two vectors
%    newtonnu    - find the mean anomaly
%
%  references    :
%    vallado       2007, 121, alg 9, ex 2-5
%
% [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe (r,v,mu);
% ------------------------------------------------------------------------------

function [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper] = rv2coe_vec(r,v,mu)
%         constmath;
        small = 0.00000001;
        infinite  = 999999.9;
        undefined = 999999.1;
        twopi  = 2.0 * pi;
        halfpi = pi * 0.5;
%         constastro;
%         mu         = 398600.4418;      % km3/s2

        %Preallocation
        num_el = size(r,1);
        a = zeros(num_el,1);
        p = zeros(num_el,1);
        ecc = zeros(num_el,1);
        incl = zeros(num_el,1);
        omega = zeros(num_el,1);
        argp = zeros(num_el,1);
        nu = zeros(num_el,1);
        m = zeros(num_el,1);
        arglat = zeros(num_el,1);
        truelon = zeros(num_el,1);
        lonper = zeros(num_el,1);
        
        % -------------------------  implementation   -----------------
        magr = sqrt(sum(r.^2,2));
        magv = sqrt(sum(v.^2,2));
        
        % ------------------  find h n and e vectors   ----------------
        [hbar] = cross_vec( r,v );
        magh = sqrt(sum(hbar.^2,2));
        check_small = magh>small;
        idx_small = find(check_small);
        sum_check_small = sum(check_small);
        
        nbar = zeros(sum_check_small,3);
        nbar(:,1) = -hbar(idx_small,2);
        nbar(:,2) = hbar(idx_small,1);
        magn = sqrt(sum(nbar.^2,2));
        magr_check = magr(idx_small);
        magv_check = magv(idx_small);
        magh_check = magh(idx_small);
        r_check = r(idx_small,:);
        v_check = v(idx_small,:);
        c1 = magv_check.*magv_check - mu./magr_check;
        rdotv = dot( r_check,v_check,2 );
        ebar = [(c1.*r_check(:,1) - rdotv.*v_check(:,1)),(c1.*r_check(:,2) - rdotv.*v_check(:,2)),(c1.*r_check(:,3) - rdotv.*v_check(:,3))]/mu;
        ecci = sqrt(sum(ebar.^2,2));

            % ------------  find a e and semi-latus rectum   ----------
        sme = ( magv_check.*magv_check*0.5  ) - ( mu ./magr_check);
        check_sme = abs(sme)>small;
        a(idx_small(check_sme)) = -mu  ./ (2.0 *sme(check_sme));
        a(idx_small(~check_sme)) = infinite;
        p(idx_small) = magh_check.*magh_check/mu;

            % -----------------  find inclination   -------------------
        hk = hbar(idx_small,3)./magh_check;
        incli = acos( hk );

            % --------  determine type of orbit for later use  --------
            % ------ elliptical, parabolic, hyperbolic inclined -------
        typeorbit = zeros(sum_check_small,1); %[0,1,2,3]=[ei,ee,ce,ci], ei is default
        check_ecc = ecci < small;
        check_notecc = ~check_ecc;
        idx_ecc = find(check_ecc);
        idx_notecc = find(check_notecc);
        check_in = incli(idx_ecc)<small | abs(incli(idx_ecc)-pi)<small;
                % ----------------  circular equatorial ---------------
        typeorbit(idx_ecc(check_in)) = 2;
                    % --------------  circular inclined ---------------
        typeorbit(idx_ecc(~check_in)) = 3;
        % - elliptical, parabolic, hyperbolic equatorial --
        typeorbit(idx_notecc(incli(idx_notecc)<small | abs(incli(idx_notecc)-pi)<small)) = 1;
            % ----------  find longitude of ascending node ------------
        check_magn = magn > small;

        temp = nbar(check_magn,1) ./ magn(check_magn);
        check_temp = abs(temp) > 1.0;
        temp(check_temp) = sign(temp(check_temp));
        omega_temp = acos( temp );
        check_nbar = nbar(check_magn,2) < 0.0;
        omega_temp(check_nbar) = twopi - omega_temp(check_nbar);
        omega(idx_small(check_magn)) = omega_temp;
        omega(idx_small(~check_magn)) = undefined;
            % ---------------- find argument of perigee ---------------
        check_ei = typeorbit==0;

        argp_temp = angl_vec( nbar(check_ei,:),ebar(check_ei,:));
        check_ebar = ebar(check_ei,3) < 0.0;
        argp_temp(check_ebar) = twopi - argp_temp(check_ebar);
        argp(idx_small(check_ei)) = argp_temp;
        argp(idx_small(~check_ei)) = undefined;
            % ------------  find true anomaly at epoch    -------------
        check_e = typeorbit<1.5;
        
        nu_temp =  angl_vec( ebar(check_e,:),r_check(check_e,:));
        check_rdotv = ( rdotv(check_e) < 0.0  );
        nu_temp(check_rdotv) = twopi - nu_temp(check_rdotv);
        nu(idx_small(check_e)) = nu_temp;
        nu(idx_small(~check_e)) = undefined;
            % ----  find argument of latitude - circular inclined -----
        check_ci = typeorbit==3;

        arglat_temp =  angl_vec( nbar(check_ci,:),r_check(check_ci,:));
        check_r = ( r_check(check_ci,3) < 0.0  );
        arglat_temp(check_r) = twopi - arglat_temp(check_r);
        arglat(idx_small(check_ci)) = arglat_temp;
        m(idx_small(check_ci)) = arglat_temp;
        arglat(idx_small(~check_ci)) = undefined;
            % -- find longitude of perigee - elliptical equatorial ----
        check_ee = typeorbit(idx_notecc)==1;
        idx_ee = idx_notecc(check_ee);
        
        temp = ebar(idx_ee,1)./ecci(idx_ee);
        check_temp = abs(temp) > 1.0;
        temp(check_temp) = sign(temp(check_temp));
        lonper_temp = acos( temp );
        check_ebar = ebar(idx_ee,2) < 0.0;
        lonper_temp(check_ebar) = twopi - lonper_temp(check_ebar);
        check_incl = incli(idx_ee) > halfpi;
        lonper_temp(check_incl) = twopi - lonper_temp(check_incl);
        lonper(idx_small(idx_ee)) = lonper_temp;
        lonper(idx_small(idx_notecc(~check_ee))) = undefined;
            % -------- find true longitude - circular equatorial ------
        check_magr = magr_check & typeorbit==2;
        
        temp = r_check(check_magr,1)./magr_check(check_magr);
        check_temp = abs(temp) > 1.0;
        temp(check_temp) = sign(temp(check_temp));
        truelon_temp = acos( temp );
        check_r = r_check(check_magr,2) < 0.0;
        truelon_temp(check_r) = twopi - truelon_temp(check_r);
        check_incl = incli(check_r) > halfpi;
        truelon_temp(check_incl) = twopi - truelon_temp(check_incl);
        truelon(idx_small(check_magr)) = truelon_temp;
        m(idx_small(check_magr)) = truelon_temp;
        truelon(idx_small(~check_magr)) = undefined;
            % ------------ find mean anomaly for all orbits -----------
        [~,m(idx_small(check_e))] = newtonnu_vec(ecci(check_e), nu(idx_small(check_e)) );
        ecc(idx_small)  = ecci;
        incl(idx_small) = incli;

        idx_notsmall = find(~check_small);
        p(idx_notsmall)       = undefined;
        a(idx_notsmall)       = undefined;
        ecc(idx_notsmall)     = undefined;
        incl(idx_notsmall)    = undefined;
        omega(idx_notsmall)   = undefined;
        argp(idx_notsmall)    = undefined;
        nu(idx_notsmall)      = undefined;
        m(idx_notsmall)       = undefined;
        arglat(idx_notsmall)  = undefined;
        truelon(idx_notsmall) = undefined;
        lonper(idx_notsmall)  = undefined;

