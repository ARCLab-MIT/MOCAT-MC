%
% ------------------------------------------------------------------------------
%
%                           function eq2rv
%
%  this function finds the classical orbital elements given the equinoctial
%    elements.
%
%  author        : david vallado                  719-573-2600    9 jun 2002
%
%  revisions
%    vallado     - fix elliptical equatorial orbits case         19 oct 2002
%    vallado     - add constant file use                         29 jun 2003
%
%  inputs          description                    range / units
%    af          -
%    ag          -
%    n           - mean motion                    rad
%    meanlon     - mean longitude                 rad
%    chi         -
%    psi         -
%
%  outputs       :
%    r           - eci position vector            km
%    v           - eci velocity vector            km/s
%
%  locals        :
%    temp        - temporary variable
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
%  coupling      :
%
%  references    :
%    vallado       2001, xx
%    chobotov            30
%
% [r,v] = eq2rv( af, ag, n, meanlon, chi, psi);
% ------------------------------------------------------------------------------

function [r,v] = eq2rv( af, ag, meanlon, n, chi, psi);

        % -------------------------  implementation   -----------------
        constmath;
        constastro;

        arglat  = 999999.1;
        lonper  = 999999.1;
        truelon = 999999.1;

        a = (mu/n^2)^(1.0/3.0);

        ecc = sqrt (af^2 + ag^2);

        p = a * (1.0 - ecc*ecc);

        incl = 2.0 * atan( sqrt(chi^2 + psi^2) );

        % -------- setup retrograde factor ----------------------------
        fr = 1;
        % -------- set this so it only affects i = 180 deg orbits!! ---
        if abs(incl-pi) < small
            fr = -1;
          end

        omega = atan2( chi, psi);

        argp = atan2( fr*ag,af ) - atan2( chi,psi );

        if ( ecc < small )
            % ----------------  circular equatorial  ------------------
            if (incl<small) | ( abs(incl-pi)< small )
                argp = 0.0;
                omega= 0.0;
%                truelon = nu;
              else
                % --------------  circular inclined  ------------------
                argp= 0.0;
%                arglat = nu;
              end
          else
            % ---------------  elliptical equatorial  -----------------
            if ( ( incl<small) | (abs(incl-pi)<small) )
%                argp = lonper;
                omega= 0.0;
              end
          end

        m = meanlon - omega - argp;
        m = rem (m+twopi,twopi);

        [e0,nu] = newtonm ( ecc,m );

        % ----------  fix for elliptical equatorial orbits ------------
        if ( ecc < small )
           % ----------------  circular equatorial  ------------------
           if (incl<small) | ( abs(incl-pi)< small )
               argp    = undefined;
               omega   = undefined;
               truelon = nu;
             else
               % --------------  circular inclined  ------------------
               argp  = undefined;
               arglat= nu;
             end
           nu   = undefined;
         else
           % ---------------  elliptical equatorial  -----------------
           if ( ( incl < small) | (abs(incl-pi) < small) )
               lonper = argp;
               argp    = undefined;
               omega   = undefined;
             end
         end

        % -------- now convert back to position and velocity vectors
        [r,v] = coe2rv(p,ecc,incl,omega,argp,nu,arglat,truelon,lonper);

