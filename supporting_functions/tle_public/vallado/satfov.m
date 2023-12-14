% ------------------------------------------------------------------------------
%
%                           function satfov
%
%  this function finds parameters relating to a satellite's fov.
%
%  author        : david vallado                  719-573-2600   27 may 2002
%
%  revisions
%                -
%
%  inputs          description                    range / units
%    incl        - inclination                    rad
%    az          - azimuth                        rad
%    slatgd      - geodetic latitude of sat       rad
%    slon        - longitude of sat               rad
%    salt        - altitudeof satellite           er
%    tfov        - total field of view            rad
%    etactr      - ctr where sensor looks         rad
%
%  outputs       :
%    fovmax      - maximum field of view          rad
%    totalrng    -
%    rhomax      -
%    rhomin      -
%    tgtlat      -
%    tgtlon      -
%
%  locals        :
%    r           -
%    etahopriz   -
%    rhohoriz    -
%    gamma       -
%    rho         -
%    fovmin      -
%    lat         -
%    lon         -
%    maxlat      -
%    minlkat     -
%    i           - index
%
%  coupling      :
%    path        - finds tgt location given initial location, range, and az
%
%  references    :
%    vallado       2001, 776-781, eq 11-8 to eq 11-13, ex 11-1
%
% [totalrng, rhomax, rhomin,tgtlat,tgtlon] = ...
%  satfov ( incl,az, slatgd, slon, salt,tfov,etactr,fovmax );
% ------------------------------------------------------------------------------

function [totalrng, rhomax, rhomin,tgtlat,tgtlon] = ...
         satfov ( incl,az, slatgd, slon, salt,tfov,etactr,fovmax );

        rad2deg    =    180.0/pi;

        % -------------------------  implementation   -----------------
        % ------- find satellite parameters and limiting cases --------
        r       = 1.0  + salt;
        etahoriz= asin(1.0 /r);
        rhohoriz= r*cos(etahoriz);

        % ---------------- find ground range angle --------------------
        fovmax= tfov*0.5  + etactr;
        gamma = pi - asin( r*sin(fovmax) );   % must use larger angle
        rho   = cos( gamma ) + r*cos(fovmax);
        rhomax= asin( rho*sin(fovmax) );

        % -------- for minimum, if the sensor looks off axis ----------
        if ( abs(etactr) > 0.00001  )
            fovmin  = etactr - tfov*0.5;
            gamma   = pi - asin( r*sin(fovmin) );  % use larger
            rho     = cos( gamma ) + r*cos(fovmin);
            rhomin  = asin( rho*sin(fovmin) );
            totalrng= rhomax - rhomin
          else
            % --------------------- nadir pointing --------------------
            fovmin  = 0.0;
            rhomin  = 0.0;
            totalrng= 2.0 *rhomax  % equal sided
          end

        % -------------- find location of center of fov ---------------
        if ( abs(etactr) > 0.00001  )
            [lat, lon] = pathm( slatgd,slon,rhomin + totalrng*0.5 ,az )
          else
            lat= slatgd
            lon= slon
          end

        % ----- loop around the new circle with the sensor range ------
        for i= 0 , 72
            az= i*5.0 /rad2deg;
            [tgtlat,tgtlon] = pathm( lat,lon,totalrng*0.5 ,az );
            if ( i == 0 )
                maxlat= tgtlat;
              end
            if ( i == 36 )
                minlat= tgtlat;
              end
          end

