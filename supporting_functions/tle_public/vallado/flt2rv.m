%
% ----------------------------------------------------------------------------
%
%                           function flt2rv.m
%
%  this function transforms  the flight elements - latgc, lon, fpav, az,
%    position and velocity magnitude into an eci position and velocity vector.
%
%  author        : david vallado                  719-573-2600   17 jun 2002
%
%  revisions
%    vallado     - fix extra terms in rtasc calc                  8 oct 2002
%
%  inputs          description                    range / units
%    rmag        - eci position vector magnitude  km
%    vmag        - eci velocity vector magnitude  km/sec
%    latgc       - geocentric latitude            rad
%    lon         - longitude                      rad
%    fpa         - sat flight path angle          rad
%    az          - sat flight path az             rad
%    ttt         - julian centuries of tt         centuries
%    jdut1       - julian date of ut1             days from 4713 bc
%    lod         - excess length of day           sec
%    xp          - polar motion coefficient       arc sec
%    yp          - polar motion coefficient       arc sec
%    terms       - number of terms for ast calculation 0,2
%    ddpsi,ddeps - corrections for fk5 to gcrf    rad
%
%  outputs       :
%    r           - eci position vector            km
%    v           - eci velocity vector            km/s
%
%  locals        :
%    fpav        - sat flight path anglefrom vert rad
%
%  coupling      :
%    none        -
%
%  references    :
%    vallado       2001, xx
%    chobotov            67
%
% [r,v] = flt2rv ( rmag,vmag,latgc,lon,fpa,az,ttt,jdut1,lod,xp,yp,terms,ddpsi,ddeps );
% ----------------------------------------------------------------------------

function [r,v] = flt2rv ( rmag,vmag,latgc,lon,fpa,az,ttt,jdut1,lod,xp,yp,terms,ddpsi,ddeps );

        twopi = 2.0*pi;

        small        = 0.00000001;

        % -------- form position vector
        recef(1) = rmag*cos(latgc)*cos(lon);
        recef(2) = rmag*cos(latgc)*sin(lon);
        recef(3) = rmag*sin(latgc);
        recef=recef';

        % -------- convert r to eci
        vecef = [0;0;0];
        aecef = [0;0;0];
        [r,v,a] = ecef2eci(recef,vecef,aecef,ttt,jdut1,lod,xp,yp,terms,ddpsi,ddeps);

        % ------------- calculate rtasc and decl ------------------
        temp= sqrt( r(1)*r(1) + r(2)*r(2) );

% v needs to be defined herexxxxxxxxx
        if ( temp < small )
            rtasc= atan2( v(2) , v(1) );
          else
            rtasc= atan2( r(2) , r(1) );
          end
        decl= asin( r(3)/rmag );

        % -------- form velocity vector
        fpav = pi*0.5 - fpa;
        v(1)= vmag*( cos(rtasc)*(-cos(az)*sin(fpav)*sin(decl) + ...
                     cos(fpav)*cos(decl)) - sin(az)*sin(fpav)*sin(rtasc) );
        v(2)= vmag*( sin(rtasc)*(-cos(az)*sin(fpav)*sin(decl) + ...
                     cos(fpav)*cos(decl)) + sin(az)*sin(fpav)*cos(rtasc) );
        v(3)= vmag*( cos(az)*cos(decl)*sin(fpav) + cos(fpav)*sin(decl) );

        r = r';
        v = v';

