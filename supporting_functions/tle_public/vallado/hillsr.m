% ------------------------------------------------------------------------------
%
%                           function hillsr
%
%  this function calculates various position information for hills equations.
%
%  author        : david vallado                  719-573-2600    1 mar 2001
%
%  inputs          description                    range / units
%    r           - initial position vector of int m
%    v           - initial velocity vector of int m / s
%    alt         - altitude of tgt satellite      m
%    dts         - desired time                   s
%
%  outputs       :
%    rinit       - final position vector of int   m
%    vinit       - final velocity vector of int   m / s
%
%  locals        :
%    nt          - angular velocity times time    rad
%    omega       -
%    sinnt       - sine of nt
%    cosnt       - cosine of nt
%    radius      - magnitude of range vector       km
%
%  coupling      :
%
%
%  references    :
%    vallado       2007, 397, alg 47, ex 6-14
%
%function [rint, vint] = hillsr( r,v, alt,dts );
% ------------------------------------------------------------------------------

function [rint, vint] = hillsr( r,v, alt,dts );

        % --------------------  implementation   ----------------------
        constastro;
        constmath;
        
        radius= re + alt; % in km
        omega = sqrt( mu / (radius*radius*radius) );
        nt    = omega*dts;
        cosnt = cos( nt );
        sinnt = sin( nt );

        % --------------- determine new positions  --------------------
        rint(1)= ( v(1)/omega ) * sinnt - ...
                 ( (2.0*v(2)/omega) + 3.0*r(1) ) * cosnt + ...
                 ( (2.0*v(2)/omega) + 4.0*r(1) );
        rint(2)= ( 2.0*v(1)/omega ) * cosnt + ...
                 ( (4.0*v(2)/omega) + 6.0*r(1) ) * sinnt + ...
                 ( r(2) - (2.0*v(1)/omega) ) - ...
                 ( 3.0*v(2) + 6.0*omega*r(1) )*dts;
        rint(3)= r(3)*cosnt + (v(3)/omega)*sinnt;

        % --------------- determine new velocities  -------------------
        vint(1)= v(1)*cosnt + (2.0*v(2)+3.0*omega*r(1))*sinnt;
        vint(2)= -2.0*v(1)*sinnt + (4.0*v(2) ...
                  +6.0*omega*r(1))*cosnt - (3.0*v(2)+6.0*omega*r(1));
        vint(3)= -r(3)*omega*sinnt + v(3)*cosnt;



