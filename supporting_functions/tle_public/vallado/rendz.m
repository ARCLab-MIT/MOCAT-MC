% ------------------------------------------------------------------------------
%
%                           procedure rendz
%
%  this procedure calculates parameters for a hohmann transfer rendezvous.
%
%  author        : david vallado                  719-573-2600    1 mar 2001
%
%  inputs          description                    range / units
%    rcs1        - radius of circular orbit int   er
%    rcs2        - radius of circular orbit tgt   er
%    einit       - ecc of first orbit
%    efinal      - ecc of final orbit
%    nuinit      - true anomaly of first orbit    0 or pi rad
%    nufinal     - true anomaly of final orbit    0 or pi rad
%    phasei      - initial phase angle (tgt-int)  +(ahead) or -(behind) rad
%    numrevs     - number of revs to wait
%    ktgt        -
%    kint        -
%
%  outputs       :
%    phasef      - final phase angle              rad
%    waittime    - wait before next intercept opp tu
%    deltav      - change in velocity             er/tu
%
%  locals        :
%    dttutrans   - time of flight of trans orbit  tu
%    atrans      - semimajor axis of trans orbit  er
%    angveltgt   - angular velocity of target     rad / tu
%    angvelint   - angular velocity of int        rad / tu
%    leadang     - lead angle                     rad
%
%  coupling      :
%    power       - raise a base to a power
%
%  references    :
%    vallado       2007, 364, alg 44, alg 45, ex 6-8, ex 6-9
%function [ phasef,waittime,deltav] = rendz(rcs1,rcs3,phasei,einit,efinal,nuinit,nufinal,ktgt,kint);
% ----------------------------------------------------------------------------- }

function [ phasef,waittime,deltav] = rendz(rcs1,rcs3,phasei,einit,efinal,nuinit,nufinal,ktgt,kint);
     twopi   =  6.28318530717959;
     mu = 1.0;  % canonical; 

     atrans    = (rcs1 + rcs3) / 2.0;
     dttutrans = pi * sqrt( atrans * atrans * atrans / mu );
     angvelint = sqrt( mu / (rcs1 * rcs1 * rcs1) );
     angveltgt = sqrt( mu / (rcs3 * rcs3 * rcs3) );
     vint      = sqrt( mu/rcs1 )

     % ---------- check for satellites in the same orbits ----------- }
     if abs( angvelint - angveltgt ) < 0.000001  
         periodtrans= ( ktgt * twopi + phasei ) / angveltgt;
         atrans     = (periodtrans/(twopi * kint)) ^(2.0/3.0)
         rp         = 2.0 * atrans - rcs1;
         if rp < 1.0 
             fprintf(1,' error - the transfer orbit intersects the earth ' );
           end; 
         vtrans  = sqrt( ((2.0 * mu)/rcs1) - (mu/atrans) )
         deltav  = 2.0 * (vtrans-vint);
         waittime= 0.0;
% fprintf(1,'tpha ',periodtrans:11:7,' vint ',vint:11:7,' phi ',phasei*rad:11:7 );}
 phasef= phasei;
 waittime= periodtrans;
 leadang= 0.0;
       else
         leadang = angveltgt * dttutrans;
         phasef  = leadang - pi
         waittime= ( phasef - phasei + 2.0 * pi * ktgt ) / ( angvelint - angveltgt );

         a1  = (rcs1 * (1.0 + einit * cos(nuinit))) / (1.0 - einit * einit );
         a2  = ( rcs1 + rcs3 ) / 2.0;
         a3  = (rcs3 * (1.0 + efinal * cos(nufinal))) / (1.0 - efinal * efinal );
         sme1= -mu / (2.0 * a1);
         sme2= -mu / (2.0 * a2);
         sme3= -mu / (2.0 * a3);
     % -----------------  find delta v at point a  ------------------ }
         vinit = sqrt( 2.0 * ( (mu/rcs1) + sme1 ) );
         vtransa= sqrt( 2.0 * ( (mu/rcs1) + sme2 ) );
         deltava= abs( vtransa - vinit );

     % -----------------  find delta v at point b  ------------------ }
         vfinal = sqrt( 2.0 * ( (mu/rcs3) + sme3 ) );
         vtransb= sqrt( 2.0 * ( (mu/rcs3) + sme2 ) );
         deltavb= abs( vfinal - vtransb );
         deltav= deltava + deltavb;
       end;


