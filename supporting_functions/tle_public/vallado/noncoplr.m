% -----
%   vallado       2007, 370, alg 46, ex 6-10
% function [ttrans,tphase,dvphase,dvtrans1,dvtrans2,aphase ] = noncoplr(phasei,aint,atgt,ktgt,kint,arglatint,nodeint,truelon,deltai);
%------ }

function [ttrans,tphase,dvphase,dvtrans1,dvtrans2,aphase ] = noncoplr(phasei,aint,atgt,ktgt,kint,arglatint,nodeint,truelon,deltai);
      twopi =  6.28318530717959;
      rad   = 57.29577951308230;
      mu = 1.0;  % cannonical

     angvelint = sqrt( mu / (aint * aint * aint) );
     angveltgt = sqrt( mu / (atgt * atgt * atgt) );
     atrans   = (aint + atgt) * 0.5;
     ttrans = pi * sqrt( (atrans * atrans * atrans) / mu );

     deltatnode = phasei / angvelint;

     lead = angveltgt * ttrans

     omeganode = angveltgt * deltatnode;

     phasenew = nodeint + pi - (truelon + omeganode);
 
     leadnew = pi + phasenew;

     tphase= (leadnew - lead + twopi * ktgt) / angveltgt;

     aphase = (mu * (tphase/(twopi * kint)) ^ 2 ) ^(1.0/3.0)

     % -----------------  find deltav's  ----------------- }
     vint= sqrt(mu/aint)
     vphase= sqrt(2.0 * mu/aint - mu/aphase);
     dvphase= vphase - vint;

     vtrans1= sqrt(2.0 * mu/aint - mu/atrans);
     dvtrans1= vtrans1 - vphase;

     vtrans2= sqrt(2.0 * mu/atgt - mu/atrans);
     vtgt= sqrt(mu/atgt)
     dvtrans2= sqrt(vtgt * vtgt + vtrans2 * vtrans2 - 2.0 * vtgt * vtrans2 * cos(deltai));

     ttotal = deltatnode + ttrans + tphase