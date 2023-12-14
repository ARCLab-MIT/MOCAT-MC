%     -----------------------------------------------------------------
%
%                              Ex7_2.m
%
%  this file demonstrates example 7-2.
%
%                          companion code for
%             fundamentals of astrodynamics and applications
%                                 2007
%                            by david vallado
%
%     (w) 719-573-2600, email dvallado@agi.com
%
%     *****************************************************************
%
%  current :
%             9 oct 07  david vallado
%                         original
%  changes :
%             9 oct 07  david vallado
%                         original baseline
%
%     *****************************************************************

   constmath;
   constastro;

   % gaussian angles only (vallado book, p. 420 new effort)
   jd1 = jday(2007, 8, 20, 11, 40, 0.000 );
   jd2 = jday(2007, 8, 20, 11, 50, 0.000 );
%   jd3 = jday(2007, 8, 20, 12,  0, 0.000 );
   jd3 = jday(2007, 8, 20, 12, 20, 0.000 );  % try unequally spaced obs

   % topocentric values
   rtasc1 = -0.4172870/rad;   % right ascension - first sighting
   rtasc2 = 55.0931551/rad;   % right ascension - second sighting
%   rtasc3 = 98.7739537/rad;   % right ascension - third sighting
   rtasc3 = 134.2826693/rad;  % right ascension - third sighting
   decl1 = 17.4626616/rad;    % declination - first sighting
   decl2 = 36.5731946/rad;    % declination - second sighting
%   decl3 = 31.1314513/rad;    % declination - third sighting
   decl3 = 12.0351097/rad;    % declination - third sighting

   % site position
   latgd = 40.0/rad;
   lon   = -110.0/rad;
   alt   = 2.0;                     % km
 
   % at 8-20-07 11:50,
   r2ans = [5963.6422128   5722.1777645   6660.2466242];
   v2ans = [  -4.3643202      4.6055371      1.5157093];
   [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe (r2ans,v2ans);
   fprintf(1,'ans coes %11.4f%11.4f%13.9f%13.7f%11.5f%11.5f%11.5f%11.5f\n',...
           p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad );

   year  =  2007;
   mon   =   8;
   day   =  20;
   hr    =  11;
   min   =  50;
   sec   =   0.0000;
   dut1  =  0.0;
   dat   = 33;
   xp    =  0.0;
   yp    =  0.0;
   lod   =  0.0;
   timezone= 0;
   terms = 0;
   ddpsi = 0.0;
   ddeps = 0.0;

   utc = sec;
   ut1 = utc+dut1;
   tai = utc+dat;
   tt  = tai+32.184;
   jdut1 = jday(year,mon,day,hr,min,ut1);
   jdtt  = jday(year,mon,day,hr,min,tt);
   ttt   =  (jdtt-2451545.0)/36525.0;
   fprintf(1,'year %5i ',year);
   fprintf(1,'mon %4i ',mon);
   fprintf(1,'day %3i ',day);
   fprintf(1,'hr %3i:%2i:%8.6f\n',hr,min,sec );
   fprintf(1,'dut1 %8.6f s',dut1);
   fprintf(1,' dat %3i s',dat);
   fprintf(1,' xp %8.6f "',xp);
   fprintf(1,' yp %8.6f "',yp);
   fprintf(1,' lod %8.6f s\n',lod);

   % -------------- convert each site vector from ecef to eci -----------------
   [rs,vs] = site ( latgd,lon,alt ); % in ecef
   a = [0;0;0];   % dummy acceleration variable for the ecef2eci routine
   [year,mon,day,hr,min,sec] = invjday(jd1);
   [ut1, tut1, jdut1, utc, tai, tt, ttt, jdtt, tdb, ttdb, jdtdb ] ...
             = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );
   [rsite1,vseci,aeci] = ecef2eci(rs,vs,a,ttt,jdut1,lod,xp,yp,2,ddpsi,ddeps);

   [year,mon,day,hr,min,sec] = invjday(jd2);
   [ut1, tut1, jdut1, utc, tai, tt, ttt, jdtt, tdb, ttdb, jdtdb ] ...
             = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );
   [rsite2,vseci,aeci] = ecef2eci(rs,vs,a,ttt,jdut1,lod,xp,yp,2,ddpsi,ddeps);

   [year,mon,day,hr,min,sec] = invjday(jd3);
   [ut1, tut1, jdut1, utc, tai, tt, ttt, jdtt, tdb, ttdb, jdtdb ] ...
             = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );
   [rsite3,vseci,aeci] = ecef2eci(rs,vs,a,ttt,jdut1,lod,xp,yp,2,ddpsi,ddeps);



   [r2,v2] = anglesdr ( decl1,decl2,decl3,rtasc1,rtasc2, ...
                        rtasc3,jd1,jd2,jd3, rsite1,rsite2,rsite3 );
%   [r2,v2] = anglesg ( decl1,decl2,decl3,rtasc1,rtasc2, ...
%                        rtasc3,jd1,jd2,jd3, rsite1,rsite2,rsite3 );


   % -------------- write out answer --------------
   fprintf(1,'v2     %11.7f   %11.7f  %11.7f ',v2);
   fprintf(1,'v2 %11.7f   %11.7f  %11.7f \n',v2/velkmps);
   fprintf(1,'v2 ans %11.7f   %11.7f  %11.7f ',v2ans);
   fprintf(1,'v2 %11.7f   %11.7f  %11.7f \n',v2ans/velkmps);
   fprintf(1,'r2     %11.7f   %11.7f  %11.7f   %11.7f  %11.7f  %11.7f \n',r2/re, r2);
   fprintf(1,'r2 ans %11.7f   %11.7f  %11.7f   %11.7f  %11.7f  %11.7f \n',r2ans/re, r2ans);


   [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe (r2,v2);
   fprintf(1,'          p km       a km      ecc      incl deg     raan deg     argp deg      nu deg      m deg      arglat   truelon    lonper\n');
   fprintf(1,'    coes %11.4f%11.4f%13.9f%13.7f%11.5f%11.5f%11.5f%11.5f%11.5f%11.5f%11.5f\n',...
              p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad, ...
              arglat*rad,truelon*rad,lonper*rad );
       
   [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe (r2ans,v2ans);
   fprintf(1,'ans coes %11.4f%11.4f%13.9f%13.7f%11.5f%11.5f%11.5f%11.5f\n',...
          p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad );


