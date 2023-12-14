%     -----------------------------------------------------------------
%
%                              ex3_1415.m
%
%  this file tests the reduction functions.
%
%                          companion code for
%             fundamentals of astrodynamics and applications
%                                 2007
%                            by david vallado
%
%     (h)               email dvallado@msn.com
%     (w) 719-573-2600, email dvallado@agi.com
%
%     *****************************************************************
%
%  current :
%            30 mar 07  david vallado
%                         original
%  changes :
%            13 feb 07  david vallado
%                         original baseline
%
%     *****************************************************************


        recef = [-1033.4793830;  7901.2952754;  6380.3565958;];    
        vecef = [-3.225636520;  -2.872451450;   5.531924446;];
        aecef = [0.001;0.002;0.003];

        year=2004;
        mon = 4;
        day = 6;
        hr =  7;
        min= 51;
        sec= 28.386009;

        dut1 = -0.4399619;
        dat  = 32;
        xp   = -0.140682;  % "
        yp   =  0.333309;
        lod  =  0.0015563;
        ddpsi = -0.052195;  % "
        ddeps = -0.003875;
        timezone=0;
        order = 106;
        eqeterms = 2; % use the extra eqeq terms in j2000
        opt = 'a'; % specify the iau00 approach

        fprintf(1,'test program for reduction functions \n\n');

        fprintf(1,'input data \n\n');
        fprintf(1,' year %5i ',year);
        fprintf(1,' mon %4i ',mon);
        fprintf(1,' day %3i ',day);
        fprintf(1,' %3i:%2i:%8.6f\n ',hr,min,sec );
        fprintf(1,' dut1 %8.6f s',dut1);
        fprintf(1,' dat %3i s',dat);
        fprintf(1,' xp %8.6f "',xp);
        fprintf(1,' yp %8.6f "',yp);
        fprintf(1,' lod %8.6f s\n',lod);
        fprintf(1,' ddpsi %8.6f " ddeps  %8.6f\n',ddpsi,ddeps);
        fprintf(1,' order %3i  eqeterms %31  opt %3s \n',order,eqeterms,opt );
        fprintf(1,'units are km and km/s and km/s2\n' );

        ddpsi = ddpsi*pi / (180*3600);  % rad
        ddeps = ddeps*pi / (180*3600);
        timezone=0;

        % -------- convtime    - convert time from utc to all the others
        fprintf(1,'convtime results\n');
        [ut1, tut1, jdut1, utc, tai, tt, ttt, jdtt, tdb, ttdb, jdtdb ] ...
        = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );
        fprintf(1,'ut1 %8.6f tut1 %16.12f jdut1 %18.11f ',ut1,tut1,jdut1 );
        [h,m,s] = sec2hms( ut1 );
        fprintf(1,'hms %3i %3i %8.6f \n',h,m,s);
        fprintf(1,'utc %8.6f ',utc );
        [h,m,s] = sec2hms( utc );
        fprintf(1,'hms %3i %3i %8.6f \n',h,m,s);
        fprintf(1,'tai %8.6f',tai );
        [h,m,s] = sec2hms( tai );
        fprintf(1,'hms %3i %3i %8.6f \n',h,m,s);
        fprintf(1,'tt  %8.6f ttt  %16.12f jdtt  %18.11f ',tt,ttt,jdtt );
        [h,m,s] = sec2hms( tt );
        fprintf(1,'hms %3i %3i %8.6f \n',h,m,s);
        fprintf(1,'tdb %8.6f ttdb %16.12f jdtdb %18.11f\n',tdb,ttdb,jdtdb );

        % -------- precess     - transformation matrix for precession
        fprintf(1,'precession matrix \n');
        [prec,psia,wa,ea,xa] = precess ( ttt, '80' );
        fprintf(1,'%16.11f %16.11f %16.11f\n',prec );

        % -------- nutation    - transformation matrix for nutation
        fprintf(1,'nutation matrix \n');
        [deltapsi,trueeps,meaneps,omega,nut] = nutation(ttt,0,0);
        fprintf(1,'%16.11f %16.11f %16.11f\n',nut );

        % -------- sidereal    - transformation matrix for sidereal time
        fprintf(1,'sidereal time matrix \n');
        [st,stdot] = sidereal(jdut1,deltapsi,meaneps,omega,lod,2);
        fprintf(1,'%16.11f %16.11f %16.11f\n',st );

        % -------- polarm      - transformation matrix for polar motion
        fprintf(1,'polar motion matrix \n');
        [pm] = polarm(xp,yp,ttt,'80');
        fprintf(1,'%16.11f %16.11f %16.11f\n',pm );

        % -------- truemean    - transformation matrix for teme
        fprintf(1,'truemean matrix \n');
        [nutteme] = truemean ( ttt,order,eqeterms,opt );
        fprintf(1,'%16.11f %16.11f %16.11f\n',nutteme );


        
        
        fprintf(1,'\n\n ============== convert various coordinate systems from ecef =================== \n');
        fprintf(1,'ITRF          IAU-76/FK5   %14.7f %14.7f %14.7f',recef );
        fprintf(1,' v %14.9f %14.9f %14.9f',vecef );
        fprintf(1,' a %14.9f %14.9f %14.9f\n',aecef );

        % -------- pef transformations
        [rpef,vpef,apef] = ecef2pef  ( recef,vecef,aecef, xp, yp, ttt );
        fprintf(1,'PEF           IAU-76/FK5   %14.7f %14.7f %14.7f',rpef );
        fprintf(1,' v %14.9f %14.9f %14.9f',vpef );
        fprintf(1,' a %14.9f %14.9f %14.9f\n',apef );
        
        % -------- tod transformations
        [rtod,vtod,atod] = ecef2tod  ( recef,vecef,aecef,ttt,jdut1,lod,xp,yp,2,ddpsi,ddeps );
        fprintf(1,'TOD 2 w corr  IAU-76/FK5   %14.7f %14.7f %14.7f',rtod );
        fprintf(1,' v %14.9f %14.9f %14.9f',vtod );
        fprintf(1,' a %14.9f %14.9f %14.9f\n',atod );
        [rtod20,vtod20,atod20] = ecef2tod  ( recef,vecef,aecef,ttt,jdut1,lod,xp,yp,2,0,0 );
        fprintf(1,'TOD 2 wo corr IAU-76/FK5   %14.7f %14.7f %14.7f',rtod20 );
        fprintf(1,' v %14.9f %14.9f %14.9f',vtod20 );
        fprintf(1,' a %14.9f %14.9f %14.9f\n',atod20 );
        [rtod00,vtod00,atod00] = ecef2tod  ( recef,vecef,aecef,ttt,jdut1,lod,xp,yp,0,0,0 );
        fprintf(1,'TOD 0 wo corr IAU-76/FK5   %14.7f %14.7f %14.7f',rtod00 );
        fprintf(1,' v %14.9f %14.9f %14.9f',vtod00 );
        fprintf(1,' a %14.9f %14.9f %14.9f\n',atod00 );
        
        % -------- mod transformations
        [rmod,vmod,amod] = ecef2mod  ( recef,vecef,aecef,ttt,jdut1,lod,xp,yp,2,ddpsi,ddeps );
        fprintf(1,'MOD 2 w corr  IAU-76/FK5   %14.7f %14.7f %14.7f',rmod );
        fprintf(1,' v %14.9f %14.9f %14.9f',vmod );
        fprintf(1,' a %14.9f %14.9f %14.9f\n',amod );
        [rmod20,vmod20,amod20] = ecef2mod  ( recef,vecef,aecef,ttt,jdut1,lod,xp,yp,2,0,0 );
        fprintf(1,'MOD 2 wo corr IAU-76/FK5   %14.7f %14.7f %14.7f',rmod20 );
        fprintf(1,' v %14.9f %14.9f %14.9f',vmod20 );
        fprintf(1,' a %14.9f %14.9f %14.9f\n',amod20 );
        [rmod00,vmod00,amod00] = ecef2mod  ( recef,vecef,aecef,ttt,jdut1,lod,xp,yp,0,0,0 );
        fprintf(1,'MOD 0 wo corr IAU-76/FK5   %14.7f %14.7f %14.7f',rmod00 );
        fprintf(1,' v %14.9f %14.9f %14.9f',vmod00 );
        fprintf(1,' a %14.9f %14.9f %14.9f\n',amod00 );

        % -------- eci transformations
        [recig,vecig,aecig] = ecef2eci(recef,vecef,aecef,ttt,jdut1,lod,xp,yp,2,ddpsi,ddeps);
        fprintf(1,'GCRF 2 w corr IAU-76/FK5   %14.7f %14.7f %14.7f',recig );
        fprintf(1,' v %14.9f %14.9f %14.9f',vecig );
        fprintf(1,' a %14.9f %14.9f %14.9f\n',aecig );
        % now do check back to itrf to make sure you get the same vector
        [recef2w,vecef2w,aecef2w] = eci2ecef(recig,vecig,aecig,ttt,jdut1,lod,xp,yp,2,ddpsi,ddeps);
        fprintf(1,'eci-ecef 2, wc             %14.7f %14.7f %14.7f',recef2w );
        fprintf(1,' v %14.9f %14.9f %14.9f',vecef2w );
        fprintf(1,' a %14.9f %14.9f %14.9f\n',aecef2w );

        [reci20,veci20,aeci20] = ecef2eci(recef,vecef,aecef,ttt,jdut1,lod,xp,yp,2,0,0);
        fprintf(1,'GCRF 2wo corr IAU-76/FK5   %14.7f %14.7f %14.7f',reci20 );
        fprintf(1,' v %14.9f %14.9f %14.9f',veci20 );
        fprintf(1,' a %14.9f %14.9f %14.9f\n',aeci20 );
        [recef,vecef,aecef] = eci2ecef(reci20,veci20,aeci20,ttt,jdut1,lod,xp,yp,2,0,0);
        fprintf(1,'ecef 2 wo corr             %14.7f %14.7f %14.7f',recef );
        fprintf(1,' v %14.9f %14.9f %14.9f',vecef );
        fprintf(1,' a %14.9f %14.9f %14.9f\n',aecef );
        
        [reci00,veci00,aeci00] = ecef2eci(recef,vecef,aecef,ttt,jdut1,lod,xp,yp,0,0,0);
        fprintf(1,'GCRF 0wo corr IAU-76/FK5   %14.7f %14.7f %14.7f',reci00 );
        fprintf(1,' v %14.9f %14.9f %14.9f',veci00 );
        fprintf(1,' a %14.9f %14.9f %14.9f\n',aeci00 );
        [recef,vecef,aecef] = eci2ecef(reci00,veci00,aeci00,ttt,jdut1,lod,xp,yp,0,0,0);
        fprintf(1,'ecef 0 wo corr             %14.7f %14.7f %14.7f',recef );
        fprintf(1,' v %14.9f %14.9f %14.9f',vecef );
        fprintf(1,' a %14.9f %14.9f %14.9f\n',aecef );

        
        
        % -------- ecef2teme    - transform ecef to teme vectors
        [rteme,vteme,ateme] = ecef2teme  ( recef,vecef,aecef,ttt,jdut1,lod,xp,yp );
        fprintf(1,'TEME                       %14.7f %14.7f %14.7f',rteme );
        fprintf(1,' v %14.9f %14.9f %14.9f',vteme );
        fprintf(1,' a %14.9f %14.9f %14.9f\n',ateme );
        [receft,veceft,aeceft] = teme2ecef  ( rteme,vteme,ateme,ttt,jdut1,lod,xp,yp );
        fprintf(1,'teme - ecef                %14.7f %14.7f %14.7f',receft );
        fprintf(1,' v %14.9f %14.9f %14.9f',veceft );
        fprintf(1,' a %14.9f %14.9f %14.9f\n',aeceft );
%        % -------- teme2eci    - transform teme to eci vectors
%        [rj2,vj2,aj2] = teme2eci(rteme,vteme,ateme,ttt,4,0,opt);
%        fprintf(1,'teme-eci       %14.7f %14.7f %14.7f',rj2 );
%        fprintf(1,' v %14.9f %14.9f %14.9f',vj2 );
%        fprintf(1,' a %14.9f %14.9f %14.9f\n',aj2 );


        tm = fk4;
        r1950 = tm' * reci20;
        v1950 = tm' * veci20;
        fprintf(1,'FK4                        %14.7f %14.7f %14.7f',r1950 );
        fprintf(1,' v %14.9f %14.9f %14.9f \n',v1950 );

        reci20i = tm * r1950;
        veci20i = tm * v1950;
        fprintf(1,'FK4                        %14.7f %14.7f %14.7f',reci20i );
        fprintf(1,' v %14.9f %14.9f %14.9f \n',veci20i );
        
        % -------- pef2eci     - transform pef to eci vectors
        [reci,veci,aeci] = pef2eci(rpef,vpef,apef,ttt,jdut1,lod,2,0,0);
        fprintf(1,'pef-eci 20                 %14.7f %14.7f %14.7f',reci );
        fprintf(1,' v %14.9f %14.9f %14.9f',veci );
        fprintf(1,' a %14.9f %14.9f %14.9f\n',aeci );

        % -------- tod2eci     - transform tod to eci vectors
        [reci0,veci0,aeci0] = tod2eci(rtod00,vtod00,atod00,ttt,0,0);
        fprintf(1,'tod-eci wo corr            %14.7f %14.7f %14.7f',reci0 );
        fprintf(1,' v %14.9f %14.9f %14.9f',veci0 );
        fprintf(1,' a %14.9f %14.9f %14.9f\n',aeci0 );

        % -------- mod2eci     - transform mod to eci vectors
        [recig,vecig,aecig] = mod2eci(rmod,vmod,amod,ttt);
        fprintf(1,'mod-eci                    %14.7f %14.7f %14.7f',recig );
        fprintf(1,' v %14.9f %14.9f %14.9f',vecig );
        fprintf(1,' a %14.9f %14.9f %14.9f\n',aecig );

        
        [rtod0,vtod0,atod0] = eci2tod(reci,veci,aeci,ttt,0,0);
        fprintf(1,'eci-TOD wo corr            %14.7f %14.7f %14.7f',rtod0 );
        fprintf(1,' v %14.9f %14.9f %14.9f',vtod0 );
        fprintf(1,' a %14.9f %14.9f %14.9f\n',atod0 );

        [rpef20,vpef20,apef20] = eci2pef(reci20,veci20,aeci20,ttt,jdut1,lod,2,0,0);
        fprintf(1,'eci-PEF 20,                %14.7f %14.7f %14.7f',rpef20 );
        fprintf(1,' v %14.9f %14.9f %14.9f',vpef20 );
        fprintf(1,' a %14.9f %14.9f %14.9f\n',apef20 );

        [rpef00,vpef00,apef00] = eci2pef(reci00,veci00,aeci00,ttt,jdut1,lod,0,0,0);
        fprintf(1,'eci-PEF 00                 %14.7f %14.7f %14.7f',rpef00 );
        fprintf(1,' v %14.9f %14.9f %14.9f',vpef00 );
        fprintf(1,' a %14.9f %14.9f %14.9f\n',apef00 );

        [rtodwc,vtodwc,atodwc] = eci2tod(recig,vecig,aecig,ttt,ddpsi,ddeps);
        fprintf(1,'eci-TOD wc                 %14.7f %14.7f %14.7f',rtodwc );
        fprintf(1,' v %14.9f %14.9f %14.9f',vtodwc );
        fprintf(1,' a %14.9f %14.9f %14.9f\n',atodwc );

        [rpef,vpef,apef] = eci2pef(recig,vecig,aecig,ttt,jdut1,lod,2,ddpsi,ddeps);
        fprintf(1,'eci-pef 2 wc               %14.7f %14.7f %14.7f',rpef );
        fprintf(1,' v %14.9f %14.9f %14.9f',vpef );
        fprintf(1,' a %14.9f %14.9f %14.9f\n',apef );



%        % -------- eci2teme    - transform eci to teme vectors
%        [rteme,vteme,ateme] = eci2teme(reci,veci,aeci,ttt,4,0,opt);
%        fprintf(1,'order   4  terms   0  opt %3s \n',opt );
%        fprintf(1,'eci-teme %14.7f %14.7f %14.7f',rteme );
%        fprintf(1,' v %14.9f %14.9f %14.9f',vteme );
%        fprintf(1,' a %14.9f %14.9f %14.9f\n',ateme );



        % --------------------------- now do iau2000
        % -------- eci2ecef    - transform iau2000 eci to ecef vectors
%    recef=[ -1033.4793830  7901.2952758  6380.3565953];
%    vecef=[    -3.225636520  -2.872451450   5.531924446];
        fprintf(1,'\n\n\n\n now do iau2000 processing., start from the itrf vector: \n');
        fprintf(1,'eci-ecef                   %14.7f %14.7f %14.7f',recef );
        fprintf(1,' v %14.9f %14.9f %14.9f',vecef );
        fprintf(1,' a %14.9f %14.9f %14.9f\n',aecef );

        % -------- ecef2eci    - transform iau2000 ecef to eci vectors
        [recigg,vecigg,aecig] = iau00f2i ( recef,vecef,aecef,ttt,jdut1,lod,xp,yp,'a' );
        fprintf(1,'GCRF          IAU-2006     %14.7f %14.7f %14.7f',recigg );
        fprintf(1,' v %14.9f %14.9f %14.9f',vecigg );
        fprintf(1,' a %14.9f %14.9f %14.9f\n',aecig );

        [recefg,vecefg,aecefg] = iau00i2f  ( recigg,vecigg,aecig,ttt,jdut1,lod,xp,yp,'a' );
        fprintf(1,'itrf                       %14.7f %14.7f %14.7f',recefg );
        fprintf(1,' v %14.9f %14.9f %14.9f',vecefg );
        fprintf(1,' a %14.9f %14.9f %14.9f\n',aecefg );

        % -------- ecef2eci    - transform iau2000 ecef to eci vectors
        [recigg,vecigg,aecig] = iau00f2i ( recef,vecef,aecef,ttt,jdut1,lod,xp,yp,'b' );
        fprintf(1,'GCRF          IAU-2000 B   %14.7f %14.7f %14.7f',recigg );
        fprintf(1,' v %14.9f %14.9f %14.9f',vecigg );
        fprintf(1,' a %14.9f %14.9f %14.9f\n',aecig );

        [recefg,vecefg,aecefg] = iau00i2f ( recigg,vecigg,aecig,ttt,jdut1,lod,xp,yp,'b' );
        fprintf(1,'itrf                       %14.7f %14.7f %14.7f',recefg );
        fprintf(1,' v %14.9f %14.9f %14.9f',vecefg );
        fprintf(1,' a %14.9f %14.9f %14.9f\n',aecefg );

        % -------- ecef2eci    - transform iau2000 ecef to eci vectors
        [recigg,vecigg,aecig] = iau00f2i ( recef,vecef,aecef,ttt,jdut1,lod,xp,yp,'c' );
        fprintf(1,'GCRF          IAU-2006 CIO %14.7f %14.7f %14.7f',recigg );
        fprintf(1,' v %14.9f %14.9f %14.9f',vecigg );
        fprintf(1,' a %14.9f %14.9f %14.9f\n',aecig );

        [recefg,vecefg,aecefg] = iau00i2f  ( recigg,vecigg,aeci,ttt,jdut1,lod,xp,yp,'c' );
        fprintf(1,'itrf                       %14.7f %14.7f %14.7f',recefg );
        fprintf(1,' v %14.9f %14.9f %14.9f',vecefg );
        fprintf(1,' a %14.9f %14.9f %14.9f\n',aecefg );


        fprintf(1,'differences \n\n');
        fprintf(1,' gcrf-fk5  %14.9fm \n',mag(recig-reci)*1000 );
        fprintf(1,' fk5-mod   %14.9fm \n',mag(reci-rmod)*1000 );
        fprintf(1,' mod-tod   %14.9fm \n',mag(rmod-rtod)*1000 );
        fprintf(1,' tod-pef   %14.9fm \n',mag(rtod-rpef)*1000 );
        fprintf(1,' pef-ecef  %14.9fm \n',mag(rpef-recef)*1000 );
        fprintf(1,' ggcrf-fk5 %14.9fm \n',mag(recigg-reci)*1000 );

