%     -----------------------------------------------------------------
%
%                              Ex5_1.m
%
%  this file demonstrates example 5-1.
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
%             7 jun 07  david vallado
%                         original
%  changes :
%            13 feb 07  david vallado
%                         original baseline
%
%     *****************************************************************

        constmath;

        jd = jday(2006,4,2,  0,0,0.0)
        fprintf(1,'jd  %11.9f \n',jd );

        [rsun,rtasc,decl] = sun ( jd );
        fprintf(1,'sun  rtasc %14.6f deg decl %14.6f deg\n',rtasc*rad,decl*rad );
        fprintf(1,'sun MOD %11.9f%11.9f%11.9f au\n',rsun );
        fprintf(1,'sun MOD %14.4f%14.4f%14.4f km\n',rsun*149597870.0 );

        rsunaa = [0.9776872 0.1911521  0.0828717]*149597870.0; % astronomical alm value into km
        fprintf(1,'rs aa ICRF %11.9f %11.9f %11.9f km \n',rsunaa);

        da = rsun*149597870.0-rsunaa;
        fprintf(1,'delta mod,eci  %11.9f %11.9f %11.9f %11.4f km\n',da, mag(da) );

        ttt= ( jd - 2451545.0  )/ 36525.0;
        vmod = [0 0 0]';
        amod = [0 0 0]';
        [reci,veci,aeci] = mod2eci  ( rsun',vmod,amod,ttt );

        db = reci*149597870.0-rsunaa';
        fprintf(1,'delta ~eci,eci  %11.9f %11.9f %11.9f %11.4f km %11.4f au \n',db, mag(db),mag(db)/149597870.0 );

        [hms] = hms2rad( 0,44,33.42 );
        [dms] = dms2rad( 4,47,18.3 );
        fprintf(1,'hms ast alm rtasc %11.9f decl %11.9f \n',hms*rad,dms*rad );

