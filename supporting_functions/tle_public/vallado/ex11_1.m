%     -----------------------------------------------------------------
%
%                              Ex11_1.m
%
%  this file demonstrates example 11-1.
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

        % --------  satfov calculations
        incl = 40.0/rad;
        az   = 40.0/rad;
        slatgd = 50.0 /rad;
        slon   = 40.0 / rad;
        salt   = 800.0;  % km
        tfov   = 25.0 / rad;
        etactr = 0.0 / rad;
        fovmax = 60.0 / rad;

        [totalrng, rhomax, rhomin,tgtlat,tgtlon] = ...
             satfov ( incl,az, slatgd, slon, salt,tfov,etactr,fovmax );



