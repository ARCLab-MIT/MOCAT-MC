%     -----------------------------------------------------------------
%
%                              ex6_5.m
%
%  this file demonstrates example 6-5.
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
%            25 nov 08  david vallado
%                         original
%  changes :
%            25 nov 08  david vallado
%                         original baseline
%
%     *****************************************************************

      fprintf(1,'-------------------- problem ex 6-5 \n');
      rad = 180.0 / pi;
      re = 6378.137;  
      mu = 1.0;  % canonical

      iinit= 55.0 / rad;
      ecc = 0.0;
      deltaomega = 45.0 / rad;
      vinit = 5.892311 / 7.905365;
      fpa = 0.0 / rad;
      incl = 0.0 / rad;

      [ifinal,deltav ] = nodeonly(iinit,ecc,deltaomega,vinit,fpa,incl);

      fprintf(1,'node only changes \n');
      fprintf(1,' ifinal %11.7f \n',ifinal );
      fprintf(1,' deltav %11.7f \n',deltav );


