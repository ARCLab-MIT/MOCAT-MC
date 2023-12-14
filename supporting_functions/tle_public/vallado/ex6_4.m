%     -----------------------------------------------------------------
%
%                              ex6_4.m
%
%  this file demonstrates example 6-4.
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

      rad = 180.0 / pi;
      re = 6378.137;  
      mu = 1.0;  % canonical

      fprintf(1,'-------------------- problem ex 6-4 \n');
      deltai = 15.0 / rad;
      vinit = 5.892311 / 7.905365;  % er/tu
      fpa = 0.0/rad;

      [deltavionly] = ionlychg(deltai,vinit,fpa);

      fprintf(1,'inclination only changes \n');
      fprintf(1,' deltavionly  %11.7f \n\n',deltavionly );


      deltai = 15.0 / rad;
      e = 0.3;
      p = 17858.7836 / re;
      nu = 330.0 / rad;      
      a = p / (1.0 - e * e);
      r = p / (1.0 + e * cos(nu) );
      vinit  = sqrt( (2.0* mu)/r - (mu/a) );
      fpa = (e * sin(nu)) / (1.0 + e * cos(nu));

      [deltavionly] = ionlychg(deltai,vinit,fpa);

      fprintf(1,'inclination only changes \n');
      fprintf(1,' deltavionly  %11.7f \n\n',deltavionly );

