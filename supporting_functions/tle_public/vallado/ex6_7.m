%     -----------------------------------------------------------------
%
%                              ex6_7.m
%
%  this file demonstrates example 6-7.
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
      fprintf(1,'-------------------- problem ex 6-7 \n');
      rinit  = (re + 191.0)/re;
      rfinal = (re + 35780.0)/re;
      einit  = 0.0;
      efinal = 0.0;
      iinit = 28.5/rad;
      ifinal=  0.0/rad;
      deltai = ifinal - iinit;

      [deltai1,deltava,deltavb,gam1,gam2] = combined( rinit,rfinal,einit,efinal,deltai );

      fprintf(1,'combined maneuver \n');
      fprintf(1,' deltava  %11.7f \n',deltava );
      fprintf(1,' deltavb  %11.7f \n',deltavb );
      fprintf(1,' deltai1  %11.7f \n',deltai1 * rad );
      fprintf(1,' gam1  %11.7f \n',gam1 * rad );
      fprintf(1,' gam2  %11.7f \n',gam2 * rad );

