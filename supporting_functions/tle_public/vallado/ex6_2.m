%     -----------------------------------------------------------------
%
%                              ex6_2.m
%
%  this file demonstrates example 6-2.
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
      fprintf(1,'-------------------- problem ex 6-2 \n');
      rinit  = (re + 191.3411)/re;
      rb     = (re + 503873.0)/re;
      rfinal = (re + 376310.0)/re;
      einit  = 0.0;
      efinal = 0.0;
      nuinit = 0.0/rad;
      nufinal= 180.0/rad;

      fprintf(1,'initial position \n');
      fprintf(1,' rinit  %11.7f \n',rinit);
      fprintf(1,' rfinal %11.7f \n',rfinal);
      fprintf(1,' einit   %11.7f \n',einit);
      fprintf(1,' efinal  %11.7f \n',efinal);
      fprintf(1,' nuinit  %11.7f \n',nuinit * rad);
      fprintf(1,' nufinal %11.7f \n',nufinal * rad);

      [deltava,deltavb,deltavc,dttu ] = biellip(rinit,rb,rfinal,einit,efinal,nuinit,nufinal);

      fprintf(1,'bi-elliptic answers \n');
      fprintf(1,' deltava  %11.7f \n',deltava );
      fprintf(1,' deltavb  %11.7f \n',deltavb );
      fprintf(1,' deltavc  %11.7f \n',deltavc );
      fprintf(1,' dttu  %11.7f tu %11.7f sec \n',dttu,dttu*13.44685206374);

