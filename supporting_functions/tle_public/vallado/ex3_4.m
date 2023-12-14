%     -----------------------------------------------------------------
%
%                              Ex3_4.m
%
%  this file demonstrates example 3-4.
%
%                          companion code for
%             fundamentals of astrodyanmics and applications
%                                 2007
%                            by david vallado
%
%     (h)               email dvallado@msn.com
%     (w) 719-573-2600, email dvallado@agi.com
%
%     *****************************************************************
%
%  current :
%            13 feb 07  david vallado
%                         original
%  changes :
%            13 feb 07  david vallado
%                         original baseline
%
%     *****************************************************************

       % -------- jday         - find julian date
       year = 1996;
       mon = 10;
       day = 26;
       hr  = 14;
       min = 20;
       sec = 0.00;
       fprintf(1,'\n--------jday test \n' );
       fprintf(1,'year %4i ',year);
       fprintf(1,'mon %4i ',mon);
       fprintf(1,'day %3i ',day);
       fprintf(1,'hr %3i:%2i:%8.6f\n ',hr,min,sec );

       jdut1= jday(year,mon,day,hr,min,sec);

       fprintf(1,'jd %18.10f \n',jdut1);


