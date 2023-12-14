%
%    vallado       2007, 352-359, ex 6-7
%function [deltai1,deltava,deltavb,gam1,gam2] = combined( rinit,rfinal,einit,efinal,deltai );
%

function [deltai1,deltava,deltavb,gam1,gam2] = combined( rinit,rfinal,einit,efinal,deltai );

     % --------------------  initialize values   ------------------- }
     mu = 1.0; % cannonical
     fprintf(1,'rinit %11.7f %11.7f  rfinal %11.7f  %11.7f \n',rinit, rinit , rfinal, rfinal  );
     a1 = rinit;  % assume circular orbit
     e2 = einit*einit;
     if abs( e2-1.0 ) > 0.000001 
         a2  = (rinit+rfinal) / 2.0;
         sme2= -1.0 / (2.0*a2);
       else
         a2  = 999999.9;  % undefined for parabolic orbit }
         sme2= 0.0;
       end;
fprintf(1,'a1 %11.7f %11.7f  e2 %11.7f \n',a1, a1 , e2 );
fprintf(1,'a2 %11.7f %11.7f \n',a2, a2  );

     sme1= -1.0 / (2.0*a1);

     % -----------------  find delta v at point a  ----------------- }
     vinit = sqrt( 2.0*( (mu/rinit) + sme1 ) );
     vtransa= sqrt( 2.0*( (mu/rinit) + sme2 ) );
%     fpa2a= atan( ( e2*sin(nu2a) ) / ( 1.0 + e2*cos(nu2a) ) );
%     fpa1 = atan( ( einit*sin(nuinit) ) / ( 1.0 + einit*cos(nuinit) ) );
%     deltava= sqrt( vtransa*vtransa + vinit*vinit - 2.0*vtransa*vinit* ...
%                     ( sin(fpa2a)*sin(fpa1)+cos(fpa2a)*cos(fpa1)*cos(deltai)));

     % -----------------  find delta v at point b  ----------------- }
     vfinal = sqrt( mu/rfinal );
vkmps = 7.905365719014;
fprintf(1,'vinit %11.7f %11.7f  vfinal %11.7f  %11.7f \n',vinit, vinit*vkmps, vfinal, vfinal*vkmps );
     vtransb= sqrt( 2.0*( (mu/rfinal) + sme2 ) );
%     fpa2b= atan( ( e2*sin(nu2b) ) / ( 1.0 + e2*cos(nu2b) ) );
%     fpa3 = atan( ( efinal*sin(nufinal) ) / ( 1.0 + efinal*cos(nufinal) ) );


fprintf(1,'vtransa %11.7f %11.7f  vtransb %11.7f  %11.7f \n',vtransa, vtransa*vkmps, vtransb, vtransb*vkmps );
    
%     deltavb= sqrt( vtransb*vtransb + vfinal*vfinal - 2.0*vtransb*vfinal* ...
                     %( sin(fpa2b)*sin(fpa3)+cos(fpa2b)*cos(fpa3)*cos(deltai)));

     ratio = rfinal/rinit;
     s = 1.0/deltai * atan(sin(deltai)/(ratio^1.5+cos(deltai) ) )            
     deltai1 = s*deltai;
     deltai2 = (1.0-s)*deltai;
      
     deltava= sqrt( vinit^2 + vtransa^2 - 2.0*vinit*vtransa*cos(deltai1) );               
     deltavb= sqrt( vfinal^2 + vtransb^2 - 2.0*vfinal*vtransb*cos(deltai2) );                       
                     
                  
     gam1 = acos( -(vinit^2+deltava^2-vtransa^2 ) / (2.0*vinit*deltava) );
     gam2 = acos( -(vtransb^2+deltavb^2-vfinal^2 ) / (2.0*vtransb*deltavb) );

 

