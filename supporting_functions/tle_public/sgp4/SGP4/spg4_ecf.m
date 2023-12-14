function  [satrec, xsat_ecf, vsat_ecf, gst]=spg4_ecf(satrec,tsince);
%USAGE:   [satrec, xsat_ecf, vsat_ecf, gst]=spg4_ecf(satrec,tsince);
%     Convert spg4 eci output to ecf

[satrec, xsat_eci, vsat_eci] = sgp4(satrec,tsince);
%This segment converts eci coordinates to ecf
%Compute Greenwich Apparent Siderial Time
gst=gstime(satrec.jdsatepoch+tsince/1440);
%Now rotate the coordinates
CGAST = cos(gst); SGAST = sin(gst);
xsat_ecf(1)= xsat_eci(1)*CGAST+xsat_eci(2)*SGAST;
xsat_ecf(2)=-xsat_eci(1)*SGAST+xsat_eci(2)*CGAST;
xsat_ecf(3)= xsat_eci(3);
%Apply rotation to convert velocity vector from ECI to ECEF coordinates
OMEGAE = 7.29211586D-5;  %Earth rotation rate in rad/s
vsat_ecf(1)= vsat_eci(1)*CGAST+vsat_eci(2)*SGAST+OMEGAE*xsat_ecf(2);
vsat_ecf(2)=-vsat_eci(1)*SGAST+vsat_eci(2)*CGAST-OMEGAE*xsat_ecf(1);
vsat_ecf(3)= vsat_eci(3);
return