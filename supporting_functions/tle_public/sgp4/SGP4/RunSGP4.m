%
%Set paths
clear
%set in sgp4init
global tumin mu radiusearthkm xke j2 j3 j4 j3oj2  

global opsmode

dtr=pi/180;
min_per_day=60*24;
path2TLE='c:\users\workplace\mfiles\ScintillationBook_mfiles\SPG4_MatLab_Code\DemoData';


%***************Get TLE Orbital Elements**********************
if 0
    TLEfiles=dir([path2TLE,'\*.TLE']);
    if isempty(TLEfiles)
        error('No TLE file')
    else
        fid=fopen(TLEfiles(1).name,'r');
        TLE=fscanf(fid,'%69c%69c');
        longstr1=TLE(1:69);
        longstr2=TLE(72:140);
        clear TLE
        fclose(fid);
    end
else
    load([path2TLE,'\SCION-TLE.mat']);
    longstr1=line1;
    longstr2=line2;
end
fprintf('USING 2-Line Elements: \n')
fprintf('%s \n',longstr1)
fprintf('%s \n',longstr2)
  

%****************Get Station Location**************************
stationfiles=dir([path2TLE,'\*.mat']);
if isempty(stationfiles)
   error('No stationfiles file')
else
   load(fullfile(path2TLE,stationfiles(1).name));
   origin_llh=[lct.station.rx_latitude*dtr;...
               lct.station.rx_longitude*dtr;...
               lct.station.rx_altitude];
   rx_name=lct.station.rx_name;
end

%********************Initialize spg4****************************
satrec = twoline2rvMOD(longstr1,longstr2);

fprintf('\n')
fprintf('Satellite ID %5i \n',satrec.satnum)

        
%********************Run SPG4 for Multiple Orbits****************        
    
dt=2/60;  %10 sec
npts=ceil(101/dt);
tsince=[0:npts-1]*dt;
if (satrec.epochyr < 57)
    Eyear= satrec.epochyr + 2000;
else
    Eyear= satrec.epochyr + 1900;
end
[Emon,Eday,Ehr,Emin,Esec] = days2mdh(Eyear,satrec.epochdays);
UTsec=Ehr*3600+Emin*60+Esec;
%timehack=[Eyear,Emon,Eday,Ehr,Emin,Esec];
%GAST=siderealtime(timehack,0);
gst = gstime(satrec.jdsatepoch);
fprintf(' YEAR MO  DAY UTSEC \n')
fprintf('%5i %2i %4i %5.2f gst=%6.4f rad \n',Eyear,Emon,Eday,UTsec,gst);

gst=zeros(1,npts);
xsat_eci=zeros(3,npts);
vsat_eci=zeros(3,npts);
xsat_ecf=zeros(3,npts);
vsat_ecf=zeros(3,npts);
for n=1:npts
    [satrec, xsat_eci(:,n), vsat_eci(:,n)] = sgp4(satrec,tsince(n));
    %This segment converts eci coordinates to ecf
    %Compute Greenwich Apparent Siderial Time
    gst(n)=gstime(satrec.jdsatepoch+tsince(n)/min_per_day);
    %Now rotate the coordinates
    CGAST = cos(gst(n)); SGAST = sin(gst(n));
    xsat_ecf(1,n)= xsat_eci(1,n)*CGAST+xsat_eci(2,n)*SGAST;
    xsat_ecf(2,n)=-xsat_eci(1,n)*SGAST+xsat_eci(2,n)*CGAST;
    xsat_ecf(3,n)= xsat_eci(3,n);
    %Apply rotation to convert velocity vector from ECI to ECEF coordinates
    OMEGAE = 7.29211586D-5;  %Earth rotation rate in rad/s
    vsat_ecf(1,n)= vsat_eci(1,n)*CGAST+vsat_eci(2,n)*SGAST+OMEGAE*xsat_ecf(2,n);
    vsat_ecf(2,n)=-vsat_eci(1,n)*SGAST+vsat_eci(2,n)*CGAST-OMEGAE*xsat_ecf(1,n);
    vsat_ecf(3,n)= vsat_eci(3,n);
end

%Scale state vectors to mks units
xsat_ecf=xsat_ecf*1000;  %m
vsat_ecf=vsat_ecf*1000;  %mps

if 1
iPos=find(xsat_ecf(3,:)>=0);
iNeg=find(xsat_ecf(3,:) <0);
figure
plot3(xsat_ecf(1,iPos),xsat_ecf(2,iPos),xsat_ecf(3,iPos),'r.')
hold on
plot3(xsat_ecf(1,iNeg),xsat_ecf(2,iNeg),xsat_ecf(3,iNeg),'b.')
hold on
grid on
end

sat_llh=ecf2llhT(xsat_ecf);
sat_tcs=llh2tcsT(sat_llh,origin_llh);
sat_elev=atan2(sat_tcs(3,:),sqrt(sat_tcs(1,:).^2+sat_tcs(2,:).^2));


figure
subplot(2,1,1)
plot(sat_llh(2,:)/dtr,sat_llh(1,:)/dtr,'r.')
hold on
plot(origin_llh(2)/dtr,origin_llh(1)/dtr,'cp')
grid on
ylabel('Latitude--deg')
axis([-200 200 -15 15])

subplot(2,1,2)
plot(sat_llh(2,:)/dtr,sat_elev/dtr,'r.')
grid on
ylabel('Elevation from SaoLouis--deg')
xlabel('Longitude--deg')
axis([-200 200 -100 100])



    

