%Demonstrate Calculation of Beacon Satellite Propagation Geometry
%
%    Calculates the propagation geometry and magnetic field
%    parameters from spg10sys orbit and IRGF11 B-field predictions.  For 
%    satellites other than object 32765, tsince parameters and sampling
%    intervals should be adjusted
%
%    The outputs are stored in a .mat file for subsequent display compuation
%
%Required Libraries: GPS_CoordinateXforms, SPG4, IGRF, utilities 
%

%Written by
%Charles L. Rino
%Rino Consulting
%crino@mindspring.com
%Version: September 28, 2010

dtr=pi/180;
min_per_day=60*24;

%***************Get TLE Orbital Elements**********************
path2data='..\SGP4\DemoData';
TLEfiles=dir([path2data,'\*.mat']);
if isempty(TLEfiles)
    error('No TLE file')
else
    %DemoData currently has only a single set of elements
    %ADD library and selection options here for more general use
    load(fullfile(path2data,TLEfiles(2).name));
    longstr1=line1;
    longstr2=line2;
    fprintf('USING 2-Line Elements: \n')
    fprintf('%s \n',longstr1)
    fprintf('%s \n',longstr2)
    load(fullfile(path2data,TLEfiles(1).name));  %
end

%****************Get Station Location**************************
stationfiles=dir([path2data,'\*.mat']);
if isempty(stationfiles)
   error('No station files file')
else
   load([path2data,'\',stationfiles(1).name]);
   origin_llh=[lct.station.rx_latitude*dtr;...
               lct.station.rx_longitude*dtr;...
               lct.station.rx_altitude];
   rx_name=lct.station.rx_name;
end

%********************Initialize spg4****************************
%set in sgp4init
%global tumin mu radiusearthkm xke j2 j3 j4   
%global opsmode
satrec = twoline2rvMOD(longstr1,longstr2);
rx_name=sscanf(rx_name,'%s');
fprintf('\n')
fprintf('Satellite ID %5i \n',satrec.satnum)
fprintf('Station %s: Lon=%6.4f Lat %6.4f Alt=%6.2f m \n',...
            rx_name,origin_llh(2)/dtr,origin_llh(1)/dtr,origin_llh(3))

%********************Run SPG4 for Multiple Orbits****************            
%NOTE:  dt, npts, and tsince offset tailored to object 32765 parameters 
dt=10/60;  %10 sec
npts=5000;
tsince_offset=10000;
%Expand dt, npts, and tsince_offset for broader selection capabilities
tsince=tsince_offset+[0:npts-1]*dt;
if (satrec.epochyr < 57)
    Eyear= satrec.epochyr + 2000;
else
    Eyear= satrec.epochyr + 1900;
end
[Emon,Eday,Ehr,Emin,Esec] = days2mdh(Eyear,0);
UTsec=Ehr*3600+Emin*60+Esec;
fprintf('EPOCH: YEAR MO  DAY UTSEC \n')
fprintf('      %5i %2i %4i %5.2f \n\n',Eyear,Emon,Eday,UTsec);

xsat_ecf=zeros(3,npts);
vsat_ecf=zeros(3,npts);
for n=1:npts
   [satrec, xsat_ecf(:,n), vsat_ecf(:,n)]=spg4_ecf(satrec,tsince(n));
end
%Scale state vectors to mks units
xsat_ecf=xsat_ecf*1000;  %m
vsat_ecf=vsat_ecf*1000;  %mps

sat_llh=ecf2llhT(xsat_ecf);            %ECF to geodetic (llh)  
sat_tcs=llh2tcsT(sat_llh,origin_llh);  %llh to tcs at origin_llh
sat_elev=atan2(sat_tcs(3,:),sqrt(sat_tcs(1,:).^2+sat_tcs(2,:).^2));
%Identify visible segments: 
notVIS=find(sat_tcs(3,:)<0);
VIS=setdiff([1:npts],notVIS);
sat_llh(:,notVIS)=NaN;
sat_tcs(:,notVIS)=NaN;

%Extract start and end times for visible segment: 
t_start=tsince(notVIS(diff(notVIS)>1));
t_end  =tsince(VIS(diff(VIS)>1));
if length(t_end)<length(t_start)
   t_end  =[t_end,tsince(VIS(end))];
elseif length(t_start)<length(t_end)
   t_start=[tsince(notVIS(1)),t_end];
end
t_mid  =(t_start+t_end)/2;

figure
plot(tsince(:)/min_per_day,sat_elev(:)/dtr,'r')
hold on
plot(t_start/min_per_day,zeros(length(t_start)),'b^')
hold on
plot(t_end/min_per_day,zeros(length(t_end)),'b^')
for npass=1:length(t_mid)
    hold on
    npass_str=num2str(npass);
    text(t_end(npass)/min_per_day,10,npass_str)
    [Emon,Eday,Ehr,Emin,Esec] = days2mdh(Eyear,satrec.epochdays+t_start(npass)/min_per_day);
    fprintf('PASS#%3i START: %5i %2i %4i %5.2f ',npass,Emon,Eday,Ehr,Emin);
    [Emon,Eday,Ehr,Emin,Esec] = days2mdh(Eyear,satrec.epochdays+t_end(npass)/min_per_day);
    fprintf('END: %5i %2i %4i %5.2f \n',Emon,Eday,Ehr,Emin);
end
grid on
xlabel('UT--days')
ylabel('elevation--deg')
title(['Satellite ID: ',num2str(satrec.satnum),'--Station: ',rx_name])

%Select Pass for summary
nPass=7;
h_intercept=300000;   %300 km

SatID.satnum=satrec.satnum;
SatID.year=Eyear;
SatID.mon =Emon;
SatID.day =Eday;
plotTitle=['Object ',num2str(satrec.satnum),' Pass ',num2str(nPass)];

t_start=t_start(nPass);
t_end=  t_end(nPass);
npts=ceil((t_end-t_start)*60);
tsince_min=linspace(t_start,t_end,npts);
satGEOM_struct=satGEOM(satrec,Eyear,tsince_min,origin_llh,h_intercept,plotTitle);

SummarizeSatelliteGeometry(satGEOM_struct,SatID,origin_llh,...
                                      tsince_min,plotTitle);
