%Test tcs2height
%  USEFUL CONVERSIONS
%  1 m  =  .3048 ft
%  1 mi (statute)  = 1.607 km
%  1 mi (nautical) = 1.852 km 
rtd  = 180/pi;
dtr  = pi/180;
Re   = 6378166;
ft2m = 0.3048;
%%%%   Dam Neck Radar   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rlat= 36.7788;                        %36 46.73'  N
Rlon=-75.9573;                        %75 57.44'  W
Rheight=22;                           %Approx height of antenna (phase center?)
origin =[Rlat; Rlon; 0]*dtr;
%**********************************************************************************
%Generate a tcs file
nsamp=2000;
x=linspace(0,300000,nsamp);
y=x;
z=-100*ones(size(x));
rng=sqrt(x.^2+y.^2+z.^2);
tcs=[x; y; z;];

%Generate height from surface to tcs (negative => below surface)
height=tcs2heightT(tcs,origin);      %Call with origin => WGS84 ellipsoid
height0=tcs2heightT(tcs);            %Call with no second argument or empty => Spherical Earth
figure
subplot(2,1,1)
plot(rng/1000,height/1000,'r')
hold on
plot(rng/1000,height0/1000,'m')
grid on
ylabel('height--km')
legend('WGS84','Sphere')

subplot(2,1,2)
plot(rng/1000,height-height0,'r')
grid on
ylabel('diff--m')
xlabel('range--km')

