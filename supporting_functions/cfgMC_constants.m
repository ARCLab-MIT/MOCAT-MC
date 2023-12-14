% Set conversion units
cfgMC.DAY2MIN = 60*24;
cfgMC.DAY2SEC = cfgMC.DAY2MIN*60;    
cfgMC.YEAR2DAY = 365.2425;  % days(years(1)) https://www.mathworks.com/help/matlab/ref/duration.years.html
cfgMC.YEAR2MIN = cfgMC.YEAR2DAY * cfgMC.DAY2MIN;    
cfgMC.rad = pi/180;

% GLOBAL VARIABLES
whichconst = 84;
[tumin, mu_const, radiusearthkm, xke, j2, j3, j4, j3oj2] = getgravc( whichconst );
cfgMC.tumin = tumin;
cfgMC.mu_const = mu_const;
cfgMC.radiusearthkm = radiusearthkm;
cfgMC.j2 = j2;
cfgMC.omega_earth = 2*pi/(cfgMC.DAY2SEC);