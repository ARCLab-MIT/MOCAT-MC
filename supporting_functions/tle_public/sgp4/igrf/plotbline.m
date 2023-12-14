% PLOTBLINE Plot the one of Earth's magnetic field lines using the IGRF.
% 
% Plots a globe and a single magnetic field lines starting at lat_start and
% lon_start. The line will extend a length of distance km. The plot will
% spin if spin is true and will continue to spin until the user hits
% CTRL+C.
% 
% If the user does not have the Mapping Toolbox, a primitive globe
% consisting just of latitude and longitude lines is drawn with the equator
% and prime meridian thicker lines.

clear;
close all;

font = 'Times New Roman';
axis_font = 12;
title_font = 12;

time = datenum([2000 1 1]);
lat_start = -60; % Geodetic latitude in degrees.
lon_start = 180; % Geodetic longitude in degrees.
alt_start = 0; % Altitude in km.
distance = 90e3; % km.
nsteps = abs(distance)/1;
spin = false;

% Get the magnetic field line points.
[lat, lon, alt] = igrfline(time, lat_start, lon_start, alt_start, ...
    'geod', distance, nsteps);
lat = lat(alt > -1); lon = lon(alt > -1); alt = alt(alt > -1);
lon(lon > 180) = lon(lon > 180) - 360;

% Plot the magnetic field lines.
figure;
hold on;

% If the mapping toolbox is not available, just plot an ellipsoid with
% latitude and longitude lines.
if ~license('test', 'MAP_Toolbox')
    
    % WGS84 parameters.
    a = 6378.137; f = 1/298.257223563;
    b = a*(1 - f); e2 = 1 - (b/a)^2; ep2 = (a/b)^2 - 1;
    
    % Latitude lines of the ellipsoid to plot.
    latitudelines = -90:30:90;
    phi = (0:1:360)*pi/180;
    [LATLINES, PHI] = meshgrid(latitudelines*pi/180, phi);
    RLATLINES = sqrt(2)*a*b./sqrt((b^2 - a^2)*cos(2*LATLINES) + a^2 + b^2);
    [xlat, ylat, zlat] = sph2cart(PHI, LATLINES, RLATLINES);
    requator = sqrt(2)*a*b./sqrt((b^2 - a^2)*ones(size(phi)) + a^2 + b^2);
    [xeq, yeq, zeq] = sph2cart(phi, 0, requator);
    
    % Longitude lines of the ellipsoid to plot.
    longitudelines = 0:30:360;
    theta = (-90:1:90)*pi/180;
    [LONLINES, THETA] = meshgrid(longitudelines*pi/180, theta);
    RLONLINES = sqrt(2)*a*b./sqrt((b^2 - a^2)*cos(2*THETA) + a^2 + b^2);
    [xlon, ylon, zlon] = sph2cart(LONLINES, THETA, RLONLINES);
    rprime = sqrt(2)*a*b./sqrt((b^2 - a^2)*cos(2*theta) + a^2 + b^2);
    [xpm, ypm, zpm] = sph2cart(0, theta, rprime);

    % % If you had a vector of coast lines like coast in the Mapping
    % % toolbox, you could plot that too.
    % c = load('coast');
    % tcoast = c.lat*pi/180;
    % phicoast = c.long*pi/180;
    % rcoast = sqrt(2)*a*b./sqrt((b^2 - a^2)*cos(2*tcoast) + a^2 + b^2);
    % [xcoast, ycoast, zcoast] = sph2cart(phicoast, tcoast, rcoast);
    % plot3(xcoast, ycoast, zcoast, 'Color', [0 0.5 0]);
    
    % Convert lla to xyz.
    [x, y, z] = geod2ecef(lat, lon, alt*1e3); % geod coord
    x = x/1e3; y = y/1e3; z = z/1e3;          % geod coord
    % [x, y, z] = sph2cart(lon*pi/180, lat*pi/180, alt); % geoc coord
    
    % Make the plots.
    plot3(x, y, z, 'r');
    plot3(xlat, ylat, zlat, 'b');
    plot3(xeq, yeq, zeq, 'b', 'LineWidth', 5);
    plot3(xlon, ylon, zlon, 'b');
    plot3(xpm, ypm, zpm, 'b', 'LineWidth', 5);
    axis equal;
    
else
    
    % Just plot a globe.
    load topo;
    axesm('globe', 'Geoid', 6371.2)
    meshm(topo, topolegend); demcmap(topo);
    % [x, y, z] = sph2cart(lon*pi/180, lat*pi/180, alt*1e3); % geoc coord
    % [lat, lon, alt] = ecef2geod(x, y, z); alt = alt/1e3;   % geoc coord
    plot3m(lat, lon, alt, 'r');
    
end

% Set the plot background to black.
set(gcf, 'color', 'k');
axis off;
title(['Magnetic Field Line at ' datestr(time)], 'FontName', font, ...
    'FontSize', title_font, 'Color', 'w');

% Spin the plot indefinitely.
index = 0;
while spin
    view(mod(index, 360), 23.5); % Earth's axis tilts by 23.5 degrees
    drawnow;
    pause(0.1);
    index = index - 5;
end