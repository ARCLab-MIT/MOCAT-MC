%Set path for Satellite Orbit Computation
currentDirectory=pwd;
[path2Codes,~]=fileparts(currentDirectory);
addpath(fullfile(path2Codes,'utilities'));
addpath(fullfile(path2Codes,'SGP4'));
%NOTE:  If path to GPS_CoordinateXforms and IGRF are known
%       replace replace ff's with paths and eliminate pickfile 
ff=uigetdir('','GPS_CoordinateXforms');
addpath(ff);
ff=uigetdir('','IGRF');
addpath(ff);