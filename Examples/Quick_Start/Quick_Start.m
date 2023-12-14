% Quick Start
clc;clear;

% add folders and subfolders: supporting_functions and supporting_data
addpath(genpath('../../supporting_data/')); 
addpath(genpath('../../supporting_functions'));

% initial condition file
ICfile = '2020.mat'

% MOCAT MC configuration
seed = 1;% random number generator seed

disp('MC configuration starting...');
cfgMC = setup_MCconfig(seed,ICfile);
fprintf('Seed %i\n', seed);

% MOCAT MC evolution
fprintf('Initial Population:  %i sats\n', size(cfgMC.mat_sats,1));
fprintf('Launches per year: %i\n', size(cfgMC.repeatLaunches,1));
disp('Starting main_mc...');
[nS,nD,nN,nB]=main_mc(cfgMC,seed);

% MOCAT MC postprocess: ratio of satellite in all space objects
ratio = nS/(nS+nD+nN+nB);
fprintf('Quick Start under no launch scenario done!\n')
fprintf('Satellite ratio in all space objects after evolution: %f\n', ratio)
