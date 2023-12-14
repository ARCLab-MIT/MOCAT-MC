% Scenario No Launch
clc;clear;

% add folders and subfolders: supporting_functions and supporting_data
addpath(genpath('../../supporting_data/')); 
addpath(genpath('../../supporting_functions'));

% initial condition file
ICfile = '2020.mat'

t = 1:5:361;
deorbit_list_m = zeros(3,730);


for idx = 1:3
    seed = idx;% seed = 10; random number generation base

    disp('MC configuration starting...');
    cfgMC = setup_MCconfig(seed,ICfile);
    fprintf('Seed %i\n', seed);

    fprintf('Initial Population:  %i sats\n', size(cfgMC.mat_sats,1));
    fprintf('Launches per year: %i\n', size(cfgMC.repeatLaunches,1));
    disp('Starting main_mc...');

    [nS,nD,nN,nB,deorbitlist_r]=main_mc(cfgMC,seed);
    deorbit_list_m(idx,:) = deorbitlist_r;
end


figure
plot(linspace(2020,2030,730),deorbit_list_m(1,:),"LineWidth",2)
hold on
plot(linspace(2020,2030,730),deorbit_list_m(2,:),"LineWidth",2)
plot(linspace(2020,2030,730),deorbit_list_m(3,:),"LineWidth",2)

legend("Seed 1", "Seed 2", "Seed 3")
xlabel("Time (Year)")
ylabel("Decayed Objects")
title("Population Evolution")
xlim([2020 2030])
hold off

