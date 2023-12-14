% Configuration for MC run
% setup_MCconfig and initSim
% -----------------------------------
% setup_MCconfig:
% configures constants, scenarios parameters, propagation
% time, launch models, propagator selection, collision mode, atmosphere
% data, animation, and output files.
% 
% Di Wu, 10/17/2023
%------------------------------------------------

function cfgMC = setup_MCconfig(rngseed,ICfile)
    % Setting up scenario:
    % Inputs:
    %       rngseed: a random seed
    %       ICfile: initial conditions file containing initial condition matrix 
    %       which follows standard mat_sats (matlab satellite) matrix form
    % idx_a = 1; idx_ecco = 2; idx_inclo = 3; idx_nodeo = 4; idx_argpo = 5; idx_mo = 6; idx_bstar = 7; idx_mass = 8; idx_radius = 9;
    % idx_error = 10; idx_controlled = 11; idx_a_desired = 12; idx_missionlife = 13; idx_constel = 14; idx_date_created = 15; idx_launch_date = 16;
    % idx_r = [17 18 19]; idx_v = [20 21 22]; idx_objectclass = 23; idx_ID = 24;
    % Outputs: cfgMC (configure MC) struct for function main_mc.m
    % Di Wu, 10/17/2023
    %---------------------------------------------------------------------
    % other standard MIT setup:
    % No explosions, No collision avoidance
    % 8 year payload lifetime, 95% PMD
    % Launches: 2001-2009 repeated specified by launchRepeatYrs (~ 75 / yr via Fig 2.14 in ESA 2022)
    %---------------------------------------------------------------------

    % constants for MC configure
    cfgMC_constants;
    
    % SCENARIO PARAMETERS
    cfgMC.PMD = 0.95;                   % POST MISSION DISPOSAL (active sats only); ESA: 55% RB, 40% PL, 90% Const
    cfgMC.alph = 0.01;                  % COLLISION AVOIDANCE failure probability with one sat active
    cfgMC.alph_a = 0;                   % COLLISION AVOIDANCE failure probability with both sat active
    cfgMC.orbtol = 5;                   % orbit control tolerance for controlled satellites [km]
    cfgMC.step_control = 2;             % orbit control tolerance checking timesteps
    cfgMC.P_frag = 0;                   % EXPLOSION PROBABILITY per day of Rocket Body Fragmentation (if P_frag=0, explosions are not considered!)
    cfgMC.P_frag_cutoff = 18;           % EXPLOSION PROBABILITY age at which objects don't explode.
    cfgMC.altitude_limit_low = 200;     % SHELL lower limit of altitude [km] 
    cfgMC.altitude_limit_up = 2000;     % SHELL upper limit of altitude [km]
    cfgMC.missionlifetime = 8;          % PAYLOADS operational life [years]

    % SET PROPAGATION TIMES
    t0_prop = 0;                                    % initial PROPAGATION time [min]
    nyears = 1;                                     % length of PROPAGATION [years]
    tf_prop = cfgMC.YEAR2MIN * nyears;              % length of PROPAGATION [min]
    cfgMC.dt_days = 5;                              % CUBE METHOD and PROPAGATION sampling time [days]
    DeltaT = cfgMC.dt_days*cfgMC.DAY2MIN;           % CUBE METHOD and PROPAGATION sampling time [min]
    cfgMC.tsince = t0_prop:DeltaT:t0_prop+tf_prop;  % PROPAGATION time list
    cfgMC.n_time = length(cfgMC.tsince);            % length of PROPAGATION time list
    
    % LAUNCHES
    Simulation = 'TLE';                     % 'TLE'
    launch_model = 'no_launch';             % random, matsat, no_launch, data, Somma      
    % 0: random launch via poisson distribution (see initSim below)
    % 1: repeat launches between years X and Y (ESA style launch)
    
    cfgMC.launchRepeatYrs = [2018,2022];    % Min/max year of obj to repeatedly launch (inclusive)
                                            % Only used if TLElaunchRepeat == 1
    cfgMC.launchRepeatSmooth = 0;           % [0/1] average out the above so yearly launch rate remains the same

    % PREPARE INITIAL CONDITION POPULATION
    fillin_physical_parameters;             % Fill in missing satellite physical parameters
    
    %-----------------------
    % Initialize INITIAL CONDITION POPULATION and LAUNCHES
    cfgMC = initSim(cfgMC, Simulation, launch_model, ICfile);   
    %-----------------------
    % Initialize SHELL information specified by paramSSEM
    paramSSEM.N_shell = 36;
    paramSSEM.h_min = cfgMC.altitude_limit_low; 
    paramSSEM.h_max = cfgMC.altitude_limit_up;
    paramSSEM.R02 = linspace(paramSSEM.h_min,paramSSEM.h_max,paramSSEM.N_shell+1); % R02 is the boundary list of all shells
    paramSSEM.re = radiusearthkm; % km
    cfgMC.paramSSEM = paramSSEM;
    %------------------------

    % PROPAGATOR
    cfgMC.use_sgp4 = false;             % only 'false' is currently supported

    % COLLISION
    cfgMC.skipCollisions = 0;           % if 1, collision step is skipped in main_mc
    cfgMC.max_frag = inf;
    
    % CUBE METHOD
    cfgMC.CUBE_RES = 50;                    % CUBE METHOD resolution for the size of cube
    cfgMC.collision_alt_limit = 45000;      % Ignoring satellites above 45000km for collision evaluation
    
    % ATMOSPHERE
    fillin_atmosphere;                  % Fill in atmosphere setup

    % Animation 
    cfgMC.animation = 'no';             % yes to live plot of the simulation, no otherwise

    % Save output file
    cfgMC.save_diaryName = '';          % save commandline output text to this output;  to skip saving: ''
    cfgMC.save_output_file = 0;         
        % 0: don't save;
        % 1: save entire workspace
        % 2: sats_info, config, params
        % 3: summary (S_MC, N_MC, ...), config
        % 4: summary and collision stats (frag_info)
        % 5: just collision stats (frag_info)
        % 10: save mat_sats every 'saveMSnTimesteps' timesteps
        %       savevars = {'MCconfig','param','paramSSEM','matsatsperN'};
        % 11: flag10 & frag_info5 (collision stats per timestep)
        %      frag_info5: p1 mass, p2 mass (0.1 kg), alt10km, debAn, debBn [Mx5] uint16
    cfgMC.saveMSnTimesteps = 146;  % every ~2 yrs

    % foldername = 'Simulation_Record';
    % % Check if the folder already exists
    % if exist("foldername", 'dir') == 0     % Folder does not exist, so create it
    %     mkdir(foldername);
    % end
    % filename_save = sprintf('%s/TLEIC_year%i_rand%i.mat',foldername,cfgMC.time0.Year,rngseed);    % filename to save information to
    filename_save = sprintf('TLEIC_year%i_rand%i.mat',cfgMC.time0.Year,rngseed);    % filename to save information to
    cfgMC.filename_save = filename_save;
    cfgMC.n_save_checkpoint = inf; % cfgMC.saveMSnTimesteps * 10; % save the results every n_save_checkpoint step
end













