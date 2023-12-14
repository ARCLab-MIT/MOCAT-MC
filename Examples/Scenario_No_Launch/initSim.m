% Configuration for MC run
% setup_MCconfig and initSim
% -----------------------------------
% initSim:
% configures initial condition in mat_sats and launches in
% repeatLaunches. repeatLaunches is named this way so as to emphasize the
% repeatition feature throughout the propagation.
% 
% Di Wu, 10/17/2023
%------------------------------------------------

function cfgMCout = initSim(cfg, Simulation, launch_model, ICfile)
%   Input: 
%       launch_model: random, matsat (repeat launch based on matsat), no_launch, data, Somma
%
% Di Wu, 10/18/2023
% ------------------------------------------------
    % parameters
    radiusearthkm = cfg.radiusearthkm;
    getidx; %index script

    % TLEs-based initial population + random generation of the launch rate
    try
        fn = which(ICfile);             % find file for initial condition
        load(fn,'mat_sats');            % load variable mat_sats
        load(fn,'time0');               % load variable time0 
        fprintf('%i satellite entries on %s loaded \nfrom %s\n', size(mat_sats,1),time0,fn)

        % Perigees and Apogees
        a_all = mat_sats(:,idx_a)*radiusearthkm; 
        e_all = mat_sats(:,idx_ecco);
        ap_all = a_all.*(1-e_all);
        aa_all = a_all.*(1+e_all);
    catch
        fprintf('current path: %s\n', pwd); 
        error('initial population (.mat) not found in path');
    end

    % Exclude objects out of the desired altitude_limit_up and altitude_limit_low
    mat_sats(ap_all>(cfg.altitude_limit_up + radiusearthkm) | ap_all<(cfg.altitude_limit_low + radiusearthkm),:) = [];

    % Fill in missing DISCOS data if specified
    switch cfg.fillMassRadius 
        case 0      % don't fill in missing DISCOS data (many objects with 0 radius and/or mass)
            [g1,g2,g3] = getZeroGroups(mat_sats);  % just get indexes
        case 1      % 1: ESA's method -- assume spherical aluminum depending on RCS size (S/M/L) 
            [mat_sats,g1,g2,g3] = fillMassRadiusESA(mat_sats); % ESA method for adding in missing physical parameters
        case 2      % 2: resampling method
            [mat_sats,g1,g2,g3] = fillMassRadiusResample(mat_sats); % Resample method
        otherwise
            error('fillMassRadius must be 0, 1 or 2');
    end
        
    % Mark old payloads as derelict (turn off 'controlled' flag)
    payloadInds = mat_sats(:,idx_objectclass) == 1;
    mat_sats(payloadInds, idx_controlled) = 1;      % set controlled flag on for all payloads
    derelict_threshold_year = time0.Year - cfg.missionlifetime; % threshold based on mission lifetime and initial condition year
    ind_derelicts = jd2date(mat_sats(:, idx_launch_date)) < derelict_threshold_year;  % find derelicts
    mat_sats(ind_derelicts, idx_controlled) = 0;    % set controlled flag off for derelicts

    % Set launch_model
    switch launch_model
        case {'no_launch','no'} % no launch
            repeatLaunches = []; % empty launches
            launchMC_step = zeros(0,23);
            additional_launches = zeros(0,23);
            ind_launch = [];
            ind_launch_add = [];
        case {'random','matsat'} % random launch via poisson distribution or repeat launches between years X and Y (ESA style launch)
            launchMC_step = zeros(0,23);
            additional_launches = zeros(0,23);
            ind_launch = [];
            ind_launch_add = [];
            if strcmp(launch_model,'matsat')
                idx_launch_date = 16; idx_mass = 8;
                jds = mat_sats(:,idx_launch_date);       % grab Julian Dates of matsats
                ind_inlaunchwindow = (jd2date(jds) <= cfg.launchRepeatYrs(2) ...
                    & jd2date(jds) >= cfg.launchRepeatYrs(1));
                %-------------------------------------------------------------
                % create REPEATLUNACHES based on past launch performance 
                % within the launch window specified by launchRepeatYrs(1) and launchRepeatYrs(2)
                repeatLaunches = mat_sats(ind_inlaunchwindow,:);
                %-------------------------------------------------------------
                % constellations
                if exist('cfg.constellationFile','var') && ~isempty(cfg.constellationFile) % with constellation
                    ind_constellation = repeatLaunches(:,idx_mass) == 260 | repeatLaunches(:,idx_mass) == 148;
                    repeatLaunches(ind_constellation,:) = [];
                    fprintf('TLE launch repeat selected \n\t %i objects total over %i years ... \n\t omitting %i constellations from repeatlaunches ...\n\t constellations are launched via cfg.constellation from %s\n',...
                    size(repeatLaunches,1), diff(cfg.launchRepeatYrs)+1, sum(ind_constellation),cfg.constellationFile);
                else % without constellation
                    fprintf('TLE launch repeat selected (%i objects total over %i years)\n',...
                    sum(ind_inlaunchwindow), diff(cfg.launchRepeatYrs)+1);
                end
            end
        case {'data','Somma'}
            if strcmp(launch_model,'data')
                filename = 'Pop_Launches_Interface/SSEM_Results/results_3_5e6_5years_1capacity.mat'; % equilibrium solution
            elseif strcmp(launch_model,'Somma')
                filename = 'Pop_Launches_Interface/distribution_launched_objects_altitude_somma.csv';
            end
            [launchMC_step,additional_launches,ind_launch,ind_launch_add] = ...
                prepare_launch_profile_vec(filename, Simulation, cfg.tsince, nyears,...
                cfg.altitude_limit_low, cfg.altitude_limit_up, launch_frequency, ...
                DAY2MIN,YEAR2MIN, paramSSEM);
        otherwise
            error('launch_model must be [data] [Somma] [random] [matsat]')
    end

    % turn on constellation flag for Starlink and Oneweb
        ind_constellation = mat_sats(:,idx_mass) == 260 | mat_sats(:,idx_mass) == 148;
        mat_sats(ind_constellation,idx_constel) = 1;
        if ~isempty(repeatLaunches) % only when repeatLaunches is not empty
        ind_constellation = repeatLaunches(:,idx_mass) == 260 |repeatLaunches(:,idx_mass) == 148 ;
        repeatLaunches(ind_constellation,idx_constel) = 1;
        end
        if ~isempty(additional_launches)
        ind_constellation = additional_launches(:,idx_mass) == 260 | additional_launches(:,idx_mass) == 148 ;
        additional_launches(ind_constellation,idx_constel) = 1;
        end

    % add in mission lifetime
        mat_sats(mat_sats(:,idx_controlled)==1, idx_missionlife) = cfg.missionlifetime;
        if ~isempty(repeatLaunches) % only when repeatLaunches is not empty
        repeatLaunches(repeatLaunches(:,idx_controlled)==1, idx_missionlife) = cfg.missionlifetime;
        end
        if ~isempty(additional_launches)
        additional_launches(additional_launches(:,idx_controlled)==1, idx_missionlife) = cfg.missionlifetime;
        end

    % add in desired_a
        mat_sats(mat_sats(:,idx_controlled)==1, idx_a_desired) = mat_sats(mat_sats(:,idx_controlled)==1,idx_a);
        if ~isempty(repeatLaunches) % only when repeatLaunches is not empty
        repeatLaunches(repeatLaunches(:,idx_controlled)==1, idx_missionlife) = repeatLaunches(repeatLaunches(:,idx_controlled)==1,idx_a);
        end
        if ~isempty(additional_launches)
        additional_launches(additional_launches(:,idx_controlled)==1, idx_missionlife) = additional_launches(additional_launches(:,idx_controlled)==1,idx_a);
        end

    % recalculate Bstar  
    if cfg.physicalBstar
        % Recalculate Bstar as 1/2 * Cd * A/m * rho0;  rho0 = 0.157e6 kg/m^2/RE
        %   Physical def: C_0 = 1/2 * Cd * A/m * rho_0
        %   TLE's use     C_0 = Bstar / 0.157e6 * rho;  [] done in propagation code
        % see line 89 in analytic_propagation_vec.m
        mat_sats(:, idx_bstar) = 0.5 * 2.2 * mat_sats(:, idx_radius).^2 ./ mat_sats(:, idx_mass) * 0.157;
        if ~isempty(repeatLaunches) % only when repeatLaunches is not empty
        repeatLaunches(:, idx_bstar) = 0.5 * 2.2 * repeatLaunches(:, idx_radius).^2 ./ repeatLaunches(:, idx_mass) * 0.157;    
        end
    end

    % outputs
    cfgMCout = cfg;
    cfgMCout.a_all = a_all;
    cfgMCout.ap_all = ap_all;
    cfgMCout.aa_all = aa_all;
    cfgMCout.mat_sats =  mat_sats;
    cfgMCout.repeatLaunches = repeatLaunches;
    cfgMCout.time0 = time0;
    cfgMCout.launchMC_step = launchMC_step;
    cfgMCout.additional_launches = additional_launches;
    cfgMCout.ind_launch = ind_launch;
    cfgMCout.ind_launch_add = ind_launch_add;
    cfgMCout.launch_model = launch_model;

end