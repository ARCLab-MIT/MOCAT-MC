% Configuration for MC run

function cfgMC = setupTLE_LNT_historic(rngseed, TLEyear)
    % SCENARIO: Initial population is LEO TLEs from YYYY.mat (>=2000)
    % https://conference.sdo.esoc.esa.int/proceedings/sdc6/paper/199
    
    % Start with current TLE's (18k objects)
    % LNT's LC is an input variable to main_mc_LNT.m

    if nargin < 1
        rngseed = 1;
    end

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

    %%%%%%%%%%% above to be moved out (but also needed below)
    
    % SCENARIO PARAMETERS
    cfgMC.PMD = 0.4;                    % post mission disposal (active sats only); ESA: 55% RB, 40% PL, 90% Const
    cfgMC.PMDconstel = 0.9;             % in orbcontrol_constel.m
    cfgMC.PMDrb = 0.55;                 % in main_mc_LNT
    cfgMC.alph = 0.01;                 % active control failure probability with one sat active
    cfgMC.rad2pLNTfunc = @(x) 1./(1+exp(-25*(x - 0.3))); % curve for alph effectiveness for LNT objects; input: radius
            % rad2pLNTfunc = @(x) 1./(1+exp(-25*(x - 0.3)));
            % plot([0:0.01:1], (1-rad2pLNTfunc([0:0.01:1]))*(1-alph)+alph); ylim([0,1])
    cfgMC.alph_a = 0;                   % active control failure probability with both active
    cfgMC.orbtol = 5;                 % orbital tolerance for controlled satellites [km]
    cfgMC.step_control = 2;             % timesteps to check for orbit control tolerance
    cfgMC.P_frag = 0;                   % probability per day of Rocket Body Fragmentation (if P_frag=0, explosions are not considered!)
    cfgMC.P_frag_cutoff = 18;          % age at which objects don't explode.
    cfgMC.altitude_limit_low = 200;     % lower limit of altitude [km] 
    cfgMC.altitude_limit_up = 2000;     % upper limit of altitude [km] % IADC: perigee < 2000
    
    cfgMC.missionlifetime = 8;          % [years] operational life for payloads

    
    % SET TIMES
    cfgMC.time0 = datetime(TLEyear,1,1);
    t0_prop = 0;                            % [min]
    nyears = 200;                           % [years]
    tf_prop = cfgMC.YEAR2MIN * nyears;      % [min]
    cfgMC.dt_days = 5;                      % [days] sampling time for CUBE method and propagation
    DeltaT = cfgMC.dt_days*cfgMC.DAY2MIN;     % [min] sampling time for CUBE method and propagation
    cfgMC.tsince = t0_prop:DeltaT:t0_prop+tf_prop; % [min] set the length of the simulation, 200 years is typical;
    cfgMC.n_time = length(cfgMC.tsince);
    
    % LAUNCHES
    Simulation = 'TLE';                     % 'TLE'
    cfgMC.total_launch_per_year = 0;        % only for 'TLE' with TLElaunchRepeat = 0 for 'random' launch_model; used in launches_current_step_vec.m 'random' and future_traffic_model_vec; 
    cfgMC.launch_increase_per_year = 0;     % increase in launch rate per year since t0
    launch_frequency = 13;                  % [days]   frequency > DeltaT; only for MOCAT or SOMMA
    TLElaunchRepeat = 1;                    % 0: 'random' launch via poisson distribution (see initSim below)
                                            % 1: ESA style launch (repeat  launches between years X and Y) ==> matsats
    cfgMC.launchRepeatYrs = [2001,2000];    % Min/max year of obj to repeatedly launch (inclusive)
                                            % Only used if TLElaunchRepeat == 1
    cfgMC.launchRepeatSmooth = 0;           % [0/1] average out the above so yearly launch rate remains the same
    % cfgMC.launchMinPerShell = 10;           % [integer >= 0] min number of launches (payload) per shell as specified in paramSSEM.R02

        % Constellation
        cfgMC.constellationFile = ''; 
%         cfgMC.constellationSheet = 'forexport';  % all 54 constellations
        cfgMC.constellationSheet = '';  % all 54 constellations
        cfgMC = setupConstellation(cfgMC);  % populate cfgMC.constellation (table)

    % Modify initial population
        % Fill in missing satellite physical parameters
        cfgMC.fillMassRadius = 2;     % 0: don't fill in missing DISCOS data (many objects with 0 radius and/or mass)
                                      % 1: ESA's method -- assume spherical aluminum depending on RCS size (S/M/L) 
                                      % 2: resampling method
        cfgMC.initpopMultiplier = 1;  % multiplier for initial population 
        cfgMC.physicalBstar = 1;        % [0/1] recalculate B* as Bstar = 1/2 * Cd * A/m * 0.157e6
                                        % calc in initSim; from setupTLEssem_uniform.m

    cfgMC = initSim(cfgMC, Simulation, TLElaunchRepeat);   % Set population and launch rate
    

    % PROPAGATOR
    cfgMC.use_sgp4 = false;             % only 'false' is currently supported

    % COLLISION
    cfgMC.skipCollisions = 0;           % if 1, collision step is skipped in main_mc
    cfgMC.max_frag = inf;
    
    % Cube method Resolution level, note the statement below from J.C.Liou 2003
    % "a dimension of 1% or less of the average semimajor axis of objects in the system is sufficient."
    % 0.5% is used by "Assessing Collision Algorithms for the NewSpace Era"
%     cfgMC.CUBE_RES = mean(cfgMC.a_all(cfgMC.ap_all<(cfgMC.altitude_limit_up+radiusearthkm),1))*0.01; 
    cfgMC.CUBE_RES = 50;  % km
    cfgMC.collision_alt_limit = 45000; %Ignoring satellites above 45000km for collision evaluation
    
    % ATMOSPHERIC MODEL (pre-computed JB2008)
    density_profile = 'JB2008'; % The options are 'static' or 'JB2008'
    cfgMC.density_profile = density_profile;
    if strcmpi(cfgMC.density_profile,'JB2008')
        cfgMC = initJB2008(cfgMC);
    end

    % Animation
    cfgMC.animation = 'no';             % yes to live plot of the simulation, no otherwise

    % Save output file
    cfgMC.save_diaryName = '';          % save commandline output text to this output;  to skip saving: ''
    cfgMC.save_output_file = 11;         
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

    if cfgMC.save_output_file >= 3 && strcmp(Simulation,'TLE')
        % paramSSEM structure needed to summarize data into bins for TLE simulation
        % This needed for Fast_MC2SSEM_population.m
        % Below snippet from prepare_pop_in.m
        paramSSEM.N_shell = 36;
        paramSSEM.h_min = cfgMC.altitude_limit_low;  % [km] defined in setup above
        paramSSEM.h_max = cfgMC.altitude_limit_up;
        paramSSEM.R02 = linspace(paramSSEM.h_min,paramSSEM.h_max,paramSSEM.N_shell+1);
        paramSSEM.re = radiusearthkm; % km
        cfgMC.paramSSEM = paramSSEM;
    end
    tv = datestr(now, 'yyyymmddTHHMMSS');
    filename_save = sprintf('TLE_%s_%i.mat',tv,rngseed);
    cfgMC.filename_save = filename_save;
    cfgMC.n_save_checkpoint = inf; % cfgMC.saveMSnTimesteps * 10; % save the results every n_save_checkpoint step
end


function cfgMCout = initJB2008(cfgMC)
    % NB: it works just from March 2020 as a starting date of the propagation
    load('dens_jb2008_032020_022224.mat')
    
    dens_times=zeros(length(dens_highvar.month),1);
    
    for k=1:length(dens_highvar.month)
        dens_times(k,1)= juliandate(datetime(dens_highvar.year(k),dens_highvar.month(k),0));
    end
    
    [dens_times2,alt2] = meshgrid(dens_times,dens_highvar.alt);

    cfgMCout = cfgMC;
    cfgMCout.param.dens_times = dens_times2;
    cfgMCout.param.alt = alt2;
    cfgMCout.param.dens_value = dens_highvar.dens;
end

function cfgMCout = initSim(cfg, Simulation, TLElaunchRepeat)
    whichconst = 84;
    [tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2] = getgravc( whichconst );

    % TLEs-based initial population + random generation of the launch rate
    try
%       myVars = {'mat_sats','jdsatepoch_sat1','a_all','ap_all','aa_all','mass_all','radius_all'};
%         fn = which('initialized_01-2023.mat');  % 1 MB
        fn = which(sprintf('%i.mat',year(cfg.time0)));  % 1 MB
                    % old version: initialized.mat; 25 MB            
%       load(fn, myVars{:});
        if isempty(fn)
            error('%i.mat NOT FOUND FOR TLE DATA',year(cfg.time0));
        end
        load(fn);  % just 'mat_sats'
        fprintf('%i satellite entries loaded from %s\n', size(mat_sats,1),fn)

        % fill in other variables from initialized.mat
%       time0 = datetime(jdsatepoch_sat1,'convertfrom','juliandate'); % [datetime]
%       'a_all','ap_all','aa_all','mass_all','radius_all'};
%         time0 = datetime(2021,1,1);
        % compute all Perigees and Apogees
        a_all = mat_sats(:,1)*radiusearthkm;
        e_all = mat_sats(:,2);
        ap_all = a_all.*(1-e_all);
        aa_all = a_all.*(1+e_all);

    catch
        fprintf('current path: %s\n', pwd); 
        error('initial population (.mat) not found in path');
    end

%   time0 = datetime(jdsatepoch_sat1,'convertfrom','juliandate'); % [datetime]
    % Extract objects below the desired altitude
    mat_sats(ap_all>(cfg.altitude_limit_up + radiusearthkm) | ap_all<(cfg.altitude_limit_low + radiusearthkm),:) = [];
    if ~TLElaunchRepeat
        launch_model = 'random'; % generate random launches
    else
        launch_model = 'matsat'; % generate set launch via mat_sats (~ESA)
    end

    % fill in missing DISCOS data if specified
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
    getidx;
    payloadInds = mat_sats(:,idx_objectclass) == 1;
    mat_sats(payloadInds , idx_controlled) = 1;  % set controlled flag on for all payloads...
    ind_matsat_old = jd2date(mat_sats(:, idx_launch_date)) < year(cfg.time0) - cfg.missionlifetime;  % then turn off old payloads    
    mat_sats(ind_matsat_old ,idx_controlled) = 0;

    % multiply initial population
    mat_sats = multiplyInitialPop(cfg,mat_sats,g1,g2,g3);   % only cfg.mat_sat affected; resampling method used

    % Set Launch model
    switch launch_model
        case {'random','matsat'}
            launchMC_step = zeros(0,23);
%             additional_launches = zeros(0,23);
            ind_launch = [];
            ind_launch_add = [];
            if strcmp(launch_model,'matsat')
                idx_launch_date = 16; idx_mass = 8;
                jds = mat_sats(:,idx_launch_date);       % grab Julian Dates of matsats
                launchwindow = (jd2date(jds) <= cfg.launchRepeatYrs(2) ...
                    & jd2date(jds) >= cfg.launchRepeatYrs(1));
                cfg.repeatLaunches = mat_sats(launchwindow,:);
                if ~isempty(cfg.constellationFile)
                    constellationIdx = cfg.repeatLaunches(:,idx_mass) == 260 | cfg.repeatLaunches(:,idx_mass) == 148;
                    cfg.repeatLaunches(constellationIdx,:) = [];
                    fprintf('TLE launch repeat selected \n\t %i objects total over %i years ... \n\t omitting %i constellations from repeatlaunches ...\n\t constellations are launched via cfg.constellation from %s\n',...
                    size(cfg.repeatLaunches,1), diff(cfg.launchRepeatYrs)+1, sum(constellationIdx),cfg.constellationFile);
                else
                    fprintf('TLE launch repeat selected (%i objects total over %i years)\n',...
                    sum(launchwindow), diff(cfg.launchRepeatYrs)+1);
                end
            end
        case {'data','Somma'}
            if strcmp(launch_model,'data')
                filename = 'Pop_Launches_Interface/SSEM_Results/results_3_5e6_5years_1capacity.mat'; % equilibrium solution
                % filename = 'Pop_Launches_Interface/SSEM_Results/launch_rate_Somma.mat'; % Somma launch rate
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
        constellationIdx = mat_sats(:,idx_mass) == 260 | mat_sats(:,idx_mass) == 148;
        mat_sats(constellationIdx,idx_constel) = 1;
        constellationIdx = cfg.repeatLaunches(:,idx_mass) == 260 |cfg.repeatLaunches(:,idx_mass) == 148 ;
        cfg.repeatLaunches(constellationIdx,idx_constel) = 1;

    % add in mission lifetime (!!)
        mat_sats(mat_sats(:,idx_controlled)==1, idx_missionlife) = cfg.missionlifetime;
        cfg.repeatLaunches(cfg.repeatLaunches(:,idx_controlled)==1, idx_missionlife) = cfg.missionlifetime;

    % add in desired_a
        mat_sats(mat_sats(:,idx_controlled)==1, idx_a_desired) = mat_sats(mat_sats(:,idx_controlled)==1,idx_a);
        cfg.repeatLaunches(cfg.repeatLaunches(:,idx_controlled)==1, idx_missionlife) = cfg.repeatLaunches(cfg.repeatLaunches(:,idx_controlled)==1,idx_a);

    % RECALCULATE B*  
    if cfg.physicalBstar
        % Recalculate B* as 1/2 * Cd * A/m * rho0;  rho0 = 0.157e6 kg/m^2/RE
        %   Physical def: C_0 = 1/2 * Cd * A/m * rho_0
        %   TLE's use     C_0 = Bstar / 0.157e6 * rho;  [] done in propagation code
        % see line 89 in analytic_propagation_vec.m
        mat_sats(:, idx_bstar) = 0.5 * 2.2 * mat_sats(:, idx_radius).^2 ./ mat_sats(:, idx_mass) * 0.157;
        cfg.repeatLaunches(:, idx_bstar) = 0.5 * 2.2 * cfg.repeatLaunches(:, idx_radius).^2 ./ cfg.repeatLaunches(:, idx_mass) * 0.157;    
    end

    cfgMCout = cfg;
    
    cfgMCout.a_all = a_all;
    cfgMCout.ap_all = ap_all;
    cfgMCout.aa_all = aa_all;
    cfgMCout.mat_sats =  mat_sats;
    cfgMCout.time0 = cfg.time0;
    cfgMCout.launchMC_step = launchMC_step;
    cfgMCout.ind_launch = ind_launch;
    cfgMCout.ind_launch_add = ind_launch_add;
    cfgMCout.launch_model = launch_model;

end

function [outmatsat,g1,g2,g3] = fillMassRadiusESA(inmatsat)
    % Data pipeline for mat_sats from initialized_01-2023.mat
    % 1) Get TLEs from space-track.org
    % 2) 
    % 3)
    % ESA method for filling in missing info from DISCOS
    % a) use space-track SATCAT designation for RCS size (S/M/L) and assign:
    %    0.1 / 1 / 10 m^2 area, for spherical density of aluminum (2,710 kg/m3)
    SATCATstname = 'spacetrack_satcat_03_2023.csv';
    satcatfn = which(SATCATstname);
    opts = detectImportOptions(SATCATstname); 
    satcatstdata = readtable(SATCATstname,opts);
    fprintf('Using RCS info from: %s \n', satcatfn);
    
    satcatSatnum = satcatstdata.NORAD_CAT_ID;

        % MATSATS DEFINITION
        idx_a = 1; idx_ecco = 2; idx_inclo = 3; 
        idx_nodeo = 4; idx_argpo = 5; idx_mo = 6; 
        idx_bstar = 7; idx_mass = 8; idx_radius = 9;
        idx_error = 10; idx_controlled = 11; idx_a_desired = 12;   % [0,0,0]
        idx_missionlife = 13; idx_constel = 14; idx_date_created = 15;  % [0,0,nan]
        idx_launch_date = 16; idx_r = [17 18 19]; idx_v = [20 21 22]; % date: JD,  RV:0
        idx_objectclass = 23; idx_ID = 24;

    % objects without radius
    noRind = inmatsat(:,idx_radius) == 0;
    noRSatnum = inmatsat(noRind,idx_ID);
    
    [~,b,c] = intersect(noRSatnum,satcatSatnum);
    
    nlg = cellfun(@(x) strcmp(x,'LARGE'),satcatstdata(c,:).RCS_SIZE);
    nmd = cellfun(@(x) strcmp(x,'MEDIUM'),satcatstdata(c,:).RCS_SIZE);
    nsm = cellfun(@(x) strcmp(x,'SMALL'), satcatstdata(c,:).RCS_SIZE);

    fprintf('Objects without radius: %i small, %i medium, %i large objects \n ... now assigned with 0.1,1,10 m^2 area\n',...
        sum(nsm),sum(nmd),sum(nlg));
    fprintf('Summary: %i of %i objects without radius are now assigned \n',...
        sum(nsm + nmd + nlg), sum(noRind));
    
    outmatsat = inmatsat;
    noRinds = find(noRind);
    outmatsat(noRinds(nlg),idx_radius) = sqrt(10/pi);  % 10 m^2
    outmatsat(noRinds(nmd),idx_radius) = sqrt(1/pi);  % 1 m^2
    outmatsat(noRinds(nsm),idx_radius) = sqrt(0.1/pi);  % 0.1 m^2

    % categorize "unknown" objects that are still 0 radius    
%     noRSatnum(logical(nlg + nmd + nsm)) = [];  % find remainding noRsats
    noRind = outmatsat(:,idx_radius) == 0;
    noRSatnum = outmatsat(noRind,idx_ID);
    noRinds = find(noRind);

    [~,~,cc] = intersect(noRSatnum,satcatSatnum);

    % PAYLOAD -> LARGE
    % DEBRIS -> SMALL
    for ci = 1:numel(cc)
        c = cc(ci);
        curType = satcatstdata(c,:).OBJECT_TYPE{1};
        switch curType 
            case 'DEBRIS'
                outmatsat(noRinds(ci),idx_radius) = sqrt(0.1/pi);
            case 'PAYLOAD'
                outmatsat(noRinds(ci),idx_radius) = sqrt(10/pi);
            otherwise
                warning('SAT ID %i has type %s, skipping', satcatSatnum(c), curType);
        end
    end
    
    % ADD MASS (spherical aluminum for ESA; 2,710 kg/m3)
    noMind = outmatsat(:,idx_mass) == 0;
    outmatsat(noMind,idx_mass) = 4/3 * pi * outmatsat(noMind,idx_radius).^3 * 2710; % kg
    fprintf('Added mass to %i objects without mass\n', sum(noMind));

    [g1,g2,g3] = getZeroGroups(inmatsat); % fill out g1,g2,g3
    g1.gm = gmdistribution([sqrt(10/pi), 4/3 * pi * sqrt(10/pi).^3 * 2710],[0,0]); % Payload mu for [r,m]
    g2.gm = gmdistribution([sqrt(10/pi), 4/3 * pi * sqrt(10/pi).^3 * 2710],[0,0]); % RB 
    g3.gm = gmdistribution([sqrt(0.1/pi), 4/3 * pi * sqrt(0.1/pi).^3 * 2710],[0,0]); % Deb

end

function [outmatsat,g1,g2,g3] = fillMassRadiusResample(inmatsat,varargin)
    % Resample existing data to fill in missing info from DISCOS
    % Input: inmatsat -- matsat to be modified
    %        g1,g2,g3 -- optional; gX.gm defines the sampling dist (gmdistribution object)
    % object classes via objclass2int(1:12,2): 
    %   P,PMRO,Pfrag,Pdeb,RB,RBMRO,RBfrag,RBdeb,Deb,OtherDeb,Unkwn,untracked

    % MATSATS DEFINITION
    idx_a = 1; idx_ecco = 2; idx_inclo = 3; 
    idx_nodeo = 4; idx_argpo = 5; idx_mo = 6; 
    idx_bstar = 7; idx_mass = 8; idx_radius = 9;
    idx_error = 10; idx_controlled = 11; idx_a_desired = 12;   % [0,0,0]
    idx_missionlife = 13; idx_constel = 14; idx_date_created = 15;  % [0,0,nan]
    idx_launch_date = 16; idx_r = [17 18 19]; idx_v = [20 21 22]; % date: JD,  RV:0
    idx_objectclass = 23; idx_ID = 24;

% dobj = inmatsat(:,idx_objectclass);
% for stats, combine into:
%    payload {1}, RB {5}, debris {2,3,4,6,7,8,9,10,11}


% fit data, save mean and covariances - if available, use given g1,g2,g3's GM model
if nargin == 1
    [g1,g2,g3] = getZeroGroups(inmatsat);
    X = inmatsat(g1.nzno,[idx_radius, idx_mass]);
    GMModel = fitgmdist(X,1);  g1.gm = GMModel;
    X = inmatsat(g2.nzno,[idx_radius, idx_mass]);
    GMModel = fitgmdist(X,1);  g2.gm = GMModel;
    X = inmatsat(g3.nzno,[idx_radius, idx_mass]);
    GMModel = fitgmdist(X,1);  g3.gm = GMModel;
elseif ~isempty(intersect(fieldnames(varargin{1}),'gm'))
    g1 = varargin{1};  % if given, use provided GM models
    g2 = varargin{2};
    g3 = varargin{3};
% else
%     g1.gm = gmdistribution(0,0);
%     g2.gm = gmdistribution(0,0);
%     g3.gm = gmdistribution(0,0);
end

% fill in empty data from mat_sats via sampling
% -- current method: re-sample even if only one parameter is missing (radius or mass)
% -- method: sample 2x N needed, remove any negative entries, randomly choose N
cursamp = random(g1.gm, numel(union(g1.zm, g1.zr))*2);
cursamp(any(cursamp<=0,2),:) = []; % remove negative entries
outmatsat = inmatsat;
outmatsat(union(g1.zm, g1.zr),[idx_radius,idx_mass]) = cursamp(1:numel(union(g1.zm, g1.zr)),:);

% g3: RB
cursamp = random(g2.gm, numel(union(g2.zm, g2.zr))*2);
cursamp(any(cursamp<=0,2),:) = []; % remove negative entries
outmatsat = inmatsat;
outmatsat(union(g2.zm, g2.zr),[idx_radius,idx_mass]) = cursamp(1:numel(union(g2.zm, g2.zr)),:);

% g3: Debris
cursamp = random(g3.gm, numel(union(g3.zm, g3.zr))*2);
cursamp(any(cursamp<=0,2),:) = []; % remove negative entries
outmatsat = inmatsat;
outmatsat(union(g3.zm, g3.zr),[idx_radius,idx_mass]) = cursamp(1:numel(union(g3.zm, g3.zr)),:);



return;

% figure; % all non-zero objects' mass vs radius
% subplot(3,1,1);  
% plot(inmatsat(g1.nz,idx_radius),inmatsat(g1.nz,idx_mass),'o');
% title('Payload')
% subplot(3,1,2);  
% plot(inmatsat(g2.nz,idx_radius),inmatsat(g2.nz,idx_mass),'o');
% title('RB')
% subplot(3,1,3); 
% plot(inmatsat(g3.nz,idx_radius),inmatsat(g3.nz,idx_mass),'o');
% title('Debris')
% xlabel('Radius (m)'); ylabel('Mass (kg)')
% 
% mass vs radius of NON-OUTLIERS (radius Or mass)
figure;
subplot(3,1,1);  
plot(inmatsat(g1.nzno,idx_radius),inmatsat(g1.nzno,idx_mass),'o');
title(sprintf('Payload (non-outliers: %i)',numel(g1.nzno)));
    % also plot aluminum
    xl = xlim; yl = ylim; alrs = linspace(xl(1),xl(2),100);
    hold on; plot(alrs, 4/3 * pi * alrs.^3 * 2710,'r--'); ylim(yl);
    legend('Payload','Al sphere (2710 kg/m^3)')
subplot(3,1,2);  
plot(inmatsat(g2.nzno,idx_radius),inmatsat(g2.nzno,idx_mass),'o');
title(sprintf('Rocket Body (non-outliers: %i)',numel(g2.nzno)));
    xl = xlim; yl = ylim; alrs = linspace(xl(1),xl(2),100);
    hold on; plot(alrs, 4/3 * pi * alrs.^3 * 2710,'r--'); ylim(yl);
subplot(3,1,3); 
plot(inmatsat(g3.nzno,idx_radius),inmatsat(g3.nzno,idx_mass),'o');
title(sprintf('Debris (non-outliers: %i)',numel(g3.nzno)));
    xl = xlim; yl = ylim; alrs = linspace(xl(1),xl(2),100);
    hold on; plot(alrs, 4/3 * pi * alrs.^3 * 2710,'r--'); ylim(yl);
xlabel('Radius (m)'); ylabel('Mass (kg)')

% mass vs radius of NON-OUTLIERS (vs Alunimum) vs r^3 (m^3)
figure;
subplot(3,1,1);  
plot(4/3 * pi * inmatsat(g1.nzno,idx_radius).^3,inmatsat(g1.nzno,idx_mass),'o');
title(sprintf('Payload (non-outliers: %i)',numel(g1.nzno)));
    % also plot aluminum
    xl = xlim; yl = ylim; alrs = linspace(xl(1),xl(2),100);
    hold on; plot(4/3 * pi * alrs.^3, 4/3 * pi * alrs.^3 * 2710,'r--'); ylim(yl);
    xlim([0,40]); ylim([0,2000]);
%     legend('Payload','Al sphere (2710 kg/m^3)')
subplot(3,1,2);  
plot(4/3 * pi * inmatsat(g2.nzno,idx_radius).^3,inmatsat(g2.nzno,idx_mass),'o');
title(sprintf('Rocket Body (non-outliers: %i)',numel(g2.nzno)));
    xl = xlim; yl = ylim; alrs = linspace(xl(1),xl(2),100);
    hold on; plot(4/3 * pi * alrs.^3, 4/3 * pi * alrs.^3 * 2710,'r--'); ylim(yl);
subplot(3,1,3); 
plot(4/3 * pi * inmatsat(g3.nzno,idx_radius).^3,inmatsat(g3.nzno,idx_mass),'o');
title(sprintf('Debris (non-outliers: %i)',numel(g3.nzno)));
    xl = xlim; yl = ylim; alrs = linspace(xl(1),xl(2),100);
    hold on; plot(4/3 * pi * alrs.^3, 4/3 * pi * alrs.^3 * 2710,'r--'); ylim(yl);
xlabel('Volume (m^3)'); ylabel('Mass (kg)')

% 
% 
% 
% % 2D histogram (WITH outliers)
% figure;
% subplot(311); histogram2(inmatsat(g1.nz,idx_radius), inmatsat(g1.nz,idx_mass), ...
%     linspace(0,6,20),linspace(0,2e4,20), 'displaystyle','bar','ShowEmptyBins','off','FaceColor','flat');
% title(sprintf('Payload (with data: %i)',numel(g1.nz)));
% xlabel('Radius (m)'); ylabel('Mass (kg)')
% subplot(312); histogram2(inmatsat(g2.nz,idx_radius), inmatsat(g2.nz,idx_mass), ...
%     linspace(0,6,20),linspace(0,1e4,20), 'displaystyle','bar','ShowEmptyBins','off','FaceColor','flat');
% title(sprintf('Rocket Body (with data: %i)',numel(g2.nz)));
% xlabel('Radius (m)'); ylabel('Mass (kg)')
% subplot(313); h = histogram2(inmatsat(g3.nz,idx_radius), inmatsat(g3.nz,idx_mass), ...
%     linspace(0,3.5,20),linspace(0,3e3,20),'displaystyle','bar','ShowEmptyBins','off','FaceColor','flat');
% title(sprintf('Debris (with data: %i)',numel(g3.nz)));
% xlabel('Radius (m)'); ylabel('Mass (kg)')
% 
% 
% % 2D histogram (non-outliers)
% figure;
% subplot(311); histogram2(inmatsat(g1.nzno,idx_radius), inmatsat(g1.nzno,idx_mass), ...
%     linspace(0,4,30), linspace(0,1e3,20), 'displaystyle','bar','ShowEmptyBins','off','FaceColor','flat');
% title(sprintf('Payload (non-outliers: %i)',numel(g1.nzno)));
% xlabel('Radius (m)'); ylabel('Mass (kg)')
% 
% subplot(312); histogram2(inmatsat(g2.nzno,idx_radius), inmatsat(g2.nzno,idx_mass), ...
%     linspace(0,4,30), linspace(0,2e3,20), 'displaystyle','bar','ShowEmptyBins','off','FaceColor','flat');
% title(sprintf('Rocket Body (non-outliers: %i)',numel(g2.nzno)));
% xlabel('Radius (m)'); ylabel('Mass (kg)')
% subplot(313); h = histogram2(inmatsat(g3.nzno,idx_radius), inmatsat(g3.nzno,idx_mass), ...
%     linspace(0,4,30), linspace(0,200,20),'displaystyle','bar','ShowEmptyBins','off','FaceColor','flat');
% title(sprintf('Debris (non-outliers: %i)',numel(g3.nzno)));
% xlabel('Radius (m)'); ylabel('Mass (kg)')
% %     xl = xlim; yl = ylim; alrs = linspace(xl(1),xl(2),100);
% %     hold on; plot(alrs, 4/3 * pi * alrs.^3 * 2710,'r--'); ylim(yl);
% 
% 
% % 2d data fitting for radius/mass data
% figure;
% subplot(133); 
% X = inmatsat(g3.nzno,[idx_radius, idx_mass]);
% 
% GMModel = fitgmdist(X,1);
% y = zeros(size(X,1),1);
% plot(X(:,1), X(:,2),'x');  hold on
% plot(mean(X(:,1)), mean(X(:,2)),'ro');
% plot(median(X(:,1)), median(X(:,2)),'rs');
% gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(GMModel,[x0 y0]),x,y);
% g = gca;
% fcontour(gmPDF,[g.XLim g.YLim])
% % title('{\bf Scatter Plot and Fitted Gaussian Mixture Contours}')
% % legend(h,'Model 0','Model1')
% hold off
% xlabel('Radius (m)'); ylabel('Mass (kg)')
% % legend('Data','Mean','Median'); 
% % title('Payloads');
% % title('Rocket Bodies');
% % title('Debris');


end

function [g1,g2,g3] = getZeroGroups(inmatsat)
idx_mass = 8; idx_radius = 9; idx_objectclass = 23;

% mass vs radius
g1.allclass = []; % group 1: payloads; logical index
g2.allclass = []; % group 2: RBs
g3.allclass = []; % group 3: all debris

for ii = 1:12
    msinds = inmatsat(:,idx_objectclass) == ii; % all obj w/ objclass
    if ii == 1
        g1.allclass = find(msinds);                                 % all payload entries
        g1.zr = g1.allclass(inmatsat(g1.allclass,idx_radius) == 0); % index of zero radius
        g1.zm = g1.allclass(inmatsat(g1.allclass,idx_mass) == 0);   % index of zero mass
        g1.nz = setdiff(g1.allclass, union(g1.zr,g1.zm));           % index of non-zero radius and mass
        g1.nzno = g1.nz(~isoutlier(inmatsat(g1.nz,idx_radius)) ...  % non-zero, non-outlier
            & ~isoutlier(inmatsat(g1.nz,idx_mass)) );
    elseif ii == 5
        g2.allclass = find(msinds);
        g2.zr = g2.allclass(inmatsat(g2.allclass,idx_radius) == 0); % all RB entries
        g2.zm = g2.allclass(inmatsat(g2.allclass,idx_mass) == 0);
        g2.nz = setdiff(g2.allclass, union(g2.zr,g2.zm));
        g2.nzno = g2.nz(~isoutlier(inmatsat(g2.nz,idx_radius)) ...
            & ~isoutlier(inmatsat(g2.nz,idx_mass)) );
    else
        g3.allclass = union(g3.allclass, find(msinds));             % all debris entries
    end
end
g3.zr = g3.allclass(inmatsat(g3.allclass,idx_radius) == 0);
g3.zm = g3.allclass(inmatsat(g3.allclass,idx_mass) == 0);
g3.nz = setdiff(g3.allclass, union(g3.zr, g3.zm));
g3.nzno = g3.nz(~isoutlier(inmatsat(g3.nz,idx_radius)) ...
    & ~isoutlier(inmatsat(g3.nz,idx_mass)) );
end

function mat_sats = multiplyInitialPop(cfgMC,mat_sats,g1,g2,g3)
    inmatsat = mat_sats;
    mult = cfgMC.initpopMultiplier;
    % randomly choose which satellites to add to initial population (per object type)
    % note that physical parameters will be filled in by resampling (fillMassRadiusResample)
    if mult < 1  % simply select a subset of original matsats
        randind1 = randi(numel(g1.allclass), [floor(numel(g1.allclass) * mult),1]); % repeats allowed
        randind2 = randi(numel(g2.allclass), [floor(numel(g2.allclass) * mult),1]);
        randind3 = randi(numel(g3.allclass), [floor(numel(g3.allclass) * mult),1]);
        mat_sats = inmatsat([g1.allclass(randind1);g2.allclass(randind2);g3.allclass(randind3)],:); 
    elseif mult > 1 % if mult>1, ensure all of original matsats are included
        randind1 = randi(numel(g1.allclass), [floor(numel(g1.allclass) * (mult-1)),1]);
        randind2 = randi(numel(g2.allclass), [floor(numel(g2.allclass) * (mult-1)),1]);
        randind3 = randi(numel(g3.allclass), [floor(numel(g3.allclass) * (mult-1)),1]);
        extrasats = inmatsat([g1.allclass(randind1);g2.allclass(randind2);g3.allclass(randind3)],:); 
        % zero out mass and radius and fill in (resample)
        idx_mass = 8; idx_radius = 9;
        extrasats(:,[idx_mass, idx_radius]) = 0; 
        if cfgMC.fillMassRadius > 0
            addmatsats = fillMassRadiusResample(extrasats,g1,g2,g3);
        else
            addmatsats = extrasats;
        end
        % scramble argperigee and mean motion
        idx_argpo = 5; idx_mo = 6;  % range: [0,2pi]
        addmatsats(:,[idx_argpo, idx_mo]) = 2 * pi * rand(size(addmatsats,1),2);
        mat_sats = [mat_sats; addmatsats]; 
    end
end

function cfg = setupConstellation(cfg)
    if ~isempty(cfg.constellationFile)
        s = detectImportOptions(cfg.constellationFile,'Sheet',cfg.constellationSheet);
        tab = readtable(cfg.constellationFile,s);
        
        % fill empty entries
        for ind = 1:size(tab,1)
            if isnan(tab(ind,:).FirstLaunch)
                tab(ind,:).FirstLaunch = 2030;  % if no launch specified, start launch in 2030
                tab(ind,:).FinishLaunch = 2045;  % if no launch specified, end launch in 2045
            end
            if isnan(tab(ind,:).FinishLaunch)        
                tab(ind,:).FinishLaunch = tab(ind,:).FirstLaunch  + 10;  % if no end is specified, end launch in 15 years
            end
            if isnan(tab(ind,:).mass)
                tab(ind,:).mass = tab(1,:).mass;  % if no mass, copy Starlink-v1's        
            end
            if isnan(tab(ind,:).radius)
                tab(ind,:).radius = (tab(ind,:).mass / tab(1,:).mass)^(1/3) * tab(1,:).radius;  % use Starlink-v1's density
            end
        end
        % add missionlife from cfg
        tab.missionlife = ones(size(tab,1),1) * cfg.missionlifetime;
    else
        tab = [];
    end
    cfg.constellation = tab;

end
