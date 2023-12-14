%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MIT ARCLab's implementation of a space enivormental Simulation
%     MIT Orbital Capacity Tool - Monte Carlo (MOCAT-MC)
% 
% Authors: Richard Linares, Daniel Jang, Davide Gusmini, Andrea D'Ambrosio, 
%       Pablo Machuca, Peng Mun Siew
% https://github.mit.edu/arclab/orbitalrisk_MC
% 
% 08/23/2022: Initial commit
% 10/27/2022: Functionalized MOCAT-MC engine 
% 11/xx/2022: Vectorized subroutines
% 12/21/2022: main_mc.m created from main_mc_func_vec.m and main_mc_summary.m
% 
% Running Instruction:  Call the function with desired config file:
%   - main_mc('setupSomma')
%   - main_mc(setupSomma)
%   - main_mc('setupSomma', RNGseed) - run with a particular seed
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nS,nD,nN,nB,deorbitlist_r]=main_mc(MCconfig,RNGseed)
    % MCconfig is a structure of all relevant variables from a config setup file
    % usually in 'configFiles' folder
    % OR a string of the config name, e.g. 'setupSomma'
  
    % Initialize RNG seed
    count = 1; 

    if (exist('RNGseed', 'var'))
        rng(RNGseed);
        fprintf('main_mc specified with seed %i \n', RNGseed);
    elseif isfield(MCconfig, 'seed')
        rng(MCconfig.seed);
        fprintf('main_mc specified with config seed %i \n', MCconfig.seed);
    end

    % LOAD CFG 
    if ischar (MCconfig)
        MCconfig = eval(MCconfig);
    end
    loadCFG(MCconfig);  % see bottom    
    
    if isfield(MCconfig,'paramSSEM')
        param.paramSSEM = MCconfig.paramSSEM;        
    end
    
    if isfield(MCconfig,'sample_params')
        param.sample_params = MCconfig.sample_params;
    else
        param.sample_params = 0;        
    end
    
    % remove large data embedded in cfg (for saving)
    MCconfig.a_all = {};
    MCconfig.ap_all = {};
    MCconfig.aa_all = {};
    MCconfig.launchMC_step = {};

    % Assign constants to param structure for some functions
    param.req = radiusearthkm;
    param.mu = mu_const;
    param.j2 = j2;
    param.max_frag = max_frag;
    paramSSEM.species = [1,1,1,0,0,0];  % species: [S,D,N,Su,B,U];

    % Density profile
    param.density_profile = density_profile;
    
    % PREALLOCATE
    n_sats = size(mat_sats,1);
    
    % For figure display
    ind_fig = 0;
    p = [];
    p_earth = [];
    hAnnotation = [];
        
    numObjects = zeros(n_time,1); % count the number of total object vs time
    numObjects(1,1) = n_sats;
    numObjTrigger = 1e6;  % save when this number of obj is reached & bump up as it's reached
    count_coll = uint8(zeros(n_time,1)); % count the number of collisions vs time
    count_expl = uint8(zeros(n_time,1)); % count the number of explosions vs time
    count_debris_coll = uint8(zeros(n_time,1)); % number of collision debris vs time
    count_debris_expl = uint8(zeros(n_time,1)); % number of explosion debris vs time
    if save_output_file==3||save_output_file==4
        sats_info = cell(1,3);
    else
        sats_info = cell(n_time,3);         % contains info for SSEM binning
    end
    frag_info = cell(n_time,4);         % contains info on cube statistics
    frag_info5 = cell(n_time,1);        % contains info on collision statistics
    num_pmd = 0;
    num_deorbited = 0;
    launch = 0;
    out_future = {};
    count_tot_launches = 0;
    file_save_index = 0;
    S_MC = nan(n_time,numel(paramSSEM.R02)-1); D_MC = S_MC;  N_MC = S_MC;
    param_mean = nan((numel(paramSSEM.R02)-1)*4,3,n_time); param_var = param_mean;
    param_median = param_mean;
    
    % MATSATS DEFINITION
    getidx;

    param.maxID = max([max(mat_sats(:,idx_ID)),0]);
    
    %launches_current_step_vec and launches_current_step_eq_vec
    idx_launch_in_extra = [idx_ID]; %additional columns to create for new launches
    
    %prop_mit_vec
    %mat_sat_in = [a,ecco,inclo,nodeo,argpo,mo,Bstar]
    idx_prop_in = [idx_a,idx_ecco,idx_inclo,idx_nodeo,idx_argpo,idx_mo,idx_bstar,idx_controlled];
    idx_thalassa_in = [idx_a,idx_ecco,idx_inclo,idx_nodeo,idx_argpo,idx_mo,idx_bstar,idx_mass,idx_radius,idx_r,idx_v];
    %mat_sat_out = [a,ecco,inclo,nodeo,argpo,mo,errors,r_eci,v_eci];
    idx_prop_out = [idx_a,idx_ecco,idx_inclo,idx_nodeo,idx_argpo,idx_mo,idx_error,idx_r,idx_v];
    
    %orbcontrol_vec
    %mat_sat_in = [a,ecco,inclo,nodeo,argpo,mo,controlled,a_desired,missionlife,launched,r,v]
    idx_control_in = [idx_a,idx_ecco,idx_inclo,idx_nodeo,idx_argpo,idx_mo,idx_controlled,idx_a_desired,idx_missionlife,idx_launch_date,idx_r,idx_v];
    %mat_sat_out = [a,controlled,r,v];
    idx_control_out = [idx_a,idx_controlled,idx_r,idx_v];
    
    %frag_exp_SBM_vec
    %p1_in = [p1.mass,p1.radius,p1.r,p1.v,p1.objectclass]
    idx_exp_in = [idx_mass,idx_radius,idx_r,idx_v,idx_objectclass];
    
    %frag_col_SBM_vec
    %p1_in = [p1.mass,p1.radius,p1.r,p1.v,p1.objectclass]
    %p2_in = [p2.mass,p2.radius,p2.r,p2.v,p2.objectclass]
    idx_col_in = [idx_mass,idx_radius,idx_r,idx_v,idx_objectclass];

%     %Store sats_info into a file
%     sats_filename = 'sats_info_file.mat';
%     sats_file = matfile(sats_filename,'Writable',true);
%     sats_file.sats_info = sats_info; % contain info for SSEM binning
    
    objclassint_store = mat_sats(:,idx_objectclass);
    a_store = mat_sats(:,idx_a); 
    controlled_store = mat_sats(:,idx_controlled); 

    if save_output_file==3||save_output_file==4
        sats_info{1} = int8(objclassint_store);
        sats_info{2} = single(a_store);
        sats_info{3} = int8(controlled_store);    
        [popSSEM] = Fast_MC2SSEM_population(sats_info,paramSSEM);
        S_MC(1,:) = popSSEM(:,1);
        D_MC(1,:) = popSSEM(:,2);
        N_MC(1,:) = popSSEM(:,3);
    elseif save_output_file==6 
        [popSSEM,popSSEM_param_mean, popSSEM_param_var, popSSEM_param_median] = MC2SSEM_population_dist(mat_sats,paramSSEM);
        S_MC(1,:) = popSSEM(:,1);
        D_MC(1,:) = popSSEM(:,2);
        N_MC(1,:) = popSSEM(:,3);
        param_mean(:,:,1) = popSSEM_param_mean;
        param_var(:,:,1) = popSSEM_param_var;
        param_median(:,:,1) = popSSEM_param_median;
    else
        sats_info{1,1} = int8(objclassint_store);
        sats_info{1,2} = single(a_store);
        sats_info{1,3} = int8(controlled_store); 
    end
    
    % Extract SPECIE numbers of satellite, derelict, debris, rocket body
    [nS, nD, nN, nB] = categorizeObj(objclassint_store, controlled_store);

    if save_output_file > 0    % test save (check for directory and write permission)
        save([filename_save(1:end-4),'_part_1.mat'],'-v7.3','param');
        fprintf('Test saving of %s successful\n', filename_save);
    end

    fprintf('Year %i - Day %03i,\t PMD %04i,\t Deorbit %03i,\t Launches %03i,\t nFrag %03i,\t nCol %03i,\t nObjects %i (%i,%i,%i,%i)\n', ...
        year(time0),day(time0,'dayofyear'), num_pmd ,num_deorbited, ...
        length(out_future) * launch, count_expl(1), ...
        count_coll(1), numObjects(1,1), nS, nD, nN, nB);
    
    launch_data = [];

    %---------------------------------------------------------------------
    %START PROPAGATION
    for n = 2:n_time
        current_time = time0+days(tsince(n)/DAY2MIN);
        jd = juliandate(current_time);

        % LAUNCHES
        if strcmpi(launch_model,'matsat')  % repeat launches
            % repeatLaunches is a mat_sats of objects to launch again
            % launchRepeatYrs [2x1] sets the years that was used to filter above
            % plan: find n days span of repeatyrs; find dt; find rem/mod; get matsats; scramble dates
            repLauncht0 = min(repeatLaunches(:,16));
            if ~launchRepeatSmooth % not repeating launch
                repDaysRng = (launchRepeatYrs(2) - launchRepeatYrs(1) +1) * YEAR2DAY;
            else  % yearly uniform repeatign launch (smooth launch)
                repDaysRng = 1 * YEAR2DAY;
            end
            curLaunchT = mod(days(current_time - time0),repDaysRng);  % can be fractional, including days()
%            [ repeatLaunches(:,idx_launch_date) - repLauncht0  ] and [curLaunchT]  WILL be between 0 and YEAR2DAY
            repLaunchT = repeatLaunches(:,idx_launch_date) - repLauncht0;
            if curLaunchT + dt_days > YEAR2DAY   % figure out wrap-around issue
                launchInd = ( (repLaunchT  >= curLaunchT) & (repLaunchT  < YEAR2DAY) ) ...
                    | (repLaunchT < mod(curLaunchT + dt_days, YEAR2DAY)) ;
            else
                launchInd = (repLaunchT  >= curLaunchT ) & (repLaunchT < curLaunchT + dt_days);
            end
            launch = sum(launchInd) > 0;                % launch flag [0,1]
            out_future = repeatLaunches(launchInd,:);   % matsat
            currentJD = julian(year(current_time),month(current_time),day(current_time),0,0,0);
            % Scramble launch time within this time step (dt_days)
            out_future(:,idx_launch_date) = currentJD + rand(sum(launchInd),1) * dt_days; 
        elseif strcmpi(launch_model,'random') 
            [out_future,launch] = launches_current_step_vec(launch_model,time0,current_time,tsince,n,param,idx_launch_in_extra,...
                total_launch_per_year * (1+launch_increase_per_year)^floor(years(current_time - time0)), dt_days);
        elseif strcmpi(launch_model,'data') || contains(launch_model,'Somma')
            [out_future,launch] = launches_current_step_eq_vec(launch_model,time0,current_time,tsince,n,param,idx_launch_in_extra,...
                launchMC_step,ind_launch,additional_launches,ind_launch_add);
        elseif strcmpi(launch_model,'random_calendar') 
            if n<n_time
                [out_future,launch,launchMC_step] = launches_current_step_eq_vec(launch_model,time0,current_time,tsince,n,DAY2MIN,YEAR2DAY,param,idx_launch_in_extra,...
                    launchMC_step);
            else
                out_future = [];
            end
        elseif strcmp(launch_model,'no_launch')
            out_future = [];
        else
            error('no a valid launch model');
        end
        param.maxID = param.maxID+size(out_future,1);       %update maximum ID in population
        count_tot_launches = count_tot_launches+size(out_future,1);
    
    %PROPAGATION (one timestep at a time)
        n_sats = size(mat_sats,1);
        X_eci = zeros(n_sats,6);
        if ~exist('propagator','var')
            propagator = 'MIT';
        end
        
        if strcmpi(propagator,'SGP4')         % select between SGP4 or MIT's own propagator (deprecated at the moment)
            deorbit=[];
            for i = 1:n_sats
                [sats{i}, X_eci_temp, ~,~]=spg4_ecf(sats{i},tsince(n)); %tsince(n) tsince(n)

                % deorbit if spg4 outputs zero position vector or if satellite is
                % below 90 km altitude
                if norm(X_eci_temp(1,:))==0 || ((sats{i}.a*radiusearthkm)*(1-sats{i}.ecco)-radiusearthkm)<150 || norm(X_eci_temp(1:3,:))<(radiusearthkm+100) || sats{i}.error~=0 || sats{i}.a<0
                    deorbit=[deorbit;i];
                    X_eci(i,:)=0;    
                    sats{i}.r=[0 0 0];
                    sats{i}.v=[0 0 0];
                else
                    if norm(X_eci_temp(4:6,1))>30 % <<<<<<<<<<????
                        warning('satellite velocity > 30 km/s! ')
                    end
                    X_eci(i,:)=X_eci_temp;
                    sats{i}.r=X_eci_temp(1:3,1)';
                    sats{i}.v=X_eci_temp(4:6,1)';
                end
            end
        elseif strcmpi(propagator,'THALASSA') % use thalassa
            param.jd=jd;
            if n>1
                dt = 60*(tsince(n)-tsince(n-1));% units of time in seconds
            else
                dt = 60*tsince(n);% units of time in seconds
            end
            for i = 1:n_sats
                mat_sats(i,idx_prop_out) = prop_thalassa(mat_sats(i,idx_thalassa_in),dt,param);
            end
            deorbit = find(mat_sats(:,idx_r(1))==0 | (mat_sats(:,idx_a)*radiusearthkm).*(1-mat_sats(:,idx_ecco))<(150+radiusearthkm) | sqrt(sum(mat_sats(:,idx_r).^2,2))<(radiusearthkm+100) | mat_sats(:,idx_error)~=0 | mat_sats(:,idx_a)<0);           
        else % use prop_mit 
            param.jd=jd;
            
            %prop_mit_vec
            %mat_sat_in = [a,ecco,inclo,nodeo,argpo,mo,Bstar]
            %mat_sat_out = [a,ecco,inclo,nodeo,argpo,mo,errors,r_eci,v_eci];
            if n>1
                [mat_sats(:,idx_prop_out)] = prop_mit_vec(mat_sats(:,idx_prop_in),60*(tsince(n)-tsince(n-1)),param);% units of time in seconds
            else
                [mat_sats(:,idx_prop_out)] = prop_mit_vec(mat_sats(:,idx_prop_in),60*tsince(n),param);% units of time in seconds
            end

            % deorbit if spg4 outputs zero position vector or if satellite is
            % below 90 km altitude
            deorbit = find(mat_sats(:,idx_r(1))==0 | (mat_sats(:,idx_a)*radiusearthkm).*(1-mat_sats(:,idx_ecco))<(150+radiusearthkm) | sqrt(sum(mat_sats(:,idx_r).^2,2))<(radiusearthkm+100) | mat_sats(:,idx_error)~=0 | mat_sats(:,idx_a)<0);           
        end     
        % Remove all deorbited objects
        % deorbit_list = [deorbit_list;deorbit];
        num_deorbited = length(deorbit);
        mat_sats(deorbit,:) = [];
        
    %ORBIT CONTROL        
        if mod(n,step_control)==1 || step_control==1
            
            %orbcontrol_vec
            %mat_sat_in = [a,ecco,inclo,nodeo,argpo,mo,controlled,a_desired,missionlife,launched,r,v]
            %mat_sat_out = [a,controlled,r,v];
            [mat_sats(:,idx_control_out),deorbit_PMD] = orbcontrol_vec(mat_sats(:,idx_control_in),tsince(n),time0,orbtol,PMD,DAY2MIN,YEAR2DAY,param);
            
            % Remove all post-mission disposal satellites
            num_pmd = length(deorbit_PMD);
            mat_sats(deorbit_PMD,:) = [];
        else
            num_pmd = 0;
        end
    
    %EXPLOSIONS (for Rocket Body)
        n_sats = size(mat_sats,1);
        out_frag = [];
        find_rocket = find(mat_sats(:,idx_objectclass)==5); %Rocket bodies
        rand_P_exp = rand(numel(find_rocket),1); %random numbers for each rocket body
        if exist('P_frag_cutoff','var')
            find_P_exp = find(rand_P_exp<P_frag & ...
                (year(current_time) - jd2date(mat_sats(find_rocket,idx_launch_date))) < P_frag_cutoff);
        else
            find_P_exp = find(rand_P_exp<P_frag);
        end
        remove_frag = find_rocket(find_P_exp); %Pre-identify objects to remove
        for idx_P_exp_temp = numel(find_P_exp):-1:1 %reverse order in case elements need to be deleted from remove_frag variable
            idx_P_exp = remove_frag(idx_P_exp_temp);
            
            p1_all = mat_sats(idx_P_exp,:);
            p1_mass = p1_all(1,idx_mass);
            p1_objectclass = p1_all(1,idx_objectclass);

            p1_in = p1_all(1,idx_exp_in);

            %frag_exp_SBM_vec
            %p1_in = [p1.mass,p1.radius,p1.r,p1.v,p1.objectclass]
            [debris1] = frag_exp_SBM_vec(tsince(n), p1_in, param);
            param.maxID = param.maxID+size(debris1,1); %update maximum ID in population, since new objects have been added
                        
            if isempty(debris1) %if no debris is generated, do not remove object from simulation
                remove_frag(idx_P_exp_temp) = [];
            else %if some debris is generated, add it to out_frag matrix                       
                out_frag = [out_frag; debris1];
                
                fprintf(2,'Year %i - Day %03i \t Explosion, p1 type %s, %0.1f kg, nDebris %i\n', ...
                    year(current_time), day(current_time,'dayofyear'), ...
                    p1_objectclass, p1_mass, size(debris1,1));
                count_expl(n) = count_expl(n)+1;      
            end
        end
        mat_sats(remove_frag,:) = [];

    %COLLISIONS
        if ( exist('skipCollisions','var') && skipCollisions == 1 ) || isempty(mat_sats)
%             n_res = 0;
            collision_array = [];
        else
%             res = cube_vec(mat_sats(:,idx_r), CUBE_RES);          % do CUBE
%             n_res = length(res);
%             [duplicates,unique_group,duplicate_idx]=cube_vec_v2(mat_sats(:,idx_r), CUBE_RES);
%             n_res = length(unique_group);
            collision_cell = cube_vec_v3(mat_sats(:,idx_r), CUBE_RES, collision_alt_limit);          % do CUBE
            collision_array = cell2mat(collision_cell);
%             n_res = length(res);
        end 

        remove_collision = [];
        out_collision = [];

        if ~isempty(collision_array)
            p1_idx = collision_array(:,1);
            p2_idx = collision_array(:,2);
            p1_all = mat_sats(p1_idx,:);
            p2_all = mat_sats(p2_idx,:);

%             % Sanity check
%             del_r= p1_all(:,idx_r)- p2_all(:,idx_r);
%             if max(abs(del_r(:))) > CUBE_RES
%                 max(abs(del_r))
%             end

            p1_controlled = p1_all(:,idx_controlled); p1_radius = p1_all(:,idx_radius); p1_v = p1_all(:,idx_v); 
            p2_controlled = p2_all(:,idx_controlled); p2_radius = p2_all(:,idx_radius); p2_v = p2_all(:,idx_v); 
            % probability of collision
            Pij = collision_prob_vec(p1_radius, p1_v, p2_radius, p2_v, CUBE_RES);
            % probability of collision over dt: 5 days
            % P = 4/3*pi*((p1.radius^3*p2.radius^3)/CUBE_RES^3)*Pij*dt_days*DAY2SEC;
            P = zeros(size(p1_controlled));
            sum_controlled = p1_controlled + p2_controlled;
            check_0 = sum_controlled<0.5;
            find_not0 = find(~check_0);
            check_1 = sum_controlled(find_not0)<1.5;
            find_1 = find_not0(check_1);
            find_2 = find_not0(~check_1);
    %         if (p1_controlled + p2_controlled) == 1
            P(find_1) = Pij(find_1)*(alph*dt_days*DAY2SEC);
    %         elseif (p1_controlled + p2_controlled) == 2
            P(find_2) = Pij(find_2)*(alph_a*dt_days*DAY2SEC);
    %         else
            P(check_0) = Pij(check_0)*(dt_days* DAY2SEC);
    %         end
            rand_P = rand(size(p1_controlled));

            % Monte-carlo simulation - there is a collision if random number is lower than prob of collision
            find_P = find(rand_P<P);
            frag4tmp = zeros(numel(find_P),3);    % col4:   index for collision, debAn, debBn               [Mx3] uint16
            frag5tmp = zeros(numel(find_P),5);    % "col5": p1 mass, p2 mass (0.1 kg), alt10km, debAn, debBn [Mx5] uint16
            for idx_P_temp = 1:numel(find_P)
                idx_P = find_P(idx_P_temp);
                p1_mass = p1_all(idx_P,idx_mass);
                p2_mass = p2_all(idx_P,idx_mass);
                p1_objectclass = p1_all(idx_P,idx_objectclass);
                p2_objectclass = p2_all(idx_P,idx_objectclass);
                p1_r = p1_all(idx_P,idx_r);  p2_r = p2_all(idx_P,idx_r);
    
                p1_in = p1_all(idx_P,idx_col_in); 
                p2_in = p2_all(idx_P,idx_col_in); 
    
                %frag_col_SBM_vec
                %p1_in = [p1.mass,p1.radius,p1.r,p1.v,p1.objectclass]
                %p2_in = [p2.mass,p2.radius,p2.r,p2.v,p2.objectclass]
                [debris1, debris2] = frag_col_SBM_vec(tsince(n), p1_in, p2_in, param);
                param.maxID = param.maxID+size(debris1,1)+size(debris2,1); %update maximum ID in population, since new objects have been added
    
                out_collision = [out_collision; debris1; debris2];          % add debris
    
                if ~isempty(debris1) || ~isempty(debris2)
                    fprintf(2,'Year %i - Day %03i \t Collision, p1 type %i, %0.1f kg, p2 type %i, %0.1f kg, nDebrisA %i, nDebrisB %i, altkm %0.0f\n', ...
                        year(current_time), day(current_time,'dayofyear'), ...
                        p1_objectclass, p1_mass, p2_objectclass,  p2_mass, ...
                        size(debris1,1), size(debris2,1), mean([norm(p1_r), norm(p2_r)])-radiusearthkm );
                    count_coll(n) = count_coll(n)+1;
                    remove_collision = [remove_collision; p1_idx(idx_P); p2_idx(idx_P)]; % remove parent objects    
                end
                % col4:   index for collision, debAn, debBn               [Mx3] uint16
                frag4tmp(idx_P_temp,1) = idx_P;
                frag4tmp(idx_P_temp,2) = size(debris1,1);
                frag4tmp(idx_P_temp,3) = size(debris2,1);
                % "col5": p1 mass, p2 mass (0.1 kg), alt10km, debAn, debBn [Mx5] uint16
                frag5tmp(idx_P_temp,1) = p1_mass * 10;
                frag5tmp(idx_P_temp,2) = p2_mass * 10;
                frag5tmp(idx_P_temp,3) = (mean([norm(p1_r), norm(p2_r)])-radiusearthkm)/10;
                frag5tmp(idx_P_temp,4) = size(debris1,1);
                frag5tmp(idx_P_temp,5) = size(debris2,1);
            end
        end

    %DATA PROCESSING
        % note this is put at the end because we assume these object don't
        % generate collisions during their first timestep of existance
        mat_sats(remove_collision,:) = [];
        mat_sats = [mat_sats; out_future; out_frag; out_collision];
        % record launch data
        launch_data = [launch_data; out_future];

    %ACCOUNTING
        n_sats = size(mat_sats,1);
        numObjects(n,1) = n_sats;
        
        objclassint_store = mat_sats(:,idx_objectclass);
        a_store = mat_sats(:,idx_a); 
        controlled_store = mat_sats(:,idx_controlled); 
        
        if save_output_file==3||save_output_file==4
            sats_info{1} = int8(objclassint_store);
            sats_info{2} = single(a_store);
            sats_info{3} = int8(controlled_store);    
            [popSSEM] = Fast_MC2SSEM_population(sats_info,paramSSEM);
            S_MC(n,:) = popSSEM(:,1);
            D_MC(n,:) = popSSEM(:,2);
            N_MC(n,:) = popSSEM(:,3);
        elseif save_output_file==6 
            [popSSEM,popSSEM_param_mean, popSSEM_param_var, popSSEM_param_median] = MC2SSEM_population_dist(mat_sats,paramSSEM);
            S_MC(n,:) = popSSEM(:,1);
            D_MC(n,:) = popSSEM(:,2);
            N_MC(n,:) = popSSEM(:,3);
            param_mean(:,:,n) = popSSEM_param_mean;
            param_var(:,:,n) = popSSEM_param_var;
            param_median(:,:,n) = popSSEM_param_median;
        else
            sats_info{n,1} = int8(objclassint_store);
            sats_info{n,2} = single(a_store);
            sats_info{n,3} = int8(controlled_store);
        end

        

%         sats_file.sats_info(n,1) = int8(objclassint_store);
%         sats_file.sats_info(n,2) = single(a_store);
%         sats_file.sats_info(n,3) = int8(controlled_store);

        if ~isempty(collision_array) && save_output_file == 3
            frag_info5{n} = uint16(frag5tmp);
            
        elseif ~isempty(collision_array) && (save_output_file == 4 || save_output_file == 5 || save_output_file == 7)
            % save fragmentation info:  
            % col1: P, p1 satno, p2 satno                           [Nx3] single
            % col2: p1 mass, p2 mass (0.1 kg)                       [Nx2] uint16
            % col3: objclass, objclass, dv (0.1 km/s), alt10km      [Nx4] uint8
            % col4: index for collision, debAn, debBn               [Mx3] uint16
            % "col5": p1 mass, p2 mass (0.1 kg), alt10km, debAn, debBn [Mx5] uint16
            frag_info{n,1} = single([P, p1_all(:, idx_ID), p2_all(:, idx_ID)]);
            frag_info{n,2} = uint16([p1_all(:, idx_mass), p2_all(:, idx_mass)] * 10);
            frag_info{n,3} = uint8([p1_all(:, idx_objectclass), p2_all(:, idx_objectclass), ...
                                vecnorm((p1_all(:,idx_v) - p2_all(:,idx_v))')' * 10, ...
                                (mean([vecnorm(p1_all(:,idx_r),2,2), vecnorm(p2_all(:,idx_r),2,2)],2)-radiusearthkm)/10 ...
                                ]);
            frag_info{n,4} = single(frag4tmp);
            frag_info5{n} = uint16(frag5tmp);
        end
            
    
        count_debris_coll(n,1) = size(out_collision,1);
        count_debris_expl(n,1) = size(out_frag,1);



        [nS, nD, nN, nB] = categorizeObj(objclassint_store, int8(controlled_store));
        fprintf('Year %i - Day %03i,\t PMD %04i,\t Deorbit %03i,\t Launches %03i,\t nFrag %03i,\t nCol %03i,\t nObjects %i (%i,%i,%i,%i)\n', ...
            year(current_time), day(current_time,'dayofyear'), num_pmd ,num_deorbited, ...
            size(out_future,1) * launch, count_expl(n), ...
            count_coll(n), numObjects(n,1), nS, nD, nN, nB);
        
        deorbit_list(count) = num_deorbited;
        if count > 1
            deorbit_list(count) = deorbit_list(count) + deorbit_list(count-1);
        end
        count = count + 1;

    %ANIMATION
%         if strcmpi(animation,'yes')
%             [ind_fig,p,p_earth,hAnnotation] = plot_fig(mat_sats(:,idx_r),tsince,n,omega_earth,ind_fig,p,p_earth,hAnnotation);
%         end   
    
    % CHECK MEMORY CONSUMPTION
        sat_bytes = whos('sats_info').bytes;
        frag_info_bytes = whos('frag_info').bytes;
        frag_info5_bytes = whos('frag_info5').bytes;
        total_bytes = sat_bytes+frag_info_bytes+frag_info5_bytes;

    % SAVE OUTPUT 
        if save_output_file > 0 && ( (mod(n,n_save_checkpoint) == 1 || n==n_time || numObjects(n,1) > numObjTrigger || total_bytes/(1e+9) > 30) )
            if total_bytes/(1e+9) > 30
                warning('Total file size exceeded 30GB. Total file size is %s GB.',total_bytes/(1e+9));
            end

            if numObjects(n,1) > numObjTrigger
                numObjTrigger = numObjTrigger + 5e5;  % increase limit
            end
            
            % Increment file save index by 1
            file_save_index = file_save_index +1;

            switch save_output_file
                case 1  % save entire workspace
                    savevars = {'*'}; 
                    save([filename_save(1:end-4),'_part_',num2str(file_save_index),'.mat'],'-v7.3',savevars{:});
                case 2  % save sats_info, config, params
                    savevars = {'sats_info','MCconfig','param','paramSSEM'};
                    save([filename_save(1:end-4),'_part_',num2str(file_save_index),'.mat'],'-v7.3',savevars{:});
                case 3  % save summary (S_MC, N_MC, D_MC, etc)
                    savevars = {'S_MC','D_MC','N_MC','MCconfig','param','paramSSEM','frag_info5'};
%                     % preallocate
%                     S_MC = nan(n_time,numel(paramSSEM.R02)-1); D_MC = S_MC;  N_MC = S_MC;
%                     for tind = 1:n
%                         [popSSEM] = Fast_MC2SSEM_population(sats_info(tind,:),paramSSEM);
%                         S_MC(tind,:) = popSSEM(:,1);
%                         D_MC(tind,:) = popSSEM(:,2);
%                         N_MC(tind,:) = popSSEM(:,3);
%                     end
                    MCconfig.mat_sats = [];  % to save space
                    save([filename_save(1:end-4),'_part_',num2str(file_save_index),'.mat'],'-v7.3',savevars{:});
                case 4  % save summary and collision stats
                    savevars = {'S_MC','D_MC','N_MC','MCconfig','param','paramSSEM','frag_info','frag_info5'};
%                     S_MC = nan(n_time,numel(paramSSEM.R02)-1); D_MC = S_MC;  N_MC = S_MC;
%                     for tind = 1:n
%                         [popSSEM] = Fast_MC2SSEM_population(sats_info(tind,:),paramSSEM);
%                         S_MC(tind,:) = popSSEM(:,1);
%                         D_MC(tind,:) = popSSEM(:,2);
%                         N_MC(tind,:) = popSSEM(:,3);
%                     end
                    MCconfig.mat_sats = [];  % to save space
                    save([filename_save(1:end-4),'_part_',num2str(file_save_index),'.mat'],'-v7.3',savevars{:});
                case 5  % save just collision stats
                    savevars = {'MCconfig','param','paramSSEM','frag_info','frag_info5'};
                    MCconfig.mat_sats = [];  % to save space
                    save([filename_save(1:end-4),'_part_',num2str(file_save_index),'.mat'],'-v7.3',savevars{:});
                case 6  % Tracking mean/var/median of physical parameters
                    savevars = {'S_MC','D_MC','N_MC','MCconfig','param','paramSSEM','param_mean','param_var','param_median'};
                    MCconfig.mat_sats = [];  % to save space
                    save([filename_save(1:end-4),'_part_',num2str(file_save_index),'.mat'],'-v7.3',savevars{:});
                case 7  % Collision tracking (UROP)        
                    frag_info7 = cell(n_time,2); % just 2 columns needed
                    frag_info7(:,1) = frag_info(:,1);  % 1st col
                    frag_info7(:,2) = frag_info(:,4);  % 4th col
                    savevars = {'frag_info7'};
                    MCconfig.mat_sats = [];  % to save space
                    save([filename_save(1:end-4),'_part_',num2str(file_save_index),'.mat'],'-v7.3',savevars{:});
                    clear frag_info7;
                case 10 % save all, decimate matsats
                    savevars = {'MCconfig','param','paramSSEM','matsatsperN'};
                    MCconfig.mat_sats = [];  % to save space
                    save([filename_save(1:end-4),'_part_',num2str(file_save_index),'.mat'],'-v7.3',savevars{:});
                case 11 % save all, decimate matsats
                    savevars = {'MCconfig','param','paramSSEM','mat_sats','nS','nD','nN','nB'};
                    MCconfig.mat_sats = [];  % to save space
                    save([filename_save(1:end-4),'_part_',num2str(file_save_index),'.mat'],'-v7.3',savevars{:});
                otherwise
                    fprintf('save_output_file flag set to unsupported value: %i\n', save_output_file);
            end
    
            if save_output_file==3||save_output_file==4
                sats_info = cell(1,3);
            else
                sats_info = cell(n_time,3);         % contains info for SSEM binning
            end
            frag_info = cell(n_time,4);         % contains info on cube statistics
            frag_info5 = cell(n_time,1);        % contains simplified info on cube statistics
            fprintf('Variables cleared and %s saved', [filename_save(1:end-4),'_part_',num2str(file_save_index),'.mat']);
        end
    end
    % figure
    % plot(tsince/DAY2MIN/YEAR2DAY,sum(D_MC,2))
    % title('D')
    % figure
    % plot(tsince/DAY2MIN/YEAR2DAY,sum(S_MC,2))
    % title('S')
    % figure
    % plot(tsince/DAY2MIN/YEAR2DAY,sum(N_MC,2))
    % title('N')
    % figure
    % plot(tsince/DAY2MIN/YEAR2DAY,sum(D_MC,2)+sum(S_MC,2)+sum(N_MC,2))
    % title('T')
    % 
    % figure
    % plot(tsince/DAY2MIN/YEAR2DAY,D_MC)
    % title('D per shell')

    deorbitlist_r = deorbit_list;
    
    fprintf('\n === FINISHED MC RUN (main_mc.m) WITH SEED: %i === \n',RNGseed);

end

function loadCFG(cfg)
    fns = fieldnames(cfg);
    for ind = 1:numel(fns)
        if strcmp(fns{ind}, 'param')
            param = cfg.param;
            assignin('caller','param', param);
        else
            assignin('caller',fns{ind}, cfg.(fns{ind}));
        end
    end
end