function [debris1]=frag_exp_SBM_vec(ep, p1_in,param)
% Modified from fragmentation.m by djang in Oct 2022
% More edits by pmachuca
%p1_in = [p1.mass,p1.radius,p1.r,p1.v,p1.objectclass]


% inputs
%   ep: epoch
%   p1: object 1 (fields used: mass, radius, v, objectclass)
%   param: holds param.max_frag, param.mu, param.req


% PARAMETERS TO EXPELL FROM FUNCTION:
% LB (frag_exp_SBM, frag_col_SBM)



% Fragmentation (Explosion) model following 
% NASA EVOLVE 4.0 standard breakup model (2001)
% with the revision in ODQN "Proper Implementation of the 1998 NASA Breakup
% Model" (2011)

% Explosions in the SBM follows the power law of the independent variable:
% the characteristic length (Lc).  Mass is controlled by adding in large
% pieces of debris that account for most of the mass of the original object

% Quote: 
% In the region above 1 m the data shows the deposit of several large pieces
% that do not necessarily follow the power law distribution. These would realistically be larger,
% more massive components farther from the explosion center (e.g., remnants of equipment
% shelves, pressurant tanks, nozzle bells, etc.). In fact these fragments account for the bulk of
% fragment mass. Based on this understanding, the correct implementation of the 1998 NASA
% Breakup Model includes the distribution of fragments from 1 mm to 1 m following the
% power law distribution in Equation 1, with an additional two to eight large fragments after
% 1 m, keeping mass conserved

% Note that the SBM distribution is modified by ESA for MASTER-2009 and MASTER-8
% and also by other authors in other papers e.g. Cimmino (2021), etc

% Plotting the power-law distribution of explosions by SBM; compare to Fig 1 of ODQN
% In the SBM, this region is the same for all explosions between 0.1 and 1 m
%     Lcs = logspace(-2,1,1000);
%     loglog(Lcs, 6 * Lcs.^-1.6); grid on
%     ylim([1,10000]);  set(gca, 'XTickLabel',get(gca,'XTick')) ; set(gca, 'YTickLabel',get(gca,'YTick')) 
%     legend('NASA: Ncum = 6 Lc^-1.6');

LB = 0.1;       % trackable diameter dtr = 10 cm; also, characteristic length L_c later

p1_mass = p1_in(1); p1_radius = p1_in(2);
p1_r = p1_in(3:5); p1_v = p1_in(6:8);
p1_objclass = p1_in(9);

M = p1_mass;

% FIX: Add all classes for objects
% reference: MASTER-8-Final-Report (pg 65/465, tab 2.1)
% Scaling factors for the revised NASA breakup model as used for future projections in EVOLVE 4.0
% Parent Object Type                                                    s
% PROTON (SL-12) ullage motors (SOZ units)                              0.1
% Other (non-SOZ) rocket bodies                                         1.0
% EORSATs (Soviet/Russian Electronic Ocean Reconnaissance Satellites)   0.6
% Molniya type early warning satellites                                 0.1
% All Soviet/Russian battery-related events                             0.5
% All Soviet/Russian anti-satellite tests (ASAT)                        0.3
% Other payloads                                                        1.0
% NOTE: this power law was fit to RB between 600 and 1000 kg

cs = 1;      % scaling factor 's' in MASTER and NASA SBM; set at 1 for "others" always

% reference: SPACE DEBRIS, Model and Risk Analysis
% Cumulative number of debris between LB and 1 m diameter
num = floor(6 * cs * LB ^(-1.6) - 6 * cs * min([1,2*p1_radius]) ^(-1.6));       % eq(3.32);  LB is L_c;  

% NEW METHOD #3:create PDF, then sample 'num' selections
%                   but only up to 1m, then randomly sample larger objects

dd_edges = logspace(log10(LB),log10(min([1,2*p1_radius])),200);               % sampled diameters ; power law valid from LB to 1 m
log10_dd = log10(dd_edges); %log10 of diameter edge bins
dd_means = 10.^(log10_dd(1:end-1)+diff(log10_dd)/2); %mean value of each diameter edge bin (in log scale, not linear scale, since bins are generated in a logarithmic scale)

nddcdf = 6 * cs * dd_edges.^-1.6;            % CUMULATIVE distribution of n's  (eq 2.68)
ndd = max(0,-diff(nddcdf)); % diff to get PDF count for the bins (dd_edges)
floor_ndd = floor(ndd); %floor of PDF count for each bin
rand_sampling = rand(size(ndd)); %0 to 1, random number for stochastic sampling of fragment diameters
add_sampling = rand_sampling>(1-(ndd-floor_ndd)); %0 if random number is lower than 1-(decimal part of ndd), 1 if random number is above 1-(decimal part)
d_pdf = repelem(dd_means,floor_ndd+add_sampling)';   % PDF of debris objects between LB and 1 m

d = d_pdf(randperm(numel(d_pdf)));% Do not limit number of fragments to be equal to 'num'
% d = d_pdf(randi(numel(d_pdf),ceil(num),1)); % Limit number of fragments to be equal to 'num'

% figure(50)
% histogram(d,dd_edges)
% hold on
% plot(dd_means,ndd,'k')
% set(gca,'YScale','log','XScale','log')
% xlabel('d [m]')
% ylabel('Number of fragments [-]')
% title('Number of fragments vs diameter')
% 
%     figure(22); clf; subplot(321); 
%     Lcs = logspace(-2,1,100); loglog(Lcs, 6 * Lcs.^-1.6); grid on; ylim([1,10000]); 
%     hold on; loglog(dd_edges(1:end),nddcdf); xlabel('diam'); ylabel('CDF'); title('theoretical reverse CDF');
%     legend('NASA: N_{cum} = 6s Lc^{-1.6}','nddcdf from code');
%     subplot(323); histogram(d_pdf); xlabel('diam'); ylabel('count'); title('PDF to sample from')
%     subplot(325); loglog(dd_means, flip(cumsum(flip(histcounts(repelem(dd_means,round(ndd))',dd_edges)))),'-x'); xlabel('Diam (m)'); ylabel('cumulative count'); title('CDF of above (for shape)');
%     subplot(324); histogram(d); xlabel('diam'); ylabel('count'); title('PDF of sampled d');
%     subplot(326); loglog(dd_edges, 6 * dd_edges.^-1.6); hold on; loglog(dd_means, flip(cumsum(flip(histcounts(d,dd_edges)))),'-x'); xlabel('Diam (m)'); ylabel('cumulative count'); title('CDF of sampled debris');
    
% %Old code:
% nddcdf = 6 * cs * dd_edges.^-1.6;            % CUMULATIVE distribution of n's  (eq 2.68)
% ndd = -diff(nddcdf);                    % diff to get PDF count
% d_pdf = repelem(dd(1:end-1),ceil(ndd*10000))';   % PDF of debris objects between LB and 1 m
% d = d_pdf(randperm(numel(d_pdf),ceil(num)));
% % note that this dist is the same for all objects no matter mass or size of parent object

% calculate mass of objects [LB, 1 m] by d > A > Am > m
A = 0.556945*d.^(2.0047077);            % calculate area; Eq 2.72
Am = func_Am(d, p1_objclass);
m = A./Am;
% sum(m);       % about 100 kg is accounted for debris < 1 m

if sum(m) < M        % if debris mass is less than total mass,
    % Randomly select "additional two to eight large fragments after 1 m,
    % keeping mass conserved" per ODQN.  "these fragments account for the bulk
    % of fragment mass"

    % Assign remnant mass > d > A > Am > dv
    m_remSum = M - sum(m);       % remnant mass
    remDist = rand(randi([2,8]),1);         % randomly distribute by random number [2,8] <<<<<<<
    m_rem = m_remSum * remDist / sum(remDist); % randomly distribute the remnant mass
    d_rem_approx = (m_rem ./ p1_mass * p1_radius^3).^(1/3) * 2; % use same density (mass / r^3) as original object to get d_rem <<<<<<<<<<
    Am_rem = func_Am(d_rem_approx, p1_objclass); %compute A/m using approximate diameters
    A_rem = m_rem.*Am_rem; %compute area from mass and A/m
    d_rem = d_rem_approx; %use approximate diameter based on parent object
    %remove remnant objects that are too small
    idx_too_small = find(d_rem<LB & m_rem<M/1000); %identify remnant objects assigned to larger satellite that are smaller than lower bound Lc
    m_rem(idx_too_small) = []; %remove remnant objects that are too small from mass vector
    d_rem(idx_too_small) = []; %remove remnant objects that are too small from diameter vector
    A_rem(idx_too_small) = []; %remove remnant objects that are too small from area vector
    Am_rem(idx_too_small) = []; %remove remnant objects that are too small from A/m vector
else
    % method 2) sort by Lc, keep smallest objects until mass adds up <<<<<<<<<<<<<<
    [~, dord] = sort(m);
    cumsum_m = cumsum(m(dord)); %cumulative sum of fragment masses in ascending order of diameter
    lastidx = find(cumsum_m < M, 1, 'last');  % last valid index of sorted list
    cumsum_lower = max([0,cumsum_m(lastidx)]); %cumulative mass of ordered fragments that remain below M
    valididx = dord(1:lastidx); %indices of fragments fulfilling mass constraint    
    m = m(valididx);   %Select elements fulfilling mass constraint
    d = d(valididx);
    A = A(valididx);
    Am = Am(valididx); 
    
    %Check if there is mass remaining, and generate an additional fragment if needed
    m_rem = M-cumsum_lower; %remaining mass to accumulate a total mass of M
    if m_rem>M/1000 %if the remaining mass is larger than smallest fragment, add one more fragment (unless it is smaller than Lc)
        d_rem_approx = (m_rem ./ p1_mass * p1_radius^3).^(1/3) * 2; % use same density (mass / r^3) as larger satellite to get approximate d_rem <<<<<<<<<<
        Am_rem = func_Am(d_rem_approx, p1_objclass); %compute A/m using approximate diameter
        A_rem = m_rem*Am_rem; %compute area from mass and A/m
        d_rem = d_rem_approx; %use approximate diameter based on parent object
        if d_rem<LB && m_rem<M/1000 %if remaining fragment is smaller than lower bound Lc, do not generate an additional fragment
            d_rem = []; A_rem = []; Am_rem = []; m_rem = []; % no "remnants" exist
        end
    else %if there is no need for an additional fragment, do not generate one       
        d_rem = []; A_rem = []; Am_rem = []; m_rem = []; % no "remnants" exist
    end
end

% figure
% scatter(d,m,'kx')
% hold on
% scatter(d_rem,m_rem,'rx')
% scatter(p1_radius*2,p1_mass,'ro')
% legend('General frags.','p1 frags.','p1','Location','northwest')
% % plot([d;d_rem],0.556945*[d;d_rem].^2.0047077./10.^(-0.925),'k-')
% set(gca,'XScale','log','YScale','log')
% title('m vs d')
% xlabel('d')
% ylabel('m')
% 
% figure
% scatter(d,Am,'kx')
% hold on
% scatter(d_rem,Am_rem,'rx')
% scatter(p1_radius*2,pi*p1_radius^2/p1_mass,'ro')
% scatter(p1_radius*2,(2*p1_radius)^2/p1_mass,'rs')
% legend('General frags.','p1 frags.','p1 sphere','p1 cube')
% set(gca,'XScale','log','YScale','log')
% title('A/m vs d')
% xlabel('d')
% ylabel('A/m')
% 
% figure
% scatter(d.^3,m,'kx')
% hold on
% scatter(d_rem.^3,m_rem,'rx')
% scatter((p1_radius*2)^3,p1_mass,'ro')
% legend('General frags.','p1 frags.','p1','Location','northwest')
% set(gca,'XScale','log','YScale','log')
% title('m vs d^3')
% xlabel('d^3')
% ylabel('m')

% Assign dv to random directions; create samples on unit sphere
dv = func_dv([Am; Am_rem],'exp')/1000;       % km/s
u=rand(length(dv),1)*2-1;
theta=rand(length(dv),1)*2*pi;

v = sqrt(1 - u.^2);
p = [v.*cos(theta) v.*sin(theta) u];

% dv velocity vectors
try
    dv_vec = [p(:,1).*dv  p(:,2).*dv  p(:,3).*dv];
catch
    warning('error assigning dv to unit sphere')
    return
end

% create fragments
fragments = [[d; d_rem] [A; A_rem] [Am; Am_rem] [m; m_rem] dv dv_vec(:,1) dv_vec(:,2) dv_vec(:,3)];
% DEBUG: HISTOGRAM OF THESE ^
%     figure(10);clf;
%     subplot(511); histogram(d,100); xlabel('d (m)');
%     subplot(512); histogram(A,100); xlabel('A (m^2)');
%     subplot(513); histogram(Am,100); xlabel('Am (m^2/kg)');
%     subplot(514); histogram(m,100); xlabel('m (kg)');
%     subplot(515); histogram(dv * 1000,[0:10:1000]); xlabel('dv (m/s)');
%     subplot(511); title('fragment distributions', sprintf('Original mass: %0.1f kg, radius: %0.1f m',p1_mass, p1_radius));
%     pause;

if abs(sum([m; m_rem]) - M) > M*0.05
    warning('Total sum of debris mass (%0.1f kg) differs from "mass" of original objects (%0.1f kg)', ...
        sum([m; m_rem]), M);
end

if any(fragments(:,1) > p1_radius*2*1.00001)
    warning('Some fragments (%0.1f m) exceed the diameter of the original object (%0.1f m)', ...
        max(fragments(:,1)),p1_radius*2);
end

%To match mass budget, some fragments smaller than LB may have been created
%>>> remove from list of fragments, after which, mass budget may not be
%fulfilled
fragments(fragments(:,1)<LB,:) = [];

% mat_frag = [a,ecco,inclo,nodeo,argpo,mo,bstar,mass,radius,...
%             errors,controlled,a_desired,missionlife,constel,...
%             date_created,launch_date,r,v,frag_objectclass,ID_frag]; [num_sat x 24]
[debris1] = func_create_tlesv2_vec(ep, p1_r,p1_v,p1_objclass, fragments, param.max_frag, param.mu, param.req, param.maxID);   % create orbits (OEs); add mass and radius here