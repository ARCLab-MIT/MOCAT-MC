%% Merge TLE and DISCOS data - for historic TLEs (YYYY.csv)
% INSTEAD of SATCAT, just mark payloads < 8 yrs old as active...

% based on ../disco_plus_TLE.m

%% TLE data
% Treat as truth data for: NORAD ID, OE, B*, launch date, object type

% Space-track.org TLE as CSV file from query:
% https://www.space-track.org/basicspacedata/query/class/gp_history/EPOCH/2023-01-01--2023-01-03/DECAY_DATE/null-val/MEAN_MOTION/>3/orderby/NORAD_CAT_ID/format/csv
% NOTE: any more than 3 days result in error to download (too big a file?)
%       Above query: history of TLE's (GP) for epoch between 1/1/23 - 1/3/23, active sats, mean motion > 3 revs/day
% API info: https://www.space-track.org/documentation#/api
% API info on gp_history object: https://www.space-track.org/basicspacedata/modeldef/class/gp_history/format/html
% Data definition follows CCSDS Recommended Standard 502.0-B-2: https://public.ccsds.org/Pubs/502x0b2c1e2.pdf

tic
csvs = dir('*.csv'); 

load d_2023.mat;         % Discosweb Object data (up to Mar 2023)

fprintf('DISCOS data number of objects: %i \n', size(d,1))
fprintf('DISCOS data headers: \n')
disp(fieldnames(d(1).attributes))

% relevant fields: satno, objectClass, mass, shape, diameter, span, height, width

radiusearthkm = 6378.137;

getidx;

% % MATSATS DEFINITION
% idx_a = 1; idx_ecco = 2; idx_inclo = 3;
% idx_nodeo = 4; idx_argpo = 5; idx_mo = 6;
% idx_bstar = 7; idx_mass = 8; idx_radius = 9;
% idx_error = 10; idx_controlled = 11; idx_a_desired = 12;   % [0,0,0]
% idx_missionlife = 13; idx_constel = 14; idx_date_created = 15;  % [0,0,nan]
% idx_launch_date = 16; idx_r = [17 18 19]; idx_v = [20 21 22]; % date: JD,  RV:0
% idx_objectclass = 23; idx_ID = 24;


for fileind = 16:numel(csvs)
    tlename = csvs(fileind).name;
    opts = detectImportOptions(tlename);
    mattles = readtable(tlename,opts);          % includes multiple repeated objects
    [~,lastinds,~] = unique(mattles.NORAD_CAT_ID,'last');
    mattles = mattles(lastinds,:);              % overwrite

    fprintf('TLE source:  %s \n', tlename);
    fprintf('TLE data number of unique objects: %i \n', size(mattles,1))
%     fprintf('TLE data CSV headers: \n')
%     disp(opts.VariableNames)
    satIDs = mattles.NORAD_CAT_ID;
    % filter data (remove temporary sat IDs)
    fprintf('%i of %i (%0.1f%%) TLE satnos are > 80000 (temporary designations); removing... \n',...
        sum(satIDs>80000), numel(satIDs), sum(satIDs>80000)/numel(satIDs)*100 );
    mattles(mattles.NORAD_CAT_ID > 80000, :) = [];
    satIDs = mattles.NORAD_CAT_ID;


    % DISCOS data & filter (match via NORAD id etc)
    % some don't have some fields, including satno, so clean up data files
    % Deference given to TLEs (e.g. objects from DISCOS without NORAD_ID is removed; e.g. vimpel data)
    
%     % clear all; 
%     % load d_all.mat;         % Discosweb Object data (up to Aug 2022)
%     load d_2023.mat;         % Discosweb Object data (up to Mar 2023)
%     
%     fprintf('DISCOS data number of objects: %i \n', size(d,1))
%     fprintf('DISCOS data headers: \n')
%     disp(fieldnames(d(1).attributes))
%     
%     % relevant fields: satno, objectClass, mass, shape, diameter, span, height, width
    
    dc = struct2cell(d); 
    dc = dc(4,:);  % d.attributes
    ysatno = cellfun(@(x) ~isempty(x.satno), dc);  % 0: no satno; 1: yes satno
    fprintf('%i / %i DISCOS data has NORAD satellite number (%0.1f%%) \n ',sum(ysatno), numel(dc), sum(ysatno)/numel(dc)*100);
    fprintf('Removing %i DISCOS objects without NORAD IDs \n', numel(dc) - sum(ysatno));
    dd = dc(ysatno);
    ddsatnos = cellfun(@(x) x.satno, dd);
    
    % plot to show repeated satnos in discosweb data
    % figure; plot(diff(sort(ddsatnos)))  
     [~, i, j] = unique(ddsatnos,'first');
     indexToDupes = find(not(ismember(1:numel(ddsatnos),i)));
    % for ind = 1:numel(indexToDupes)
    %     ii = indexToDupes(ind);
    %     iisatno = ddsatnos(ii);
    %     dupinds = find(ddsatnos == iisatno);
    %     isequal(dd{dupinds(1)},dd{dupinds(2)}) % Test to show if data is an identical repeat; 
    %                                            % Answer: yes, all duplicates are exact duplicates
    % end
    
    fprintf('%i duplicate entries found in discos data (with satnums); removing... \n', numel(indexToDupes));
    
    dd(indexToDupes) = [];                      % remove
    ddsatnos = cellfun(@(x) x.satno, dd);       % redo satno extraction
    
    
    % TLE data (mattles)
    
    % CHECK IF TLE SATIDS ARE ALL REPRESENTED IN DISCOS DDSATNOS
    % ANSWER: NO -- some temporary IDs are assigned in TLEs (satno > 80,000, 6 digits, etc)...
    % for ind = 1:numel(satIDs)
    %     if isempty(find(ddsatnos == satIDs(ind),1))
    %         disp(satIDs(ind))
    %     end
    % end
    
    % GET MASS, RADIUS, OBJTYPE FROM DISCOS DATA to match with TLE's satIDs
    dmass = zeros(numel(satIDs),1);
    dradius = zeros(numel(satIDs),1);
    dobj = nan(numel(satIDs),1);
    missingSatnos = [];
    
    % Which discosweb data field is populated?
    [~,ia,ib] = intersect(satIDs, ddsatnos);
    indxsect = cellfun(@(x) ~isempty(x.xSectAvg), dd(ib));
    inddiameter = cellfun(@(x) ~isempty(x.diameter), dd(ib));
    indspan = cellfun(@(x) ~isempty(x.span), dd(ib));
    indwidth = cellfun(@(x) ~isempty(x.width), dd(ib));
    inddepth = cellfun(@(x) ~isempty(x.depth), dd(ib));
    indheight = cellfun(@(x) ~isempty(x.height), dd(ib));
    indmass = cellfun(@(x) ~isempty(x.mass), dd(ib));
    
    
    fprintf('Percentage of NORAD Discos objects (%i objects) with following attributes: \n', numel(ib))
    fprintf(['\txSectAvg: %0.1f%% \n\tdiameter: %0.1f%% \n\tspan: %0.1f%% \n\t' ...
        'width: %0.1f%%\n\tdepth: %0.1f%%\n\theight: %0.1f%%\n\tmass: %0.1f%%\n '], ...
        sum(indxsect)/numel(ib)*100,...
        sum(inddiameter)/numel(ib)*100,...
        sum(indspan)/numel(ib)*100,...
        sum(indwidth)/numel(ib)*100,...
        sum(inddepth)/numel(ib)*100,...
        sum(indheight)/numel(ib)*100,...
        sum(indmass)/numel(ib)*100);
    
    indNoSize = ~(indxsect | inddiameter | indspan | indwidth | inddepth | indheight);
    fprintf('\n Percentage without any size info: %0.1f%%\n', ...
        sum(indNoSize)/numel(ib)*100);
    fprintf('\n Percentage without any size NOR mass: %0.1f%%\n', ...
        sum(indNoSize & ~indmass)/numel(ib)*100);
    
    % ddsatnos(ib & indNoSize');  % list of SCC numbers of objects w/o discos size info
%     figure; histogram(dobj(indNoSize))
    
    for ind = 1:numel(satIDs)
        dind = find(ddsatnos == satIDs(ind));
        if ~isempty(dind)
            if ~isempty(dd{dind}.mass)
                dmass(ind) = dd{dind}.mass;
            end
            if ~isempty(dd{dind}.xSectAvg)
                dradius(ind) = sqrt(dd{dind}.xSectAvg/pi);  
                        % CHOOSE DISCOS RADIUS HERE ^
                        % options: depth,xSectAvg,width,xSectMax,height,diameter,span)
            end
            if ~isempty(dd{dind}.objectClass)
                try
                    dobj(ind) = objclass2int(dd{dind}.objectClass,1);
                catch
                    warning('NORAD %i: %s --> Other Debris (10)', ...
                        satIDs(ind), dd{dind}.objectClass);
                    dobj(ind) = 10;
                end
            end
        elseif satIDs(ind) == 53239 || satIDs(ind) == 54216
            dobj(ind) = 1;  % above two are Chinese payloads
        else
            missingSatnos = [missingSatnos, satIDs(ind)];
        end
    end
    
    
    
    fprintf('%i of %i (%0.1f%%) TLE objects are missing from DISCOS \n', ...
        numel(missingSatnos), numel(satIDs), numel(missingSatnos)/numel(satIDs)*100);

    % CREATE COMBINED INITIALIZED MATSATS
    
    % load('initialized.mat','mat_sats')            % to compare
    % mat_sats(1,:)
    
 
    zs = zeros(numel(satIDs),1);
    ns = nan(numel(satIDs),1);
    
    matsatsout = [mattles.SEMIMAJOR_AXIS/radiusearthkm, mattles.ECCENTRICITY,  ... % col 1:2
                    deg2rad(mattles.INCLINATION), deg2rad(mattles.MEAN_ANOMALY), ... % col 3:4
                    deg2rad(mattles.ARG_OF_PERICENTER), deg2rad(mattles.MEAN_ANOMALY), ... % col 5:6
                    mattles.BSTAR, dmass, dradius, ...                      % col 7:9
                    zs, zs, zs, zs, zs, ns, ...                             % col 10:15
                    juliandate(mattles.LAUNCH_DATE), zs, zs, zs, ...        % col 16:19
                    zs, zs, zs, dobj, mattles.NORAD_CAT_ID]; % col 20:24
    
    % NOTE: TLE'S EPOCH FOR OE IS NOT SAVED NOR SHIFTED TO BE IDENTICAL!! 
    
    % SANITY CHECK:
    %     mat_sats(1,:)
    %     matsatsout(1,:)
    
    mat_sats = matsatsout;
    savename = sprintf('%s.mat',tlename(1:end-4));
    save(savename,"mat_sats");
    fprintf('SAVED: %S \n',savename);
end
toc

%% HIST OF A FILE 
% load 2018.mat

% idx_a = 1; idx_ecco = 2; idx_inclo = 3; idx_nodeo = 4; idx_argpo = 5; idx_mo = 6; 
% idx_bstar = 7; idx_mass = 8; idx_radius = 9;
%     idx_error = 10; 

figure(10); clf;
subplot(3,5,1); % a,e
plot((mat_sats(:,idx_a)-1)*radiusearthkm, mat_sats(:,idx_ecco),'.'); 
xlabel('Alt (km)'); ylabel('ecc');
subplot(3,5,2); % inc
histogram(rad2deg(mat_sats(:,idx_inclo)));
xlabel('Inclination');
subplot(3,5,3); % nodeo
histogram(rad2deg(mat_sats(:,idx_nodeo)));
xlabel('nodeo')
subplot(3,5,4); % argpo
histogram(rad2deg(mat_sats(:,idx_argpo)));
xlabel('argpo')
subplot(3,5,5); % mo
histogram(rad2deg(mat_sats(:,idx_mo)));
xlabel('mo')
subplot(3,5,6); % b*
histogram(mat_sats(:,idx_bstar));
xlabel('bstar'); a = gca; a.YScale = 'log';
subplot(3,5,7); % mass
histogram(mat_sats(:,idx_mass));
xlabel('mass'); a = gca; a.YScale = 'log';
xlim(prctile(mat_sats(:,idx_mass),[0,99.95]));
subplot(3,5,8); % radius
histogram(mat_sats(:,idx_radius));
xlabel('radius'); a = gca; a.YScale = 'log';
xlim(prctile(mat_sats(:,idx_radius),[0,99.95]));

% idx_controlled = 11; idx_a_desired = 12; 
% idx_missionlife = 13; idx_constel = 14; idx_date_created = 15; idx_launch_date = 16;
%     idx_r = [17 18 19]; idx_v = [20 21 22]; idx_objectclass = 23; idx_ID = 24;
subplot(3,5,9); % controlled
histogram(mat_sats(:,idx_controlled));
xlabel('controlled')
subplot(3,5,10); % idx_a_desired of controlled sats
histogram(mat_sats(logical(mat_sats(:,idx_controlled)), idx_a_desired));
xlabel('Desired SMA of controlled sats')
subplot(3,5,11); % controlled
histogram(mat_sats(:,idx_launch_date));
xlabel('launch date')
subplot(3,5,12); % r vs v
r = mat_sats(:,idx_r);
v = mat_sats(:,idx_v);
plot(sqrt(r(:,1).^2 + r(:,2).^2 + r(:,3).^2), sqrt(v(:,1).^2 + v(:,2).^2 + v(:,3).^2),'.');
xlabel('r (km)'); ylabel('v (km/s)')
subplot(3,5,13); % obj class
histogram(mat_sats(:,idx_objectclass),[1:11]);
% a = gca; a.XTick = 1:11;  a.XTickLabel = objclass2int(1:11,2);
xlabel('object class');
subplot(3,5,14); % obj class
histogram(mat_sats(:,idx_ID));
xlabel('ID')










%% FROM PREVIOUS (INDIV) VERSION -- ../discos_plus_TLE.m
tlename = 'spacetrack_TLE_01_01-03_2023_MEO.csv';  
opts = detectImportOptions(tlename); 
mattles = readtable(tlename,opts);          % includes multiple repeated objects

% Filter to grab only the latest TLE per object
[~,lastinds,~] = unique(mattles.NORAD_CAT_ID,'last');
mattles = mattles(lastinds,:);              % overwrite

fprintf('TLE source:  %s \n', tlename);
fprintf('TLE data number of unique objects: %i \n', size(mattles,1))
fprintf('TLE data CSV headers: \n')
disp(opts.VariableNames)

satIDs = mattles.NORAD_CAT_ID;
figure;plot(satIDs,'x')

% filter data (remove temporary sat IDs)
fprintf('%i of %i (%0.1f%%) TLE satnos are > 80000 (temporary designations); removing... \n',...
    sum(satIDs>80000), numel(satIDs), sum(satIDs>80000)/numel(satIDs)*100 );
mattles(mattles.NORAD_CAT_ID > 80000, :) = [];
satIDs = mattles.NORAD_CAT_ID;

% histogram of B*
figure; h = histogram(mattles.BSTAR,[-0.005:0.0001:0.005]);
h.EdgeAlpha=0;  xlabel('B*')

%% Load SATCAT data for active/inactive flag
% https://celestrak.org/satcat/search.php
% https://celestrak.org/satcat/status.php

SATCATname = 'satcat_03_13_2023.csv';
opts = detectImportOptions(SATCATname); 
satcatdata = readtable(SATCATname,opts);

% fields(satcatdata)
% satcatdata.NORAD_CAT_ID
% satcatdata.OPS_STATUS_CODE
% satcatdata.OBJECT_TYPE

% 2D histogram to compare ObjectType to Status code
xcodes = unique(satcatdata.OPS_STATUS_CODE);
ytype = unique(satcatdata.OBJECT_TYPE);
X = zeros(size(satcatdata,1),1);
for ind = 1:numel(xcodes)
    X(cellfun(@(x) strcmp(x,xcodes{ind}), satcatdata.OPS_STATUS_CODE)) = ind;
end
Y = zeros(size(satcatdata,1),1);
for ind = 1:numel(ytype)
    Y(cellfun(@(x) strcmp(x,ytype{ind}), satcatdata.OBJECT_TYPE)) = ind;
end

% separate histograms
% figure; subplot(211); histogram(X,'binmethod','integers'); a = gca; a.XTickLabel = xcodes;
% title('Active status from SATCAT')
% subplot(212); histogram(Y,'binmethod','integers'); a = gca; 
% a.XTick = 1:numel(ytype);  a.XTickLabel = ytype;
% title('Object Type from SATCAT')

% 2D histograms
figure; histogram2(X,Y,'binmethod','integers','displaystyle','bar',...
    'ShowEmptyBins','off','FaceColor','flat');
a = gca; a.YTick = 1:numel(ytype);  a.YTickLabel = ytype;
a.XTick = 1:numel(xcodes);  a.XTickLabel = xcodes;
ylabel('Object type'); xlabel('Active flag')
colorbar; grid off;

% 2D simplified
% Active is any satellite with an operational status: +, P, B, S, X.
% non active: unknown, -, D
% Also, just look at the relevant TLE's (LEO, non-decayed etc) -- satIDs
[~,isatcat,~] = intersect(satcatdata.NORAD_CAT_ID, satIDs);
Xs = zeros(numel(isatcat),1);
for ind = 1:numel(xcodes)
    Xs(cellfun(@(x) strcmp(x,xcodes{ind}), satcatdata.OPS_STATUS_CODE(isatcat))) = ind;
end
Xs(Xs == 4 | Xs == 6 | Xs == 4 | Xs >= 7) = 2;  % change B,P,S,X to '+'
Xs(Xs == 5) = 4;                                    % move D down a slot
xcodess = {'?','Active','Non-active','Decayed'};
Ys = zeros(numel(isatcat),1);
for ind = 1:numel(ytype)
    Ys(cellfun(@(x) strcmp(x,ytype{ind}), satcatdata.OBJECT_TYPE(isatcat))) = ind;
end
figure; histogram2(Xs,Ys,'binmethod','integers','displaystyle','bar',...
    'ShowEmptyBins','off','FaceColor','flat');
a = gca; a.YTick = 1:numel(ytype);  a.YTickLabel = ytype;
a.XTick = 1:numel(xcodess);  a.XTickLabel = xcodess;
ylabel('Object type'); xlabel('')
colorbar; grid off;


%% Load SATCAT (space-track) data for RCS category
% https://www.space-track.org/basicspacedata/query/class/SATCAT/orderby/NORAD_CAT_ID/format/csv

% This code is here just for analysis; matsat is not affected
% Changes matsats is done in setupTLE for ESA fill-in case


SATCATstname = 'spacetrack_satcat_03_2023.csv';
opts = detectImportOptions(SATCATstname); 
satcatstdata = readtable(SATCATstname,opts);

nlg = cellfun(@(x) strcmp(x,'LARGE'), satcatstdata.RCS_SIZE);
nmd = cellfun(@(x) strcmp(x,'MEDIUM'), satcatstdata.RCS_SIZE);
nsm = cellfun(@(x) strcmp(x,'SMALL'), satcatstdata.RCS_SIZE);

fprintf('SATCAT from Space-track has %i small, %i medium, %i large objects\n',...
    sum(nsm),sum(nmd),sum(nlg));


% comparison between large/medium/small RCSs and object type

% 2D histogram to compare ObjectType to Status code
xcodes = unique(satcatstdata.RCS_SIZE);
ytype = unique(satcatstdata.OBJECT_TYPE);

X = zeros(size(satcatstdata,1),1);
for ind = 1:numel(xcodes)
    X(cellfun(@(x) strcmp(x,xcodes{ind}), satcatstdata.RCS_SIZE)) = ind;
end
Y = zeros(size(satcatstdata,1),1);
for ind = 1:numel(ytype)
    Y(cellfun(@(x) strcmp(x,ytype{ind}), satcatstdata.OBJECT_TYPE)) = ind;
end
figure; histogram2(X,Y,'binmethod','integers','displaystyle','bar',...
    'ShowEmptyBins','off','FaceColor','flat');
a = gca; a.YTick = 1:numel(ytype);  a.YTickLabel = ytype;
a.XTick = 1:numel(xcodes);  a.XTickLabel = xcodes;
ylabel('Object type'); xlabel('RCS Size')
colorbar; grid off;
title('Space-track SATCAT RCS category vs object type')

%% DISCOS data & filter (match via NORAD id etc)
% some don't have some fields, including satno, so clean up data files
% Deference given to TLEs (e.g. objects from DISCOS without NORAD_ID is removed; e.g. vimpel data)

% clear all; 
% load d_all.mat;         % Discosweb Object data (up to Aug 2022)
load d_2023.mat;         % Discosweb Object data (up to Mar 2023)

fprintf('DISCOS data number of objects: %i \n', size(d,1))
fprintf('DISCOS data headers: \n')
disp(fieldnames(d(1).attributes))

% relevant fields: satno, objectClass, mass, shape, diameter, span, height, width

dc = struct2cell(d); 
dc = dc(4,:);  % d.attributes
ysatno = cellfun(@(x) ~isempty(x.satno), dc);  % 0: no satno; 1: yes satno
fprintf('%i / %i DISCOS data has NORAD satellite number (%0.1f%%) \n ',sum(ysatno), numel(dc), sum(ysatno)/numel(dc)*100);
fprintf('Removing %i DISCOS objects without NORAD IDs \n', numel(dc) - sum(ysatno));
dd = dc(ysatno);
ddsatnos = cellfun(@(x) x.satno, dd);

% plot to show repeated satnos in discosweb data
% figure; plot(diff(sort(ddsatnos)))  
 [~, i, j] = unique(ddsatnos,'first');
 indexToDupes = find(not(ismember(1:numel(ddsatnos),i)));
% for ind = 1:numel(indexToDupes)
%     ii = indexToDupes(ind);
%     iisatno = ddsatnos(ii);
%     dupinds = find(ddsatnos == iisatno);
%     isequal(dd{dupinds(1)},dd{dupinds(2)}) % Test to show if data is an identical repeat; 
%                                            % Answer: yes, all duplicates are exact duplicates
% end

fprintf('%i duplicate entries found in discos data (with satnums); removing... \n', numel(indexToDupes));

dd(indexToDupes) = [];                      % remove
ddsatnos = cellfun(@(x) x.satno, dd);       % redo satno extraction


%% TLE data (mattles)

% CHECK IF TLE SATIDS ARE ALL REPRESENTED IN DISCOS DDSATNOS
% ANSWER: NO -- some temporary IDs are assigned in TLEs (satno > 80,000, 6 digits, etc)...
% for ind = 1:numel(satIDs)
%     if isempty(find(ddsatnos == satIDs(ind),1))
%         disp(satIDs(ind))
%     end
% end

% GET MASS, RADIUS, OBJTYPE FROM DISCOS DATA to match with TLE's satIDs
dmass = zeros(numel(satIDs),1);
dradius = zeros(numel(satIDs),1);
dobj = nan(numel(satIDs),1);
missingSatnos = [];

% Which discosweb data field is populated?
[~,ia,ib] = intersect(satIDs, ddsatnos);
indxsect = cellfun(@(x) ~isempty(x.xSectAvg), dd(ib));
inddiameter = cellfun(@(x) ~isempty(x.diameter), dd(ib));
indspan = cellfun(@(x) ~isempty(x.span), dd(ib));
indwidth = cellfun(@(x) ~isempty(x.width), dd(ib));
inddepth = cellfun(@(x) ~isempty(x.depth), dd(ib));
indheight = cellfun(@(x) ~isempty(x.height), dd(ib));
indmass = cellfun(@(x) ~isempty(x.mass), dd(ib));


fprintf('Percentage of NORAD Discos objects (%i objects) with following attributes: \n', numel(ib))
fprintf(['\txSectAvg: %0.1f%% \n\tdiameter: %0.1f%% \n\tspan: %0.1f%% \n\t' ...
    'width: %0.1f%%\n\tdepth: %0.1f%%\n\theight: %0.1f%%\n\tmass: %0.1f%%\n '], ...
    sum(indxsect)/numel(ib)*100,...
    sum(inddiameter)/numel(ib)*100,...
    sum(indspan)/numel(ib)*100,...
    sum(indwidth)/numel(ib)*100,...
    sum(inddepth)/numel(ib)*100,...
    sum(indheight)/numel(ib)*100,...
    sum(indmass)/numel(ib)*100);

indNoSize = ~(indxsect | inddiameter | indspan | indwidth | inddepth | indheight);
fprintf('\n Percentage without any size info: %0.1f%%\n', ...
    sum(indNoSize)/numel(ib)*100);
fprintf('\n Percentage without any size NOR mass: %0.1f%%\n', ...
    sum(indNoSize & ~indmass)/numel(ib)*100);

% ddsatnos(ib & indNoSize');  % list of SCC numbers of objects w/o discos size info
figure; histogram(dobj(indNoSize))

for ind = 1:numel(satIDs)
    dind = find(ddsatnos == satIDs(ind));
    if ~isempty(dind)
        if ~isempty(dd{dind}.mass)
            dmass(ind) = dd{dind}.mass;
        end
        if ~isempty(dd{dind}.xSectAvg)
            dradius(ind) = sqrt(dd{dind}.xSectAvg/pi);  
                    % CHOOSE DISCOS RADIUS HERE ^
                    % options: depth,xSectAvg,width,xSectMax,height,diameter,span)
        end
        if ~isempty(dd{dind}.objectClass)
            dobj(ind) = objclass2int(dd{dind}.objectClass,1);
        end
    elseif satIDs(ind) == 53239 || satIDs(ind) == 54216
        dobj(ind) = 1;  % above two are Chinese payloads
    else
        missingSatnos = [missingSatnos, satIDs(ind)];
    end
end



fprintf('%i of %i (%0.1f%%) TLE objects are missing from DISCOS \n', ...
    numel(missingSatnos), numel(satIDs), numel(missingSatnos)/numel(satIDs)*100);

% Data histograms
figure;   % MASS
histogram(dmass); xlabel('kg')
a = gca; a.XScale = 'log';

figure;   % RADIUS
w = cellfun(@(x) x.width, dd(ib(indwidth)));
h = histogram(w); xlabel('m'); 
h.EdgeAlpha = 0;
a = gca; a.YScale = 'log';

figure;  % ZOOMED IN
xs = cellfun(@(x) x.xSectAvg, dd(ib(indxsect)));
di = cellfun(@(x) x.diameter, dd(ib(inddiameter)));
sp = cellfun(@(x) x.span, dd(ib(indspan)));
w = cellfun(@(x) x.width, dd(ib(indwidth)));
dep = cellfun(@(x) x.depth, dd(ib(inddepth)));
ht = cellfun(@(x) x.height, dd(ib(indheight)));
ms = cellfun(@(x) x.mass, dd(ib(indmass)));
fprintf('Minimum features: xSect radius: %0.3f, diam: %0.3f, span: %0.3f, wid: %0.3f, dep: %0.3f, height: %0.3f, mass: %0.3f\n',...
    sqrt(min(xs) / pi), min(di), min(sp), min(w), min(dep), min(ht), min(ms))

h = histogram(ht,[0:0.1:5]); xlabel('m'); 
h.EdgeAlpha = 0;
% a = gca; a.YScale = 'log';

%% QUICK COMPARISON: EXISTING FULL DATA PER OBJECT TYPE
% dobj  vs  mattles.BSTAR, dmass, dradius
figure;
uobj = unique(dobj);
for oind = 1:numel(uobj)
    oo = uobj(oind);
    subplot(numel(unique(dobj)),1,find(oo==unique(dobj)));
    histogram(dradius(oo == dobj));
end



%% CREATE COMBINED INITIALIZED MATSATS

% load('initialized.mat','mat_sats')            % to compare
% mat_sats(1,:)

% radiusearthkm = cfg.radiusearthkm; 
radiusearthkm = 6378.137;

    % MATSATS DEFINITION
    idx_a = 1; idx_ecco = 2; idx_inclo = 3; 
    idx_nodeo = 4; idx_argpo = 5; idx_mo = 6; 
    idx_bstar = 7; idx_mass = 8; idx_radius = 9;
    idx_error = 10; idx_controlled = 11; idx_a_desired = 12;   % [0,0,0]
    idx_missionlife = 13; idx_constel = 14; idx_date_created = 15;  % [0,0,nan]
    idx_launch_date = 16; idx_r = [17 18 19]; idx_v = [20 21 22]; % date: JD,  RV:0
    idx_objectclass = 23; idx_ID = 24;
    
zs = zeros(numel(satIDs),1);
ns = nan(numel(satIDs),1);

matsatsout = [mattles.SEMIMAJOR_AXIS/radiusearthkm, mattles.ECCENTRICITY,  ... % col 1:2
                deg2rad(mattles.INCLINATION), deg2rad(mattles.MEAN_ANOMALY), ... % col 3:4
                deg2rad(mattles.ARG_OF_PERICENTER), deg2rad(mattles.MEAN_ANOMALY), ... % col 5:6
                mattles.BSTAR, dmass, dradius, ...                      % col 7:9
                zs, zs, zs, zs, zs, ns, ...                             % col 10:15
                juliandate(mattles.LAUNCH_DATE), zs, zs, zs, ...        % col 16:19
                zs, zs, zs, dobj, mattles.NORAD_CAT_ID]; % col 20:24

% NOTE: TLE'S EPOCH FOR OE IS NOT SAVED NOR SHIFTED TO BE IDENTICAL!! 

% SANITY CHECK:
%     mat_sats(1,:)
%     matsatsout(1,:)

mat_sats = matsatsout;
% save('initialized_01-2023_MEO.mat',"mat_sats");

%% matsat object class vs number
% also, look at radius vs mass
figure;
histogram(mat_sats(:,idx_objectclass),1:11)
a = gca;
a.XTickLabelMode = 'manual';
a.XTickLabel = objclass2int(1:11,2);
a.XTick = 1.5:11.5;
a.XTickLabelRotation = 45;
title('Number of object types in matsats'); 

figure;  % existing data
for ii = 1:12
    subplot(3,4,ii);
    msinds = mat_sats(:,idx_objectclass)==ii;
    plot(mat_sats(msinds,idx_radius),mat_sats(msinds,idx_mass),'x');
    title(sprintf('%s (%i)',objclass2int(ii,2),ii));
end

