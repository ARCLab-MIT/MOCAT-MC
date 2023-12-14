function [outmatsat,g1,g2,g3] = fillMassRadiusResample(inmatsat,varargin)
% Resample existing data to fill in missing info from DISCOS
% Input: inmatsat -- matsat to be modified
%        g1,g2,g3 -- optional; gX.gm defines the sampling dist (gmdistribution object)
% object classes via objclass2int(1:12,2):
%   P,PMRO,Pfrag,Pdeb,RB,RBMRO,RBfrag,RBdeb,Deb,OtherDeb,Unkwn,untracked

% MATSATS index DEFINITION
getidx;  

% dobj = inmatsat(:,idx_objectclass);
% for stats, combine into:
%    payload {1}, RB {5}, debris {2,3,4,6,7,8,9,10,11}

% fit data, save mean and covariances - if available, use given g1,g2,g3's GM model
[g1,g2,g3] = getZeroGroups(inmatsat);
if nargin > 1 && ~isempty(intersect(fieldnames(varargin{1}),'gm'))
    g1.gm = varargin{1}.gm;  % if supplied, use provided GM models
    g2.gm = varargin{2}.gm;
    g3.gm = varargin{3}.gm;
else  % if no GM provided, calculate
    % g1.gm = gmdistribution(0,0);
    % g2.gm = gmdistribution(0,0);
    % g3.gm = gmdistribution(0,0);
    X = inmatsat(g1.nzno,[idx_radius, idx_mass]);
    if ~isempty(X)
        GMModel = fitgmdist(X,1);  g1.gm = GMModel;
    else
        g1.gm =  gmdistribution([],[]);
    end
    X = inmatsat(g2.nzno,[idx_radius, idx_mass]);
    if ~isempty(X)
        GMModel = fitgmdist(X,1);  g2.gm = GMModel;
    else
        g2.gm = gmdistribution([],[]);
    end
        X = inmatsat(g3.nzno,[idx_radius, idx_mass]);
    if ~isempty(X)
        GMModel = fitgmdist(X,1);  g3.gm = GMModel;
    else
        g3.gm = gmdistribution([],[]);
    end
end

% fill in empty data from mat_sats via sampling
% -- current method: re-sample even if only one parameter is missing (radius or mass)
% -- method: sample 2x N needed, remove any negative entries, randomly choose N
outmatsat = inmatsat;
if 0 ~= union(g1.zm, g1.zr)
    cursamp = random(g1.gm, numel(union(g1.zm, g1.zr))*2);
    cursamp(any(cursamp<=0,2),:) = []; % remove negative entries
    outmatsat(union(g1.zm, g1.zr),[idx_radius,idx_mass]) = cursamp(1:numel(union(g1.zm, g1.zr)),:);
end
if 0 ~= union(g2.zm, g2.zr)
    cursamp = random(g2.gm, numel(union(g2.zm, g2.zr))*2);
    cursamp(any(cursamp<=0,2),:) = []; % remove negative entries
    outmatsat(union(g2.zm, g2.zr),[idx_radius,idx_mass]) = cursamp(1:numel(union(g2.zm, g2.zr)),:);
end
if 0 ~= union(g3.zm, g3.zr)
    cursamp = random(g3.gm, numel(union(g3.zm, g3.zr))*2);
    cursamp(any(cursamp<=0,2),:) = []; % remove negative entries
    outmatsat(union(g3.zm, g3.zr),[idx_radius,idx_mass]) = cursamp(1:numel(union(g3.zm, g3.zr)),:);
end
fprintf('Resampled %i, %i, %i (g1,g2,g3) objects \n', numel(union(g1.zm, g1.zr)), ...
    numel(union(g2.zm, g2.zr)), numel(union(g3.zm, g3.zr)))

return;

end