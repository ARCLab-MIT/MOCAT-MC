function fix_nodeo_from_csv(folder)
%FIX_NODEO_FROM_CSV Repair the RAAN column of the shipped initial-condition files.
% The .mat files in this folder were built by discos_plus_TLE_historic.m,
% which wrote the mean anomaly into column 4 (idx_nodeo, RAAN) instead of the
% right ascension of the ascending node. The generator is fixed now, but the
% shipped files still carry the wrong column, so this script rewrites column 4
% from the space-track OMM csv files that accompany them.
%
% Every row is matched to its csv record by (NORAD_CAT_ID, ARG_OF_PERICENTER,
% MEAN_ANOMALY); the script stops with an error if any row cannot be matched,
% so it never silently leaves a row untouched.
%
% Usage:  fix_nodeo_from_csv                      % acts on this folder
%         fix_nodeo_from_csv('/path/to/folder')   % acts on given folder
if nargin < 1 || isempty(folder)
    folder = fileparts(mfilename('fullpath'));
end

idx_nodeo = 4; idx_argpo = 5; idx_mo = 6; idx_ID = 24;

files = dir(fullfile(folder,'*.mat'));
for k = 1:numel(files)
    [~,base] = fileparts(files(k).name);
    csvfile = fullfile(folder,[base '.csv']);
    if ~exist(csvfile,'file')
        warning('no csv for %s, skipped', files(k).name);
        continue
    end
    S = load(fullfile(folder,files(k).name));
    t = readtable(csvfile,'TextType','string');

    key_csv = arrayfun(@(a,b,c) sprintf('%d|%.8f|%.8f',a,b,c), ...
        t.NORAD_CAT_ID, round(deg2rad(t.ARG_OF_PERICENTER),8), ...
        round(deg2rad(t.MEAN_ANOMALY),8), 'UniformOutput', false);
    key_mat = arrayfun(@(a,b,c) sprintf('%d|%.8f|%.8f',a,b,c), ...
        S.mat_sats(:,idx_ID), round(S.mat_sats(:,idx_argpo),8), ...
        round(S.mat_sats(:,idx_mo),8), 'UniformOutput', false);

    [tf, loc] = ismember(key_mat, key_csv);
    if ~all(tf)
        error('%s: %i rows could not be matched to %s', files(k).name, sum(~tf), csvfile);
    end

    S.mat_sats(:,idx_nodeo) = deg2rad(t.RA_OF_ASC_NODE(loc));
    save(fullfile(folder,files(k).name),'-struct','S');
    fprintf('%s: rewrote RAAN for %i objects\n', files(k).name, size(S.mat_sats,1));
end
end
