function v = lininterp1_vec(X, V, x)
% linear interpolation, given set of X and V values, and an x query
% assumes X values are in strictly increasing order
%
% Differences from matlab built-in :
%       much, much faster
%       if coordinate is exactly on the spot, doesn't look at neighbors.  e.g. interpolate([blah, blah2], [0, NaN], blah) returns 0 instead of NaN
%       extends values off the ends instead of giving NaN
%       

if length(X) ~= length(V), error('X and V sizes do not match'); end

if size(V,2)>1
    V = V';
end
if size(X,2)>1
    X = X';
end
if size(x,2)>1
    x = x';
end

% pindexx = zeros(numel(x),1);
% indexx = zeros(numel(x),1);
% slopex = zeros(numel(x),1);
% for i=1:numel(x)
%     pindexx(i) = max([NaN,find((x(i) >= X), 1, 'last')]);
%     indexx(i)  = max([NaN,find((x(i) <= X), 1, 'first')]);
% end

pindexx = discretize(x,X);
indexx = pindexx+1;

if any(isnan(pindexx))
    warning('Input time (x=[%f,%f]) is outside the range described by X(1)=%f and X(end)=%f. Assign X(1) value if x<X(1) and X(end) value if x>X(end).',min(x),max(x),X(1),X(end))
    find_nan = find(isnan(pindexx)); %indices of x outside range
    check_nan_above = x(find_nan)>X(end); %check x values above range
    idx_above = find_nan(check_nan_above); %compute indices above range
    pindexx(idx_above) = length(X)-1; %assign pindexx to lower bound of last bin
    indexx(idx_above) = length(X); %assign indexx to upper bound of last bin
    x(idx_above) = X(end); %assign x to upper bound of last bin
    idx_below = find_nan(~check_nan_above); %compute indices below range
    pindexx(idx_below) = 1; %assign pindexx to lower bound of first bin
    indexx(idx_below) = 2; %assign indexx to upper bound of first bin
    x(idx_below) = X(1); %assign x to lower bound of first bin
end

% isnanpindexx = isnan(pindexx);
% isnanindexx  = isnan(indexx);
% pindexx(isnanpindexx) = indexx(isnanpindexx);
% indexx(isnanindexx) = pindexx(isnanindexx);

% isnotequalx = pindexx~=indexx;
% Xp = X(pindexx(isnotequalx));
% slopex(isnotequalx) = (x(isnotequalx) - Xp) ./ (X(indexx(isnotequalx)) - Xp);

Xp = X(pindexx);
slopex = (x - Xp) ./ (X(indexx) - Xp);


v = V(pindexx) .* (1 - slopex) + V(indexx) .* slopex;
