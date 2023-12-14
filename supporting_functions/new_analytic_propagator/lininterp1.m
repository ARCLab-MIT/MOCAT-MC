function v = lininterp1(X, V, x)
% linear interpolation, given set of X and V values, and an x query
% assumes X values are in strictly increasing order
%
% Differences from matlab built-in :
%       much, much faster
%       if coordinate is exactly on the spot, doesn't look at neighbors.  e.g. interpolate([blah, blah2], [0, NaN], blah) returns 0 instead of NaN
%       extends values off the ends instead of giving NaN
%       

if length(X) ~= length(V), error('X and V sizes do not match'); end

pindex = find((x >= X), 1, 'last');
index = find((x <= X), 1, 'first');

if isempty(pindex)
    warning('interpolating before beginning');
    pindex = index;
    slope = 0;
elseif isempty(index)
    warning('interpolating after end');
    index = pindex;
    slope = 0;
elseif pindex == index
    slope = 0;
else
    Xp = X(pindex);
    slope = (x - Xp) / (X(index) - Xp);
end

v = V(pindex) * (1 - slope) + V(index) * slope;

end

