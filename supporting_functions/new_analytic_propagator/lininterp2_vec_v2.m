function v = lininterp2_vec_v2(X, Y, V, x, y)
% linear interpolation, given set of X, Y, and V values, and an x, y query
% assumes X and Y values are in strictly increasing order
%
% Differences from matlab built-in :
%       order of arguments switched
%       much, much faster
%       if coordinate is exactly on the spot, doesn't look at neighbors.  e.g. interpolate([blah, blah2], [0, NaN], blah) returns 0 instead of NaN
%       extends values off the ends instead of giving NaN
%       

if ((length(X) ~= size(V, 1)) || (length(Y) ~= size(V, 2))), error('[length(X), length(Y)] does not match size(V)'); end

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
% indexx(indexx>length(X)) = NaN;

if any(isnan(pindexx))
    error('Input altitude (x=[%f,%f]) is outside the range described by X(1)=%f and X(end)=%f. This should not have happened since analytical_propagator already checks for limits.',min(x),max(x),X(1),X(end))
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

if size(Y,2)>1
    Y = Y';
end
if size(y,2)>1
    y = y';
end
% pindexy = zeros(numel(y),1);
% indexy = zeros(numel(y),1);
% slopey = zeros(numel(y),1);
% for i=1:numel(y)
%     pindexy(i) = max([NaN,find((y(i) >= Y), 1, 'last')]);
%     indexy(i)  = max([NaN,find((y(i) <= Y), 1, 'first')]);
% end

pindexy = discretize(y,Y);
indexy = pindexy+1;
% indexy(indexy>length(Y)) = NaN;

if any(isnan(pindexy))
    warning('Input time (y=[%f,%f]) is outside the range described by Y(1)=%f and Y(end)=%f. Assign Y(1) value if y<Y(1) and Y(end) value if y>Y(end).',min(y),max(y),Y(1),Y(end))
    find_nan = find(isnan(pindexy)); %indices of y outside range
    check_nan_above = y(find_nan)>Y(end); %check y values above range
    idx_above = find_nan(check_nan_above); %compute indices above range
    pindexy(idx_above) = length(Y)-1; %assign pindexy to lower bound of last bin
    indexy(idx_above) = length(Y); %assign indexy to upper bound of last bin
    y(idx_above) = Y(end); %assign y to upper bound of last bin
    idx_below = find_nan(~check_nan_above); %compute indices below range
    pindexy(idx_below) = 1; %assign pindexy to lower bound of first bin
    indexy(idx_below) = 2; %assign indexy to upper bound of first bin
    y(idx_below) = Y(1); %assign y to lower bound of first bin
end

% isnanpindexy = isnan(pindexy);
% isnanindexy  = isnan(indexy);
% pindexy(isnanpindexy) = indexy(isnanpindexy);
% indexy(isnanindexy) = pindexy(isnanindexy);
% 
% isnotequaly = pindexy~=indexy;
% Yp = Y(pindexy(isnotequaly));
% slopey(isnotequaly) = (y(isnotequaly) - Yp) ./ (Y(indexy(isnotequaly)) - Yp);

Yp = Y(pindexy);
slopey = (y - Yp) ./ (Y(indexy) - Yp);

num_rows = size(V,1);
v = V(pindexx+(pindexy-1)*num_rows) .* (1 - slopex) .* (1 - slopey) + V(indexx+(pindexy-1)*num_rows) .* slopex .* (1 - slopey) ... 
  + V(pindexx+(indexy-1)*num_rows)  .* (1 - slopex) .*      slopey  + V(indexx+(indexy-1)*num_rows)  .* slopex .*      slopey;
