function v = lininterp2_vec(X, Y, V, x, y)
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
pindexx = zeros(numel(x),1);
indexx = zeros(numel(x),1);
slopex = zeros(numel(x),1);
for i=1:numel(x)
    pindexx(i) = max([NaN,find((x(i) >= X), 1, 'last')]);
    indexx(i)  = max([NaN,find((x(i) <= X), 1, 'first')]);
end

isnanpindexx = isnan(pindexx);
isnanindexx  = isnan(indexx);
pindexx(isnanpindexx) = indexx(isnanpindexx);
indexx(isnanindexx) = pindexx(isnanindexx);

isnotequalx = pindexx~=indexx;
Xp = X(pindexx(isnotequalx));
slopex(isnotequalx) = (x(isnotequalx) - Xp) ./ (X(indexx(isnotequalx)) - Xp);

if size(Y,2)>1
    Y = Y';
end
if size(y,2)>1
    y = y';
end
pindexy = zeros(numel(y),1);
indexy = zeros(numel(y),1);
slopey = zeros(numel(y),1);
for i=1:numel(y)
    pindexy(i) = max([NaN,find((y(i) >= Y), 1, 'last')]);
    indexy(i)  = max([NaN,find((y(i) <= Y), 1, 'first')]);
end

isnanpindexy = isnan(pindexy);
isnanindexy  = isnan(indexy);
pindexy(isnanpindexy) = indexy(isnanpindexy);
indexy(isnanindexy) = pindexy(isnanindexy);

isnotequaly = pindexy~=indexy;
Yp = Y(pindexy(isnotequaly));
slopey(isnotequaly) = (y(isnotequaly) - Yp) ./ (Y(indexy(isnotequaly)) - Yp);

num_rows = size(V,1);
v = V(pindexx+(pindexy-1)*num_rows) .* (1 - slopex) .* (1 - slopey) + V(indexx+(pindexy-1)*num_rows) .* slopex .* (1 - slopey) ... 
  + V(pindexx+(indexy-1)*num_rows)  .* (1 - slopex) .*      slopey  + V(indexx+(indexy-1)*num_rows)  .* slopex .*      slopey;
