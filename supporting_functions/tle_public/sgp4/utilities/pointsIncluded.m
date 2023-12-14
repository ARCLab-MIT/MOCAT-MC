function Points=pointsIncluded(xA,yA,polyRegion)
%
%    yA,xA      =   Array coordinates (column, row) physical units
%
%    polyRegion =   A 2Xnpts array of ordered coordinates that define a closed
%               polygon
%OUTPUT
%    Points =   Subset of points in Array that lie within the polygon or
%               on its boundaries
%

%     Lillian Chu
%     Vista Research, Inc.
%     December 2006
%
[n1,n2]=size(polyRegion);
if n2~=2
    error('PolyRegion Error')
elseif n1<3
    error('nodes < 3')
else
    nodes=n1;
end

M = length(yA);
N = length(xA);

minX = min(xA);
maxX = max(xA);
deltaX = (maxX-minX)/(N-1);

minY = min(yA);
maxY = max(yA);
deltaY = (maxY-minY)/(M-1);

X = (polyRegion(:,1)-minX)/deltaX + 1;
Y = (polyRegion(:,2)-minY)/deltaY + 1;

mask = poly2mask(X,Y,M,N);
Points = find(mask==1);

return