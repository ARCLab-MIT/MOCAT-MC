function  polyregion_add(polyRegion,ifig,varargin)
%         polyregion_add(polyRegion,ifig,iR,Col)
%
%    Add polyRegion outline to current image
%
%    polyRegion =   A 2Xnpts array of ordered coordinates that define a closed
%                   polygon
%    ifig       =   figure handel
%    iR         =   1  => Upper
%               =   2  => Lower region
%               =   [] => Upper & Lower
%    Col        =   Line color/style string
%OUTPUT
%    Plot 
%
figure(ifig);
Col='w';
if isempty(varargin)
    iR=[];
elseif length(varargin)==1
    iR=varargin{1};
else
    iR=varargin{1};
    Col=varargin{2};
end

[n1,n2]=size(polyRegion);
if n2~=2
    error('PolyRegion Error')
elseif n1<3
    error('nodes < 3')
else
    nodes=n1;
end
%Remove redundant node
polyRegion=polyRegion(1:n1-1,:); nodes=n1-1;
[xMin,nXmin]=min(polyRegion(:,1));
[xMax,nXmax]=max(polyRegion(:,1));
%Start with minimum x node
polyRegion=[polyRegion(nXmin:nodes,:); polyRegion(1:nXmin-1,:)];
%Make cyclic
polyRegion=[polyRegion; polyRegion(1,:)];
if isempty(iR) | iR==1
    %Upper bounding segments
    for n=1:nXmax-1
        hold on
        if polyRegion(n+1,1)~=polyRegion(n,1)
            xx=linspace(polyRegion(n,1),polyRegion(n+1,1),100);
            plot(xx,line12(xx,polyRegion(n,:),polyRegion(n+1,:)),Col)
        else
            yy=linspace(polyRegion(n,2),polyRegion(n+1,2),100);
            plot(line21(yy,polyRegion(n,:),polyRegion(n+1,:)),yy,Col)
        end
    end
end
if isempty(iR) | iR==2
    %Lower bounding segments
    for n=nXmax:nodes
        hold on
        if polyRegion(n+1,1)~=polyRegion(n,1)
            xx=linspace(polyRegion(n,1),polyRegion(n+1,1),100);
            plot(xx,line12(xx,polyRegion(n,:),polyRegion(n+1,:)),Col)
        else
            yy=linspace(polyRegion(n,2),polyRegion(n+1,2),100);
            plot(line21(yy,polyRegion(n,:),polyRegion(n+1,:)),yy,Col)
        end
    end
end
return