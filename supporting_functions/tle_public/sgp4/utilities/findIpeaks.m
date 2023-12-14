function peaks=findIpeaks(datas,thresh,iSubset,varargin)
%    Find 2D peak clusters above threshold with size truncation
%  USAGE:
%    peaks=findPeaks2(datas,thresh,iSubset,varargin)
%    peaks=findPeaks2(datas,thresh,iSubset,minHits,dXmax,dYmax) 
%  INPUTS:
%   data      = nyXnx  data array  
%   threshold = in data units
%   iSubset   = indices of points included []=> full set 
%   minHits   = minimum # of hits for cluster        (default 2) 
%   dXmax, dYmax = max extent of cluster coordinates (default 5,5)
%  OUTPUTS:
%   peaks   = structure array
%         .pos = index location in datas of brightest peak [ix,iy]
%         .Ibar= value of brightest peak
%         .Hits= # of peaks in cluster
%

if isempty(varargin)
    %Parameters for peak extraction
    dXmax=5; dYmax=5;  %Maximum size
    %NOTE: Too small => multiple hits on same peak
    %      Too large => joins adjacent peaks
    minHits=2;          %Minimum # of hits in peak
else
    minHits=varargin{1};
    dXmax  =varargin{2};
    dYmax  =varargin{3};
end

[n1 n2] = size(datas);
data=zeros(n1,n2);
data(iSubset)=datas(iSubset);

%Find peaks above threshold
iPeaks = find(data(:) > thresh(:));
Peaks = data(iPeaks);
npeaks = length(iPeaks);

%Order the Peaks and find their coordinate indices
[Peaks,iPeakS] = sort(Peaks,'descend');
iPeaks = iPeaks(iPeakS);
[iy,ix] = ind2sub([n1 n2],iPeaks);

peaks=struct;
nPeaks2test = npeaks;
iCount = 0;
while nPeaks2test > minHits
    %Find the neighbors of the current brightest peak (include the test peak in the set)
    ixtest = ix(1); iytest = iy(1);
    
    %Test X and Y distances separately
    dX = abs(ix - ixtest);
    iHitsX = find(dX < dXmax);
    dY = abs(iy - iytest);
    iHitsY = find(dY < dYmax);
    iHitsP = intersect(iHitsX, iHitsY);    %Locations in Peaks
    iHits = (ix(iHitsP)-1)*n1+iy(iHitsP);  %Locations in data
    
    %Save the current peak if it has enough neighbors
    if length(iHits)>minHits
        iCount = iCount + 1;
        peakVals = [1;1]*Peaks(iHitsP)';                  %peak Values    
        %peakLocs = [iy(iHitsP)'; ix(iHitsP)'];           %peak Locations in data
        peaks(iCount).pos = [iytest; ixtest];             %Brightest peak in current cluster 
        %peakVarTmp = (peakLocs-repmat(peaks(iCount).pos, [1 length(peakLocs)])).^2;
        %sumpeakvals = sum(peakVals(1,:));
        %peakVarTmp = peakVarTmp .* peakVals ./ sumpeakvals;
        %peaks(iCount).posvar   = max(sum(peakVarTmp,2), [1; 1]);  % ensure >= 1 cell variance
        peaks(iCount).Ibar     = data(iytest, ixtest);
        peaks(iCount).Hits     = length(peakVals);
    end;

    %Remove the current brightest peak and its neighbors
    iPeaks = setdiff(iPeaks,iHits);
    Peaks = data(iPeaks);

    %Reorder the remaining peaks
    [Peaks,iPeakS] = sort(Peaks,'descend');
    iPeaks = iPeaks(iPeakS);
    [iy,ix] = ind2sub([n1 n2],iPeaks);
    nPeaks2test = length(iPeaks);
   
end;
return