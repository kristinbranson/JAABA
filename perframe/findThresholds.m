function [binVals bins] = findThresholds(data)

numDim = size(data,2);
binVals = zeros(length(2:3:98),numDim);
bins = zeros(size(data));
numPts = size(data,1);
if numPts>1000
  sampleSize = 1000;
else
  sampleSize = numPts;
end

parfor dim = 1:numDim
  curD = data(:,dim);
  rrand = randperm(numPts);
  sel = curD(rrand(1:sampleSize));
  curVals = prctile(sel,2:3:98);
  
  curBins = sum(bsxfun(@gt,curD',curVals'))+1;
%   for ndx = 1:numPts
%     curBins(ndx) = sum(curD(ndx)>curVals)+1;
%   end
  
  binVals(:,dim) = curVals;
  bins(:,dim) = curBins;
end
