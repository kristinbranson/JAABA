function [binVals bins] = findThresholds(data,params)

numDim = size(data,2);
numBins = params.numBins;
temp = linspace(0,100,numBins+2);
prcValues = temp(2:end-1);
binVals = zeros(length(prcValues),numDim);
bins = zeros(size(data,2),size(data,1));
numPts = size(data,1);
if numPts>4000
  sampleSize = 4000;
else
  sampleSize = numPts;
end

parfor dim = 1:numDim
  curD = data(:,dim);
  rrand = randperm(numPts);
  sel = curD(rrand(1:sampleSize));
  curVals = prctile(sel,prcValues);
  curBins = sum(bsxfun(@gt,curD',curVals'))+1;
  binVals(:,dim) = curVals;
  bins(dim,:) = curBins;
end
