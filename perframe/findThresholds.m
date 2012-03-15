function [binVals bins] = findThresholds(data,params)

numDim = size(data,2);
numBins = params.numBins;
temp = linspace(0,100,numBins+2);
prcValues = temp(2:end-1);
binVals = zeros(length(prcValues),numDim);
binVals3 = zeros(1,numDim,length(prcValues));
% bins = zeros(size(data,2),size(data,1));
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
  binVals(:,dim) = curVals;
  binVals3(1,dim,:) = curVals;
%   curBins = sum(bsxfun(@gt,curD,curVals),2)+1;
%   bins(dim,:) = curBins;
end
tbins = sum(bsxfun(@gt,data,binVals3),3)+1;
bins = uint8(tbins');
