function [binVals] = findThresholds(data,params)

numDim = size(data,2);
numBins = params.numBins;
temp = linspace(0,100,numBins+2);
prcValues = temp(2:end-1);
binVals = zeros(length(prcValues),numDim);
% binVals3 = zeros(1,numDim,length(prcValues));
% bins = zeros(size(data,2),size(data,1));
numPts = size(data,1);
if numPts>4000
  sampleSize = 4000;
else
  sampleSize = numPts;
end

% tic
% 
% rrand = reshape(randsample(numPts,numDim*sampleSize,true),[sampleSize,numDim]);
% sel = sub2ind([]
% 
% toc

parfor dim = 1:numDim
  curD = data(:,dim);
%   rrand = randperm(numPts);
%   sel = curD(rrand(1:sampleSize));
  rrand = randperm(numPts,sampleSize);
  sel = curD(rrand);
  curVals = prctile(sel,prcValues);
  binVals(:,dim) = curVals;
%   binVals3(1,dim,:) = curVals;
%   curBins = sum(bsxfun(@gt,curD,curVals),2)+1;
%   bins(dim,:) = curBins;
end

% toc



% bins = uint8(ones(size(data,2),size(data,1)));
% 
% blockSize = 50000;
% numBlocks = ceil(size(data,1)/blockSize);
% 
% for ndx = 1:numBlocks
%    curb = (blockSize*(ndx-1)+1):(min(blockSize*ndx,size(data,1)));
%    bins(:,curb) = sum(bsxfun(@gt,data(curb,:),binVals3),3)'+1;
% end

% tbins = sum(bsxfun(@gt,data,binVals3),3)+1;
% bins = uint8(tbins');
