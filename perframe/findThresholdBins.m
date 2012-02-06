function bins = findThresholdBins(data,binVals)

numDim = size(data,2);
bins = zeros(size(data,2),size(data,1));
parfor dim = 1:numDim
  curD = data(:,dim);
  curVals = binVals(:,dim);
  curBins = sum(bsxfun(@gt,curD',curVals))+1;
  bins(dim,:) = curBins;
end

bins = uint8(bins);