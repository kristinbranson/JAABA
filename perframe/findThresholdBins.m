function bins = findThresholdBins(data,binVals)

numDim = size(data,2);
numBins = size(binVals,1);
if numBins+1 > intmax('uint8'),
  error('findThresholdBins:tooManyBins',...
    'bins are stored as uint8, which cannot represent %d bins.',numBins+1);
end

bins = zeros(numDim,size(data,1),'uint8');
parfor dim = 1:numDim
  curD = data(:,dim);
  curVals = binVals(:,dim);
  bins(dim,:) = uint8(sum(bsxfun(@gt,curD',curVals))+1);
end
