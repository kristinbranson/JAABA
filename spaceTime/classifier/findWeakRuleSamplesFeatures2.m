function best = findWeakRuleSamplesFeatures2(data,labels,dist,binVals,bins)

numDim = size(data,2);

best.dim = 1;
best.error = 0.5;
best.dir = 1;
best.tr = 0;

curBestErr = 0.5*ones(1,numDim);
binNo = ones(1,numDim);
bestDir = ones(1,numDim);

numS = 2500;
curSel = zeros(1,numS);
cumD = cumsum(dist);
for ndx = 1:numS
  curSel(ndx) = sum(cumD<rand(1,1))+1;
end

curBins = bins(curSel,:);
curLabels = labels(curSel);
curSelDim = find(rand(1,numDim)>0.9);

parfor dimNdx = 1:numel(curSelDim)
  dim = curSelDim(dimNdx);
  numBins = size(binVals,1)+1;
  binNdx =  curBins(:,dim)+numBins*(curLabels>0);
  allCount = histc(binNdx,.5:(2*numBins+0.5));
  allCount = allCount/sum(allCount);
  posCount = allCount(numBins+1:end-1);
  negCount = allCount(1:numBins);
    
  posLeft = 0;
  posRight = sum(posCount);
  negLeft = 0;
  negRight = sum(negCount);
  
  err = zeros(1,size(binVals,1));
  dir = ones(1,size(binVals,1));
  for ndx = 1:size(binVals,1)
    posLeft = posLeft + posCount(ndx);
    posRight = posRight - posCount(ndx);
    negLeft = negLeft + negCount(ndx);
    negRight = negRight - negCount(ndx);

    err(ndx) = posRight+negLeft-posLeft-negRight;
    if(err(ndx)<0)
      err(ndx) = - err(ndx);
      dir(ndx) = -1;
    end
  end
  err = 0.5-err/2;
  [curBestErr(dimNdx) binNo(dimNdx)]= min(err);
  bestDir(dimNdx) = dir(binNo(dimNdx));
end

[minError minDimNdx] = min(curBestErr);
minDim = curSelDim(minDimNdx);
best.error = minError;   best.dim = minDim; 
best.dir = bestDir(minDimNdx);   best.tr = binVals(binNo(minDim),minDim);

end

%{
function [predError dir] = getError(data,label,tr, dist)

dir = 1;
predLabel = 2*(data>tr)-1;
predError = sum(( (predLabel.*label) ~=1).*dist);

if(predError>0.5); predError = 1-predError; dir = -1; end
end
%}
