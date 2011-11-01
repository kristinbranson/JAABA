function best = findWeakRule(data,labels,dist,binVals,bins)

numDim = size(data,2);

best.dim = 1;
best.error = 0.5;
best.dir = 1;
best.tr = 0;


curBestErr = 0.5*ones(1,numDim);
binNo = ones(1,numDim);
bestDir = ones(1,numDim);

parfor dim = 1:numDim
  
  numBins = size(binVals,1)+1;
  binNdx =  bins(:,dim)+numBins*(labels>0);
  allCount = zeros(1,numBins*2);
  for count=1:2*numBins; 
    allCount(count) = sum(dist(binNdx==count)); 
  end; 
  posCount = allCount(numBins+1:end);
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
  [curBestErr(dim) binNo(dim)]= min(err);
  bestDir(dim) = dir(binNo(dim));
end

[minError minDim] = min(curBestErr);
best.error = minError;   best.dim = minDim; 
best.dir = bestDir(minDim);   best.tr = binVals(binNo(minDim),minDim);

end

%{
function [predError dir] = getError(data,label,tr, dist)

dir = 1;
predLabel = 2*(data>tr)-1;
predError = sum(( (predLabel.*label) ~=1).*dist);

if(predError>0.5); predError = 1-predError; dir = -1; end
end
%}
