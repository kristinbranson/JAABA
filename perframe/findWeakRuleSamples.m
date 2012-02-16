function best = findWeakRuleSamples(data,labels,dist,binVals,bins,params)

numDim = size(data,2);

best.dim = 1;
best.error = 0.5;
best.dir = 1;
best.tr = 0;

curBestErr = 0.5*ones(1,numDim);
binNo = ones(1,numDim);
bestDir = ones(1,numDim);

numS = params.numSample;
curSel = zeros(1,numS);
cumD = cumsum(dist);
randPts = sort(rand(numS,1),'ascend');
count = 1; ndx = 1;
while count<=numS
  if count > numel(cumD)
    curSel(count) = numel(cumD);
  else
    if(cumD(ndx)>=randPts(count))
      curSel(count) = ndx;
      count = count+1;
    else
      ndx = ndx+1;
    end
  end
end

curSel(curSel==0) = size(data,1);

curBins = bins(:,curSel);
curLabels = labels(curSel);
parfor dim = 1:numDim
  
  numBins = size(binVals,1)+1;
  binNdx =  curBins(dim,:)+uint8(numBins*(curLabels'>0));
  allCount = histc(binNdx,.5:(2*numBins+0.5));
  allCount = allCount/sum(allCount);
  posCount = allCount(numBins+1:end-1);
  negCount = allCount(1:numBins);
    
%{
%   posLeft = 0;
%   posRight = sum(posCount);
%   negLeft = 0;
%   negRight = sum(negCount);
%   
%   err = zeros(1,size(binVals,1));
%   dir = ones(1,size(binVals,1));
%   for ndx = 1:size(binVals,1)
%     posLeft = posLeft + posCount(ndx);
%     posRight = posRight - posCount(ndx);
%     negLeft = negLeft + negCount(ndx);
%     negRight = negRight - negCount(ndx);
% 
%     err(ndx) = posRight+negLeft-posLeft-negRight;
%     if(err(ndx)<0)
%       err(ndx) = - err(ndx);
%       dir(ndx) = -1;
%     end
%   end
%}

  dir = ones(1,numBins);
  posLeft = cumsum(posCount);
  posRight = sum(posCount)-posLeft;
  negLeft = cumsum(negCount);
  negRight = sum(negCount)-negLeft;
  binErr = posRight+negLeft-posLeft-negRight;
  negErr = binErr<0;
  binErr(negErr) = -binErr(negErr);
  dir(negErr) = -1;

  err = 0.5-binErr/2;
  [curBestErr(dim) binNo(dim)]= min(err(1:end-1));
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
