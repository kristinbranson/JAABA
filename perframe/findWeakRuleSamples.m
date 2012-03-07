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
  if ndx > numel(cumD)
    curSel(count) = numel(cumD);
    count = count+1;
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

% KB: precompute these
numBins = size(binVals,1)+1;
% for histogramming
edges = .5:numBins+.5;
% indices with positive labels
idxpos = curLabels > 0;
% always normalize histograms by Z
Z = size(curBins,2);

DOLOOP = false;
if ~DOLOOP,
% do this without looping at all
fracpos = nnz(idxpos)/Z;
fracneg = nnz(~idxpos)/Z;
posCount = histc(curBins(:,idxpos),edges,2);
posCount = posCount(:,1:end-1) / Z;
negCount = histc(curBins(:,~idxpos),edges,2);
negCount = negCount(:,1:end-1) / Z;
posLeft = cumsum(posCount,2);
posRight = fracpos - posLeft;
negLeft = cumsum(negCount,2);
negRight = fracneg - negLeft;

binErr = posRight+negLeft-posLeft-negRight;
negErr = binErr<0;
binErr(negErr) = -binErr(negErr);
dir = ones(numDim,numBins);
dir(negErr) = -1;
err = 0.5-binErr/2;
[curBestErr,binNo]= min(err(:,1:end-1),[],2);
curBestErr = curBestErr'; binNo = binNo';
bestDir = dir(sub2ind([numDim,numBins],1:numDim,binNo));

else

% pre-split up curBins for parfor
curBinsPos = curBins(:,idxpos);
curBinsNeg = curBins(:,~idxpos);

parfor dim = 1:numDim

%{
%   binNdx =  curBins(dim,:)+uint8(numBins*(curLabels'>0));
%   allCount = histc(binNdx,.5:(2*numBins+0.5));
%   allCount = allCount/sum(allCount);
%   posCount = allCount(numBins+1:end-1);
%   negCount = allCount(1:numBins);
%}
  posCount = histc(curBinsPos(dim,:),edges);
  posCount = posCount(1:end-1) / Z;
  negCount = histc(curBinsNeg(dim,:),edges);
  negCount = negCount(1:end-1) / Z;
  
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

end % end DOLOOP condition

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
