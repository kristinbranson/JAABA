function best = findWeakRule(data,labels,dist,iterNo)

persistent binVals bins;

numDim = size(data,2);

minD = min(data);
maxD = max(data);
best.dim = 1;
best.error = 0.5;
best.dir = 1;
best.tr = 0;

if(iterNo==1)
  binVals = zeros(length(2:3:98),numDim);
  bins = zeros(size(data));
  for dim = 1:numDim
    if(size(data,1)>1000)
      rrand = randperm(size(data,1));
      sel = data(rrand(1:1000),dim);
    else
      sel = data(:,dim);
    end
    binVals(:,dim) = prctile(sel,2:3:98);

    for ndx = 1:size(data,1)
     bins(ndx,dim) = sum(data(ndx,dim)>binVals(:,dim))+1;
    end
  end
  
end

curBestErr = 0.5*ones(1,numDim);
binNo = ones(1,numDim);
bestDir = ones(1,numDim);

for dim = 1:numDim
  
%   stepSize = (maxD(dim)-minD(dim))/40;
%   eps = stepSize/10;
%   for tr = (minD(dim)+eps):stepSize:(maxD(dim)-eps)
% %   for tr = trVals
%     [curErr curDir] = getError(data(:,dim),labels, tr, dist);
%     if(curErr<best.error);
%       best.error = curErr;   best.dim = dim; 
%       best.dir = curDir;     best.tr = tr;
%     end
%   end
 
 
  posCount = zeros(1,size(binVals,1)+1);
  negCount = zeros(1,size(binVals,1)+1);
  
%   for ndx = 1:size(binVals,1)+1
%     curBinNdx = bins(:,dim)==ndx;
%     posCount(ndx) = sum(dist (curBinNdx & (labels>0) ));
%     negCount(ndx) = sum(dist (curBinNdx & (labels<0) ));
%   end

  for ndx = 1:size(bins,1)
    curBinNdx = bins(ndx,dim);
    if(labels(ndx)>0)
      posCount(curBinNdx) = posCount(curBinNdx) + dist(ndx);
    else
      negCount(curBinNdx) = negCount(curBinNdx) + dist(ndx);
    end
  end
    
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

function [predError dir] = getError(data,label,tr, dist)

dir = 1;
predLabel = 2*(data>tr)-1;
predError = sum(( (predLabel.*label) ~=1).*dist);

if(predError>0.5); predError = 1-predError; dir = -1; end
end