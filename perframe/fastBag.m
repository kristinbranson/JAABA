function bmodel = fastBag(data,labels,binVals,bins,params)

numIters = 100;
wkSamples = 10;

numEx = size(data,1);
posEx = labels == 1;
negEx = ~posEx;
numPos = sum(posEx);
numNeg = sum(negEx);

if numPos<1 || numNeg<1,
  return;
end

modLabels = sign( (labels==1)-0.5);


initWt = getWeights(modLabels);

clk = tic;

scores = zeros(numEx,1);
wt = initWt;
bmodel = {};
for itt = 1:numIters
  wkRule = findWeakRuleBagged(data,modLabels,wt,binVals,bins,params,wkSamples);

  extraScores = zeros(numEx,1);
  for ndx = 1:wkSamples
    tr = wkRule.tr(ndx);
    dir = wkRule.dir(ndx);
    dim = wkRule.dim(ndx);
    if dir>0,
      tt = ((data(:,dim)> tr)-0.5)*2;
    else
      tt = ((data(:,dim)<= tr)-0.5)*2;
    end
    curError = sum( (tt.*labels).*wt);
    if(curError<0);
      curError = -curError;
      wkRule.dir(ndx) = -wkRule.dir(ndx);
    end
  
    wkRule.error(ndx) = 0.5-curError/2;
    wkRule.alpha(ndx) = 1-2*wkRule.error(ndx);
    extraScores = extraScores + tt*wkRule.alpha(ndx);
  
  end
  
  scores = scores+extraScores/wkSamples;

  bmodel{itt} = wkRule;
  
  tt = scores.*labels;
  wt = initWt./(1+exp(tt));
  wt = wt./sum(wt);

  if mod(itt,10)==0,
    etime = toc(clk);
    perct = itt/numIters*100;
    fprintf('%d%% done. Time Remaining:%.2f\n',...
      round(perct),etime*(100-perct)/perct);
  end
  
end
end


function best = findWeakRuleBagged(data,labels,dist,binVals,bins,params,wkSamples)

numDim = size(data,2);
numEx = size(data,1);

best.dim = ones(1,wkSamples);
best.error = 0.5*ones(1,wkSamples);
best.dir = ones(1,wkSamples);
best.tr = zeros(1,wkSamples);

curBestErr = 0.5*ones(1,numDim);
binNo = ones(1,numDim);
bestDir = ones(1,numDim);

numBins = size(binVals,1)+1;
% for histogramming
edges = .5:numBins+.5;
% indices with positive labels
% always normalize histograms by Z


numS = params.numSample;
curSel = randsample(numEx,numS,true,dist);

parfor dim = 1:numDim
  
  curBins = bins(dim,curSel);
  curLabels = labels(curSel);
  idxpos = curLabels > 0;
  % pre-split up curBins for parfor
  curBinsPos = curBins(idxpos);
  curBinsNeg = curBins(~idxpos);
  
  Zpos = nnz(idxpos); Zneg = nnz(~idxpos);
  posCount = histc(curBinsPos,edges);
  posCount = posCount(1:end-1) / Zpos;
  negCount = histc(curBinsNeg,edges);
  negCount = negCount(1:end-1) / Zneg;
  
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

[sortError sortDim] = sort(curBestErr,'ascend');
best.error = sortError(1:wkSamples);   best.dim = sortDim(1:wkSamples);
best.dir = bestDir(best.dim);   
for ndx = 1:wkSamples
  best.tr(ndx) = binVals(binNo(best.dim(ndx)),best.dim(ndx));
end

end

