function bmodel = fastBag(data,labels,binVals,bins,params,obj)

if nargin < 6
  obj = [];
end

numIters = 100;
wkSamples = 50;
numPerm = 300;

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
bcellmodel = {};
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
    curError = sum( ((tt.*modLabels)<0).*wt);
    if(curError>0.5);
      curError = 1-curError;
      wkRule.dir(ndx) = -wkRule.dir(ndx);
    end
  
    wkRule.error(ndx) = curError;
    wkRule.alpha(ndx) = 1-2*wkRule.error(ndx);
    extraScores = extraScores + tt*wkRule.alpha(ndx);
  
  end
  
  scores = scores+extraScores/wkSamples;

  bcellmodel{itt} = wkRule;
  
  tt = scores.*modLabels;
  wt = initWt./(1+exp(tt));
  wt = wt./sum(wt);

  if mod(itt,10)==0,
    etime = toc(clk);
    perct = itt/numIters*100;
    if ~isempty(obj),
      obj.SetStatus('Bagging: %d%% done. Time Remaining:%ds\n',...
        round(perct),round(etime*(100-perct)/perct));
    else
      fprintf('Bagging: %d%% done. Time Remaining:%d\n',...
        round(perct),round(etime*(100-perct)/perct));
    end

  end
  
end

bmodel = []; count = 1;
for ndx = 1:numel(bcellmodel)
  for bndx = 1:numel(bcellmodel{ndx}.dim)
    bmodel(count).dim = bcellmodel{ndx}.dim(bndx);
    bmodel(count).error = bcellmodel{ndx}.error(bndx);
    bmodel(count).dir = bcellmodel{ndx}.dir(bndx);
    bmodel(count).tr = bcellmodel{ndx}.tr(bndx);
    bmodel(count).alpha = bcellmodel{ndx}.alpha(bndx);
    count = count+1;
  end
end

%{
newalpha = zeros(1,numel(bmodel));
for ndx = 1:numPerm
  curidx = randperm(numel(bmodel));
  scores = zeros(numEx,1);
  curalpha = zeros(1,numel(bmodel));
  wt = initWt;
  for bndx = curidx(:)'

    tr = bmodel(bndx).tr;
    dir = bmodel(bndx).dir;
    dim = bmodel(bndx).dim;
    if dir>0,
      tt = ((data(:,dim)> tr)-0.5)*2;
    else
      tt = ((data(:,dim)<= tr)-0.5)*2;
    end
    curError = sum( (tt.*modLabels).*wt);
    if(curError<0);
      curError = 0;
    end
  
    curError = 0.5-curError/2;
    curalpha(bndx) = 1-2*curError;
    scores = scores + tt*curalpha(bndx);
    tt = scores.*modLabels;
    wt = initWt./(1+exp(tt));
    wt = wt./sum(wt);
    
  end
  newalpha = newalpha + curalpha;
  
end
newalpha = newalpha/numPerm;
for ndx = 1:numel(bmodel),
  bmodel(ndx).oldalpha = bmodel(ndx).alpha;
  bmodel(ndx).alpha = newalpha(ndx);
end
%}

end


function best = findWeakRuleBagged(data,labels,dist,binVals,bins,params,wkSamples)


numDim = size(data,2);

best.dim = 1;
best.error = 0.5;
best.dir = 1;
best.tr = 0;

% curBestErr = 0.5*ones(1,numDim);
% binNo = ones(1,numDim);
% bestDir = ones(1,numDim);

numS = params.numSample;

% always helpful to normalize
dist = dist / sum(dist);


% KB: copied from randsample
edges = min([0 cumsum(dist')],1); % protect against accumulated round-off
edges(end) = 1; % get the upper edge exact
[~,curSel] = histc(rand(numS,1),edges);

if any(curSel==0),
  warning('Sanity check: some samples are not being chosen');
  curSel(curSel==0) = size(data,1);
end

curBins = bins(:,curSel);
curLabels = labels(curSel);


% KB: precompute these
numBins = size(binVals,1)+1;
% indices with positive labels
idxpos = curLabels > 0;

% always normalize histograms by Z
Z = size(curBins,2);
Zpos = nnz(idxpos);
Zneg = Z - Zpos;

fracpos = Zpos/Z;
fracneg = Zneg/Z;

posCount = accummatrix(curBins(:,idxpos)',ones(Zpos,1),numBins)'/Z;
negCount = accummatrix(curBins(:,~idxpos)',ones(Zneg,1),numBins)'/Z;


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

% [minError minDim] = min(curBestErr); %#ok<UDIM>
% best.error = minError;   best.dim = minDim; 
% best.dir = bestDir(minDim);   best.tr = binVals(binNo(minDim),minDim);

[sortError sortDim] = sort(curBestErr , 'ascend');
best.error = sortError(1:wkSamples);   best.dim = sortDim(1:wkSamples);
best.dir = bestDir(best.dim);   
for ndx = 1:wkSamples
  best.tr(ndx) = binVals(binNo(best.dim(ndx)),best.dim(ndx));
end

end


