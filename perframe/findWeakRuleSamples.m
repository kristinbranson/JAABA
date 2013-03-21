function [best,curSel] = findWeakRuleSamples(data,labels,dist,binVals,bins,params)

numDim = size(data,2);

best.dim = 1;
best.error = 0.5;
best.dir = 1;
best.tr = 0;

% curBestErr = 0.5*ones(1,numDim);
% binNo = ones(1,numDim);
% bestDir = ones(1,numDim);

numS = params.numSample;

% KB: only sample if it is helpful
dosample = numS < numel(dist);

% always helpful to normalize
dist = dist / sum(dist);


% fot testing
% for dosample = [false,true],
%   
%   niters = 100;
%   
%   tic;
%   for iter = 1:niters,


if dosample,
  
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

else
  
  curBins = bins;
  curLabels = labels;
  curSel = true(size(labels));
end

% KB: precompute these
numBins = size(binVals,1)+1;
% indices with positive labels
idxpos = curLabels > 0;

% always normalize histograms by Z
if dosample,
  Z = size(curBins,2);
  Zpos = nnz(idxpos);
else
  Z = sum(dist);
  Zpos = sum(dist(idxpos));
end
Zneg = Z - Zpos;

fracpos = Zpos/Z;
fracneg = Zneg/Z;

if dosample,
  
  posCount = accummatrix(curBins(:,idxpos)',ones(Zpos,1),numBins)'/Z;
  negCount = accummatrix(curBins(:,~idxpos)',ones(Zneg,1),numBins)'/Z;
  
%   posCount = histc(curBins(:,idxpos),edges,2);
%   posCount = posCount(:,1:end-1) / Z;
%   negCount = histc(curBins(:,~idxpos),edges,2);
%   negCount = negCount(:,1:end-1) / Z;

else

  % weighted histogram, loop-free version

  posCount = accummatrix(curBins(:,idxpos)',dist(idxpos),numBins)';
  negCount = accummatrix(curBins(:,~idxpos)',dist(~idxpos),numBins)';
  
  
end

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

[minError minDim] = min(curBestErr); %#ok<UDIM>
best.error = minError;   best.dim = minDim; 
best.dir = bestDir(minDim);   best.tr = binVals(binNo(minDim),minDim);

%   end
%   
%   fprintf('dosample = %d, time = %f\n',dosample,toc);
%   
% end

end

%{
function [predError dir] = getError(data,label,tr, dist)

dir = 1;
predLabel = 2*(data>tr)-1;
predError = sum(( (predLabel.*label) ~=1).*dist);

if(predError>0.5); predError = 1-predError; dir = -1; end
end
%}
