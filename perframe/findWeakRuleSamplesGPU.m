function best = findWeakRuleSamplesGPU(data,labels,dist,binVals,bins,params)

numDim = size(data,2);

best.dim = 1;
best.error = 0.5;
best.dir = 1;
best.tr = 0;

curBestErr = 0.5*ones(1,numDim);
% BJA: vestigial code?
%binNo = ones(1,numDim);
%bestDir = ones(1,numDim);

numS = params.numSample;
curSel = zeros(1,numS);
cumD = cumsum(dist);
randPts = sort(rand(numS,1),'ascend');
count = 1; ndx = 1;
numel_cumD = numel(cumD);
% BJA: this while loop is slow
while count<=numS
  if ndx > numel_cumD
    curSel(count) = numel_cumD;
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
edges = gpuArray(edges);
% indices with positive labels
idxpos = curLabels > 0;
% always normalize histograms by Z
% BJA:  why though when it cancels out?
%Z = gpuArray(size(curBins,2));

fracpos = gpuArray(nnz(idxpos));  %/Z;
fracneg = gpuArray(nnz(~idxpos));  %/Z;
idxpos=gpuArray(idxpos);
curBinsP=gpuArray(curBins(:,idxpos));
curBinsN=gpuArray(curBins(:,~idxpos));

%x=uint8(floor(32*rand(1716,1200)));
%e=(0.5:31.5);
%xg=gpuArray(x);
%eg=gpuArray(e);

%disp(size(curBins(:,idxpos)));
%disp(size(edges));

%tic;
%posLeft=bsxfun(@lt,shiftdim(repmat(curBins(:,idxpos),[1 1 length(edges)]),1),edges);
%posLeft=transpose(squeeze(sum(posLeft,1)));
%posLeft=posLeft(:,2:end)-repmat(posLeft(:,1),1,size(posLeft,2)-1);
posLeft=cumsum_histc_gpu(curBinsP,edges);
%posLeft=cumsum_histc_gpu(xg,eg);
%posLeft=posLeft/Z;
%toc

%tic;
%negLeft=bsxfun(@lt,shiftdim(repmat(curBins(:,~idxpos),[1 1 length(edges)]),1),edges);
%negLeft=transpose(squeeze(sum(negLeft,1)));
%negLeft=negLeft(:,2:end)-repmat(negLeft(:,1),1,size(negLeft,2)-1);
negLeft=cumsum_histc_gpu(curBinsN,edges);
%negLeft=cumsum_histc_gpu(xg,eg);
%negLeft=negLeft/Z;
%toc

%posRight = fracpos - posLeft;
%negRight = fracneg - negLeft;

%binErr = posRight+negLeft-posLeft-negRight;
binErr = fracpos-2*posLeft+2*negLeft-fracneg;
negErr = binErr<0;
binErr(negErr) = -binErr(negErr);
dir = parallel.gpu.GPUArray.ones(numDim,numBins);
dir(negErr) = -1;
err = 0.5-binErr/2;
[curBestErr,binNo]= min(err(:,1:end-1),[],2);
curBestErr = curBestErr'; binNo = binNo';
bestDir = dir(sub2ind([numDim,numBins],1:numDim,binNo));

[minError minDim] = min(curBestErr);
best.error = gather(minError);  % BJA:  off by Z now, but not used i think.
best.dim = gather(minDim); 
best.dir = gather(bestDir(minDim));
best.tr = gather(binVals(binNo(minDim),minDim));





function N=cumsum_histc_gpu(X,EDGES)

for(i=1:length(EDGES))
  N(:,i)=sum(X<EDGES(i),2);
end

% BJA: this alternate formulations is 2x slower
%N=bsxfun(@lt,repmat(X,[1 1 length(EDGES)]),reshape(EDGES,1,1,length(EDGES)));
%N=squeeze(sum(N,2));

N=N(:,2:end)-repmat(N(:,1),1,size(N,2)-1);
