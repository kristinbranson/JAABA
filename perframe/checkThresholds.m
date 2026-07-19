function updateBins = checkThresholds(data,params,binVals)
% function updateBins = checkThresholds(data,params,binVals)
% Checks if the binVals are still accurate.

updateBins = true;
if isempty(binVals)|| size(binVals,1)~=params.numBins,
  return;
end

nex = size(data,1);
numDim = size(data,2);
numBins = params.numBins;
temp = linspace(0,100,numBins+2);
prcValues = temp(2:end-1);
mpoint = round(numel(prcValues)/2);
% selpts = [1 mpoint numel(prcValues)];
selpts = [1 numel(prcValues)];
selValues = prcValues(selpts);
binVals = permute(binVals,[3 2 1]);
curbvals = binVals(:,:,selpts);

% The comparisons are summed over examples straight away, so they are
% accumulated in blocks of dimensions. Comparing every dimension at once
% would materialize an nex x numDim x numel(selpts) array, which is larger
% than the training data itself.

% OLD code 
% kk1 = bsxfun(@lt,data,curbvals);
% kk2 = bsxfun(@le,data,curbvals);
%
% tt1 = sum(kk1,1)/nex*100;
% tt2 = sum(kk2,1)/nex*100;

tt1 = zeros(1,numDim,numel(selpts));
tt2 = zeros(1,numDim,numel(selpts));
blocksize = max(1,floor(1e8/max(nex*numel(selpts),1)));
for j0 = 1:blocksize:numDim,
  j1 = min(j0+blocksize-1,numDim);
  curdata = data(:,j0:j1);
  curb = curbvals(1,j0:j1,:);
  tt1(1,j0:j1,:) = sum(bsxfun(@lt,curdata,curb),1);
  tt2(1,j0:j1,:) = sum(bsxfun(@le,curdata,curb),1);
end

tt1 = tt1/nex*100;
tt2 = tt2/nex*100;
tt1 = squeeze(tt1);
tt2 = squeeze(tt2);
mm = abs(tt1-tt2);
dd = abs(bsxfun(@minus,tt1,selValues));

dd(mm>0.2) = 0;

if max(dd(:))<(100/numBins),
  updateBins = false;
end
