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

kk1 = bsxfun(@lt,data,curbvals);
kk2 = bsxfun(@le,data,curbvals);

tt1 = sum(kk1,1)/nex*100;
tt2 = sum(kk2,1)/nex*100;
tt1 = squeeze(tt1);
tt2 = squeeze(tt2);
mm = abs(tt1-tt2);
dd = abs(bsxfun(@minus,tt1,selValues));

dd(mm>0.2) = 0;

if max(dd(:))<(100/numBins),
  updateBins = false;
end
