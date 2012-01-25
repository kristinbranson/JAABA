function [data,units] = compute_angleonclosestfly(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);

for i1 = 1:nflies,
  fly1 = flies(i1);
  % access closestfly to ensure that dcenter is computed
  trx(fly1).closestfly_nose2ell;
  if trx.PerFrameDataExists('angleonclosestfly',n) == 0,
    % still doesn't exist, then compute it
    trx.ComputePerFrameData('closestfly_nose2ell',n);
  end
  data{i1} = trx(fly1).angleonclosestfly;
end
units = parseunits('rad');