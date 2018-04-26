% AR 3/8/2017
function [data,units] = compute_dell2ell(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);

for i1 = 1:nflies
    fly1 = flies(i1);
    % access closestfly to ensure the dell2ell is computed
    trx(fly1).closestfly_ell2ell;
    data{i1} = trx(fly1).dell2ell;
end
units = parseunits('mm')
