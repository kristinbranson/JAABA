% template computing absolute values of attributes
function [data,units] = compute_abs__template(trx,n, pFeatureName)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);  
  data{i} = abs(trx(fly).(pFeatureName));
end
units = trx.units.(pFeatureName);
%units = parseunits(pUnit);

