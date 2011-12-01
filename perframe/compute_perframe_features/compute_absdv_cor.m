% absolute value of sideways motion of center of rotation
function [data,units] = compute_absdv_cor(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  
  data{i} = abs(trx(fly).dv_cor);

end
units = parseunits('mm/s');

