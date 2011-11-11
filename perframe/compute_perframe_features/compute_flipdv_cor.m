% value of sideways motion of center of rotation with sign set by change in
% orientation sign
function [data,units] = compute_flipdv_cor(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  
  data{i} = trx(fly).dv_cor .* trx(fly).signdtheta;

end
units = parseunits('mm/s');

