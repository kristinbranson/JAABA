% change in distance to closest fly center
function [data,units] = compute_ddcenter(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  if trx(fly).nframes <= 1,
    data{i} = [];
  else
    data{i} = diff(trx(fly).dcenter,1,2) ./ trx(fly).dt;
  end
end
units = parseunits('mm/s');