% change in distance to wall
function [data,units] = compute_ddist2corner_rect(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);  
  if trx(fly).nframes <= 1,
    data{i} = [];
  else
    data{i} = diff(trx(fly).dist2corner_rect,1,2)./trx(fly).dt;
  end
end
units = parseunits('mm/s');

