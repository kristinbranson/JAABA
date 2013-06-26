function [data,units] = compute_area_mm(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  
  if trx(fly).nframes < 1,
    data{i} = zeros(1,0);
  else
    data{i} = trx(fly).area / trx.pxpermm^2;
  end

end
units = parseunits('mm^2');


