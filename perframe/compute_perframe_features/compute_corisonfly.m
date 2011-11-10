% whether the center of rotation is contained within the ellipse of the fly
function [data,units] = compute_corisonfly(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  
  if trx(fly).nframes < 2,
    data{i} = [];
  else
    [~,data{i}] = center_of_rotation2(trx,fly,false);
  end

end
units = parseunits('unit');

