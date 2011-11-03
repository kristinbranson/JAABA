% projection of the center of rotation on the minor axis
function [data,units] = compute_corfrac_min(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  
  if trx(fly).nframes < 2,
    corfrac = zeros(2,0);
  else
    [corfrac,~] = center_of_rotation2(trx,fly,false);
  end
  
  data{i} = corfrac(2,:);

end
units = parseunits('unit');

