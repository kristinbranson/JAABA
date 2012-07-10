% change in min distance from any point on the ellipse of one fly to the
% nose of another fly
function [data,units] = compute_ddell2nose(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  if trx(fly).nframes <= 1,
    data{i} = [];
  else
    data{i} = diff(trx(fly).dell2nose,1,2) ./ trx(fly).dt;
  end
end
units = parseunits('mm/s');