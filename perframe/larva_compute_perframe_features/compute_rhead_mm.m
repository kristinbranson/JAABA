% distance from radius of the head to the center
function [data,units] = compute_rhead_mm(trx,n)
larvae = trx.exp2flies{n};
numlarvae = numel(larvae);
data = cell(1,numlarvae);
for i = 1:numlarvae,
  larva = larvae(i);  
  data{i} = bsxfun(@hypot,trx(larva).yhead_mm,trx(larva).xhead_mm);
end
units = parseunits('mm');