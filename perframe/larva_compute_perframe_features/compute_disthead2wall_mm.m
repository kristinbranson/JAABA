% distance from the wall to the head
function [data,units] = compute_disthead2wall_mm(trx,n)
larvae = trx.exp2flies{n};
numlarvae = numel(larvae);
data = cell(1,numlarvae);
for i = 1:numlarvae,
  larva = larvae(i);  
  data{i} = 30-trx(larva).rhead_mm;
end
units = parseunits('mm');