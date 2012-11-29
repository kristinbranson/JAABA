%compute velang_head
function [data,units]=compute_velanghead(trx,n)
larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
velanghead=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    velanghead{1,i}=bsxfun(@atan2,trx(larva).dyhead_mm,trx(larva).dxhead_mm);
end

units=parseunits('rad');
data=velanghead;

