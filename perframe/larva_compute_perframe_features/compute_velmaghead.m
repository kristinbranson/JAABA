%compute velmag_head
function [data,units]=compute_velmaghead(trx,n)
larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
velmaghead=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    velmaghead{1,i}=bsxfun(@hypot,trx(larva).dxhead_mm,trx(larva).dyhead_mm);
end

units=parseunits('mm/s');
data=velmaghead;
 
