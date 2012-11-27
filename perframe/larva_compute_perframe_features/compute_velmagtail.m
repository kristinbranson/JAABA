%compute velmag_tail
function [data,units]=compute_velmagtail(trx,n)
larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
velmagtail=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    velmagtail{1,i}=bsxfun(@hypot,trx(larva).dxtail_mm,trx(larva).dytail_mm);
end

units=parseunits('mm/s');
data=velmagtail;