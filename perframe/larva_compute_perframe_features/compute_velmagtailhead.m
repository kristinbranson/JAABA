%compute velmag_tailhead
function [data,units]=compute_velmagtailhead(trx,n)
larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
velmagtailhead=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    velmagtailhead{1,i}=bsxfun(@hypot,trx(larva).dxtailhead_mm,trx(larva).dytailhead_mm);
end

units=parseunits('rad');
data=velmagtailhead;