%compute velmag_tailcentral
function [data,units]=compute_velmagtailcentral(trx,n)
larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
velmagtailcentral=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    velmagtailcentral{1,i}=bsxfun(@hypot,trx(larva).dxtailcentral_mm,trx(larva).dytailcentral_mm);
end

units=parseunits('mm/s');
data=velmagtailcentral;