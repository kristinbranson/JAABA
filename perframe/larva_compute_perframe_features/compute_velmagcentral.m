%compute velmag_central
function [data,units]=compute_velmagcentral(trx,n)
larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
velmagcentral=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    velmagcentral{1,i}=bsxfun(@hypot,trx(larva).dxcentral_mm,trx(larva).dycentral_mm);
end

units=parseunits('mm/s');
data=velmagcentral;