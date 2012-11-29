%compute velang_central
function [data,units]=compute_velangcentral(trx,n)
larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
velangcentral=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    velangcentral{1,i}=bsxfun(@atan2,trx(larva).dycentral_mm,trx(larva).dxcentral_mm);
end

units=parseunits('rad');
data=velangcentral;