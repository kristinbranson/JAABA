%compute orientation of tailsmcentralvector

function [data,units]=compute_tailsmcentralang(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
tailsmcentralang=cell(1,numlarvae);

for i=1:numlarvae
    larva=larvae(i);
    tailsmcentralang{1,i}=bsxfun(@atan2,trx(larva).ycentral_mm-trx(larva).ytailsm_mm,trx(larva).xcentral_mm-trx(larva).xtailsm_mm);
end
units=parseunits('rad');
data=tailsmcentralang;
