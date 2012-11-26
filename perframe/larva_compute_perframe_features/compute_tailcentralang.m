%compute orientation of tailcentralvector

function [data,units]=compute_tailcentralang(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
tailcentralang=cell(1,numlarvae);

for i=1:numlarvae
    larva=larvae(i);
    tailcentralang{1,i}=bsxfun(@atan2,trx(larva).ycentral_mm-trx(larva).ytail_mm,trx(larva).xcentral_mm-trx(larva).xtail_mm);
end
units=parseunits('rad');
data=tailcentralang;
