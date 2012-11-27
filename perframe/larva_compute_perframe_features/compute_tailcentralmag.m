%compute length of tailcentralvector

function [data,units]=compute_tailcentralmag(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
tailcentralmag=cell(1,numlarvae);

for i=1:numlarvae
    larva=larvae(i);
    tailcentralmag{1,i}=bsxfun(@hypot,trx(larva).xcentral_mm-trx(larva).xtail_mm,trx(larva).ycentral_mm-trx(larva).ytail_mm);
end

units=parseunits('mm');
data=tailcentralmag;
