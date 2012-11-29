%compute length of tailinflvector

function [data,units]=compute_tailinflmag(trx,n)


larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
tailinflmag=cell(1,numlarvae);

for i=1:numlarvae
    larva=larvae(i);
    tailinflmag{1,i}=bsxfun(@hypot,trx(larva).xinflection_mm-trx(larva).xtail_mm,trx(larva).yinflection_mm-trx(larva).ytail_mm);
end

units=parseunits('mm');
data=tailinflmag;
