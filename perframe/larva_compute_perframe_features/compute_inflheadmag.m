%compute length of inflhead vector

function [data,units]=compute_inflheadmag(trx,n)


larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
inflheadmag=cell(1,numlarvae);

for i=1:numlarvae
    larva=larvae(i);
    inflheadmag{1,i}=bsxfun(@hypot,trx(larva).xhead_mm-trx(larva).xinflection_mm,trx(larva).yhead_mm-trx(larva).yinflection_mm);
end

units=parseunits('mm');
data=inflheadmag;
