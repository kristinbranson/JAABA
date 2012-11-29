%compute orientation of inflhead vector

function [data,units]=compute_inflheadang(trx,n)


larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
inflheadang=cell(1,numlarvae);

for i=1:numlarvae
    larva=larvae(i);
    inflheadang{1,i}=bsxfun(@atan2,trx(larva).yhead_mm-trx(larva).yinflection_mm,trx(larva).xhead_mm-trx(larva).xinflection_mm);
end

units=parseunits('rad');
data=inflheadang;
