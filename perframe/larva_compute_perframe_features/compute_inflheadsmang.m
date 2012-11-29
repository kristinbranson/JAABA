%compute orientation of inflheadsm vector

function [data,units]=compute_inflheadsmang(trx,n)


larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
inflheadsmang=cell(1,numlarvae);

for i=1:numlarvae
    larva=larvae(i);
    inflheadsmang{1,i}=bsxfun(@atan2,trx(larva).yheadsm_mm-trx(larva).yinflection_mm,trx(larva).xheadsm_mm-trx(larva).xinflection_mm);
end

units=parseunits('rad');
data=inflheadsmang;