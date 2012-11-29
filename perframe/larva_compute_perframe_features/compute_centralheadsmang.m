%compute orientation of centralheadsm vector

function [data,units]=compute_centralheadsmang(trx,n)



larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
centralheadsmang=cell(1,numlarvae);

for i=1:numlarvae
    larva=larvae(i);
    centralheadsmang{1,i}=bsxfun(@atan2,trx(larva).yheadsm_mm-trx(larva).ycentral_mm,trx(larva).xheadsm_mm-trx(larva).xcentral_mm);
end

units=parseunits('rad');
data=centralheadsmang;
