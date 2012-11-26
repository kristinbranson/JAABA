%compute orientation of centralhead vector

function [data,units]=compute_centralheadang(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
centralheadang=cell(1,numlarvae);

for i=1:numlarvae
    larva=larvae(i);
    centralheadang{1,i}=bsxfun(@atan2,trx(larva).yhead_mm-trx(larva).ycentral_mm,trx(larva).xhead_mm-trx(larva).xcentral_mm);
end

units=parseunits('rad');
data=centralheadang;