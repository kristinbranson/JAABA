%compute length of centralhead vector

function [data,units]=compute_centralheadmag(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
centralheadmag=cell(1,numlarvae);

for i=1:numlarvae
    larva=larvae(i);
    centralheadmag{1,i}=bsxfun(@hypot,trx(larva).xhead_mm-trx(larva).xcentral_mm,trx(larva).yhead_mm-trx(larva).ycentral_mm);
end

units=parseunits('mm');
data=centralheadmag;
