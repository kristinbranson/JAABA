%compute length of tailheadvector

function [data,units]=compute_tailheadmag(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
tailheadmag=cell(1,numlarvae);

for i=1:numlarvae
    larva=larvae(i);
    tailheadmag{1,i}=bsxfun(@hypot,trx(larva).xhead_mm-trx(larva).xtail_mm,trx(larva).yhead_mm-trx(larva).ytail_mm);
end

units=parseunits('mm');
data=tailheadmag;
