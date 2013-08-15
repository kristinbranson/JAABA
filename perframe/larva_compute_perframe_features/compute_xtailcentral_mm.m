%compute the xcoordinate of the center in the tails coordinate system

function [data,units]=compute_xtailcentral_mm(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
xtailcentral_mm=cell(1,numlarvae);

for i=1:numlarvae
    larva=larvae(i);
    xtailcentral_mm{1,i}=trx(larva).xcentral_mm-trx(larva).xtail_mm;
end

units=parseunits('mm');
data=xtailcentral_mm;
