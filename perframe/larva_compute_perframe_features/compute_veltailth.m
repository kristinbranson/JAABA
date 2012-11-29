%compute tail speed in the tail-head direction
function [data,units]=compute_veltailth(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
veltailth=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    veltailth{1,i}=trx(larva).velmagtail.*(cos(trx(larva).velangtail).*cos(trx(larva).tailheadang(1,1:end-1))+sin(trx(larva).velangtail).*sin(trx(larva).tailheadang(1,1:end-1)));
end

units=parseunits('mm/s');
data=veltailth;