%compute tail speed in the central-head direction
function [data,units]=compute_veltailch(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
veltailch=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    veltailch{1,i}=trx(larva).velmagtail.*(cos(trx(larva).velangtail).*cos(trx(larva).centralheadang(1,1:end-1))+sin(trx(larva).velangtail).*sin(trx(larva).centralheadang(1,1:end-1)));
end

units=parseunits('mm/s');
data=veltailch;