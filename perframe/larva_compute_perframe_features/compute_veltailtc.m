%compute tail speed in the tail-central direction
function [data,units]=compute_veltailtc(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
veltailtc=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    veltailtc{1,i}=trx(larva).velmagtail.*(cos(trx(larva).velangtail).*cos(trx(larva).tailcentralang(1,1:end-1))+sin(trx(larva).velangtail).*sin(trx(larva).tailcentralang(1,1:end-1)));
end

units=parseunits('mm/s');
data=veltailtc;