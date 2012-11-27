%compute speed in the tail-central direction
function [data,units]=compute_veltc(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
veltc=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    veltc{1,i}=trx(larva).velmag_ctr.*(cos(trx(larva).velang).*cos(trx(larva).tailcentralang(1,1:end-1))+sin(trx(larva).velang).*sin(trx(larva).tailcentralang(1,1:end-1)));
end

units=parseunits('mm/s');
data=veltc;

