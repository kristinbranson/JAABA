%compute speed in the central-headsm direction
function [data,units]=compute_velchsm(trx,n)


larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
velchsm=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    velchsm{1,i}=trx(larva).velmag_ctr.*(cos(trx(larva).velang).*cos(trx(larva).centralheadsmang(1,1:end-1))+sin(trx(larva).velang).*sin(trx(larva).centralheadsmang(1,1:end-1)));
end

units=parseunits('mm/s');
data=velchsm;
