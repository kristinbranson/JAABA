%compute speed in the tailsm-central direction
function [data,units]=compute_veltsmc(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
veltsmc=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    veltsmc{1,i}=trx(larva).velmag.*(cos(trx(larva).velang).*cos(trx(larva).tailsmcentralang(1,1:end-1))+sin(trx(larva).velang).*sin(trx(larva).tailsmcentralang(1,1:end-1)));
end

units=parseunits('mm/s');
data=veltsmc;