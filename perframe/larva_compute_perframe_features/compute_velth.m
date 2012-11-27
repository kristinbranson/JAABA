%compute speed in the tail-head direction
function [data,units]=compute_velth(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
velth=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    velth{1,i}=trx(larva).velmag_ctr.*(cos(trx(larva).velang).*cos(trx(larva).tailheadang(1,1:end-1))+sin(trx(larva).velang).*sin(trx(larva).tailheadang(1,1:end-1)));
end

units=parseunits('mm/s');
data=velth;

