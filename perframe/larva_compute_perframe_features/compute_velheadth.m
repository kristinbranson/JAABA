%compute head speed in the tail-head direction
function [data,units]=compute_velheadth(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
velheadth=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    velheadth{1,i}=trx(larva).velmaghead.*(cos(trx(larva).velanghead).*cos(trx(larva).tailheadang(1,1:end-1))+sin(trx(larva).velanghead).*sin(trx(larva).tailheadang(1,1:end-1)));
end

units=parseunits('mm/s');
data=velheadth;