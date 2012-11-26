%compute head speed in the tail-central direction
function [data,units]=compute_velheadtc(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
velheadtc=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    velheadtc{1,i}=trx(larva).velmaghead.*(cos(trx(larva).velanghead).*cos(trx(larva).tailcentralang(1,1:end-1))+sin(trx(larva).velanghead).*sin(trx(larva).tailcentralang(1,1:end-1)));
end

units=parseunits('mm/s');
data=velheadtc;