%compute head speed in the central-head direction
function [data,units]=compute_velheadch(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
velheadch=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    velheadch{1,i}=trx(larva).velmaghead.*(cos(trx(larva).velanghead).*cos(trx(larva).centralheadang(1,1:end-1))+sin(trx(larva).velanghead).*sin(trx(larva).centralheadang(1,1:end-1)));
end

units=parseunits('mm/s');
data=velheadch;
