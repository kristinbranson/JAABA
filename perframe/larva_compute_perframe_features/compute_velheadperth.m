%compute head speed in the direction perpendicular to tail-head direction
function [data,units]=compute_velheadperth(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
velheadperth=cell(1,numlarvae);
tailheadangperp=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    tailheadangperp{1,i}=trx(larva).tailheadang-pi/2;
    velheadperth{1,i}=trx(larva).velmaghead.*(cos(trx(larva).velanghead).*cos(trx(larva).tailheadang(1,1:end-1)+pi/2)+sin(trx(larva).velanghead).*sin(trx(larva).tailheadang(1,1:end-1)+pi/2));
end

units=parseunits('mm/s');
data=velheadperth;