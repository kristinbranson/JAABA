%compute speed in the direction perpendicular to central-head direction
function [data,units]=compute_velperch(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
velperch=cell(1,numlarvae);
%centralheadangperp=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    %centralheadangperp{1,i}=trx(larva).centralheadang-pi/2;
    velperch{1,i}=trx(larva).velmag_ctr.*(cos(trx(larva).velang).*cos(trx(larva).centralheadang(1,1:end-1)+pi/2)+sin(trx(larva).velang).*sin(trx(larva).centralheadang(1,1:end-1)+pi/2));
end

units=parseunits('mm/s');
data=velperch;
