%compute speed in the direction perpendicular to tail-head direction
function [data,units]=compute_velperth(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
velperth=cell(1,numlarvae);
%tailheadangperp=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    %tailheadangperp{1,i}=trx(larva).tailheadang-pi/2;
    velperth{1,i}=trx(larva).velmag_ctr.*(cos(trx(larva).velang).*cos(trx(larva).tailheadang(1,1:end-1)+pi/2)+sin(trx(larva).velang).*sin(trx(larva).tailheadang(1,1:end-1)+pi/2));
end

units=parseunits('mm/s');
data=velperth;
