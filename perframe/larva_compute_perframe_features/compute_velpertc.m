%compute speed in the direction perpendicular to tail-central direction
function [data,units]=compute_velpertc(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
velpertc=cell(1,numlarvae);
%tailcentralangperp=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    %tailcentralangperp{1,i}=trx(larva).tailcentralang-pi/2;
    velpertc{1,i}=trx(larva).velmag_ctr.*(cos(trx(larva).velang).*cos(trx(larva).tailcentralang(1,1:end-1)+pi/2)+sin(trx(larva).velang).*sin(trx(larva).tailcentralang(1,1:end-1)+pi/2));
end

units=parseunits('mm/s');
data=velpertc;
