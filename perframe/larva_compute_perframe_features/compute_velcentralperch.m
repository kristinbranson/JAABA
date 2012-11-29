%compute central speed in the direction perpendicular to the central-head direction
function [data,units]=compute_velcentralperch(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
velcentralperch=cell(1,numlarvae);
centralheadangperp=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    centralheadangperp{1,i}=trx(larva).centralheadang-pi/2;
    velcentralperch{i} = trx(larva).dxcentral_mm.*cos(trx(larva).centralheadang(1:end-1)+pi/2) + ...
      trx(larva).dycentral_mm.*sin(trx(larva).centralheadang(1:end-1)+pi/2);
    % centralheadangperp not defined
    %velcentralperch{1,i}=trx(larva).velmagcentral.*(cos(trx(larva).velangcentral).*cos(trx(larva).centralheadangperp(1,1:end-1))+sin(trx(larva).velangcentral).*sin(trx(larva).centralheadangperp(1,1:end-1)));
end

units=parseunits('mm/s');
data=velcentralperch;