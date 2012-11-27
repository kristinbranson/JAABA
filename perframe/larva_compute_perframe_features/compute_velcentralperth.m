%compute central speed in the direction perpendicular to tail-head direction
function [data,units]=compute_velcentralperth(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
velcentralperth=cell(1,numlarvae);
tailheadangperp=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    tailheadangperp{1,i}=trx(larva).tailheadang-pi/2;
    velcentralperth{i} = trx(larva).dxcentral_mm.*cos(trx(larva).tailheadang(1:end-1)+pi/2) + ...
      trx(larva).dycentral_mm.*sin(trx(larva).tailheadang(1:end-1)+pi/2);
    % tailheadangperp not defined
    %velcentralperth{1,i}=trx(larva).velmagcentral.*(cos(trx(larva).velangcentral).*cos(trx(larva).tailheadangperp(1,1:end-1))+sin(trx(larva).velangcentral).*sin(trx(larva).tailheadangperp(1,1:end-1)));
end

units=parseunits('mm/s');
data=velcentralperth;