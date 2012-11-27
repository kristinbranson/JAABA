%compute central peed in the direction perpendicular to tail-central direction
function [data,units]=compute_velcentralpertc(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
velcentralpertc=cell(1,numlarvae);
tailcentralangperp=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    tailcentralangperp{1,i}=trx(larva).tailcentralang-pi/2;
    velcentralpertc{i} = trx(larva).dxcentral_mm.*cos(trx(larva).tailcentralang(1:end-1)+pi/2) + ...
      trx(larva).dycentral_mm.*sin(trx(larva).tailcentralang(1:end-1)+pi/2);
    % tailcentralangperp not defined
    %velcentralpertc{1,i}=trx(larva).velmagcentral.*(cos(trx(larva).velangcentral).*cos(trx(larva).tailcentralangperp(1,1:end-1))+sin(trx(larva).velangcentral).*sin(trx(larva).tailcentralangperp(1,1:end-1)));
end

units=parseunits('mm/s');
data=velcentralpertc;
