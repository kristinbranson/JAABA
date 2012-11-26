%compute central speed in the tail-central direction
function [data,units]=compute_velcentraltc(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
velcentraltc=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    % this is just slightly faster
    velcentraltc{i} = trx(larva).dxcentral_mm.*cos(trx(larva).tailcentralang(1:end-1)) + ...
      trx(larva).dycentral_mm.*sin(trx(larva).tailcentralang(1:end-1));
    %velcentraltc{1,i}=trx(larva).velmagcentral.*(cos(trx(larva).velangcentral).*cos(trx(larva).tailcentralang(1,1:end-1))+sin(trx(larva).velangcentral).*sin(trx(larva).tailcentralang(1,1:end-1)));
end

units=parseunits('mm/s');
data=velcentraltc;