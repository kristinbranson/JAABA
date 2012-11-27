%compute central speed in the central-head direction
function [data,units]=compute_velcentralch(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
velcentralch=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    % this is just slightly faster
    velcentralch{i} = trx(larva).dxcentral_mm.*cos(trx(larva).centralheadang(1:end-1)) + ...
      trx(larva).dycentral_mm.*sin(trx(larva).centralheadang(1:end-1));
    %velcentralch1{1,i}=trx(larva).velmagcentral.*(cos(trx(larva).velangcentral).*cos(trx(larva).centralheadang(1,1:end-1))+sin(trx(larva).velangcentral).*sin(trx(larva).centralheadang(1,1:end-1)));
end

units=parseunits('mm/s');
data=velcentralch;