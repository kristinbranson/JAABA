%compute central speed in the tail-head direction
function [data,units]=compute_velcentralth(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
velcentralth=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    % this is just slightly faster
    velcentralth{i} = trx(larva).dxcentral_mm.*cos(trx(larva).tailheadang(1:end-1)) + ...
      trx(larva).dycentral_mm.*sin(trx(larva).tailheadang(1:end-1));
    %velcentralth{1,i}=trx(larva).velmagcentral.*(cos(trx(larva).velangcentral).*cos(trx(larva).tailheadang(1,1:end-1))+sin(trx(larva).velangcentral).*sin(trx(larva).tailheadang(1,1:end-1)));
end

units=parseunits('mm/s');
data=velcentralth;