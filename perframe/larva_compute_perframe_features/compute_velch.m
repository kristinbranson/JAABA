%compute speed in the central-head direction
function [data,units]=compute_velch(trx,n)


larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
velch=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    % just a little faster
    velch{i} = trx(larva).vx_cm.*cos(trx(larva).centralheadang(1,1:end-1)) + ...
      trx(larva).vy_cm.*sin(trx(larva).centralheadang(1,1:end-1));
    % velch{1,i}=trx(larva).velmag_ctr.*(cos(trx(larva).velang).*cos(trx(larva).centralheadang(1,1:end-1))+sin(trx(larva).velang).*sin(trx(larva).centralheadang(1,1:end-1)));
end

units=parseunits('mm/s');
data=velch;

