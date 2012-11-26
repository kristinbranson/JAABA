%compute velang
function [data,units]=compute_velang(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
velang=cell(1,numlarvae);
for i=1:size(velang,2);
    larva=larvae(i);
    velang{1,i}=bsxfun(@atan2,trx(larva).vy_cm,trx(larva).vx_cm);
end

units=parseunits('rad');
data=velang;
