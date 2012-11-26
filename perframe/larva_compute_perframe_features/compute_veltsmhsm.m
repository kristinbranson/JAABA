%compute speed in tsmhsme tailsm-headsm direction
function [data,units]=compute_veltsmhsm(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
veltsmhsm=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    veltsmhsm{1,i}=trx(larva).velmag.*(cos(trx(larva).velang).*cos(trx(larva).tailsmheadsmang(1,1:end-1))+sin(trx(larva).velang).*sin(trx(larva).tailsmheadsmang(1,1:end-1)));
end

units=parseunits('mm/s');
data=veltsmhsm;
