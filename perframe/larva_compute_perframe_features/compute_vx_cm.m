%compute vx_cm
function [data,units]=compute_vx_cm(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
vx_cm=cell(1,numlarvae);

for i=1:numlarvae
    larva=larvae(i);
    vx_cm{1,i}=(trx(larva).x_mm(2:end)-trx(larva).x_mm(1:end-1))./trx(larva).dt;
end

units=parseunits('mm/s');
data=vx_cm;

