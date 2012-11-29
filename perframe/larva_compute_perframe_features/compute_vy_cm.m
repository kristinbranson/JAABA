%compute vy_cm
function [data,units]=compute_vy_cm(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
vy_cm=cell(1,numlarvae);

for i=1:numlarvae
    larva=larvae(i);
    vy_cm{1,i}=(trx(larva).y_mm(2:end)-trx(larva).y_mm(1:end-1))./trx(larva).dt;
end

units=parseunits('mm/s');
data=vy_cm;

