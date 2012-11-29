%compute dytail_mm

function [data,units]=compute_dytail_mm(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
dytail_mm=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    dytail_mm{1,i}=(trx(larva).ytail_mm(2:end)-trx(larva).ytail_mm(1:end-1))./trx(larva).dt;
end

units=parseunits('mm/s');
data=dytail_mm;