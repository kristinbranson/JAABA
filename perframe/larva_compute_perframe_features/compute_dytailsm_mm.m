%compute dytailsm_mm

function [data,units]=compute_dytailsm_mm(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
dytailsm_mm=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    dytailsm_mm{1,i}=(trx(larva).ytailsm_mm(2:end)-trx(larva).ytailsm_mm(1:end-1))./trx(larva).dt;
end

units=parseunits('mm/s');
data=dytailsm_mm;