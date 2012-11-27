%compute dxtailsm_mm

function [data,units]=compute_dxtailsm_mm(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
dxtailsm_mm=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    dxtailsm_mm{1,i}=(trx(larva).xtailsm_mm(2:end)-trx(larva).xtailsm_mm(1:end-1))./trx(larva).dt;
end

units=parseunits('mm/s');
data=dxtailsm_mm;