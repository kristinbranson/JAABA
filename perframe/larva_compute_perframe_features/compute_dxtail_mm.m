%compute dxtail_mm

function [data,units]=compute_dxtail_mm(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
dxtail_mm=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    dxtail_mm{1,i}=(trx(larva).xtail_mm(2:end)-trx(larva).xtail_mm(1:end-1))./trx(larva).dt;
end

units=parseunits('mm/s');
data=dxtail_mm;