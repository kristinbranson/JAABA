%compute dyheadsm_mm

function [data,units]=compute_dyheadsm_mm(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
dyheadsm_mm=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    dyheadsm_mm{1,i}=(trx(larva).yheadsm_mm(2:end)-trx(larva).yheadsm_mm(1:end-1))./trx(larva).dt;
end

units=parseunits('mm/s');
data=dyheadsm_mm;