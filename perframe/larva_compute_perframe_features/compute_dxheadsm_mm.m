%compute dxheadsm_mm

function [data,units]=compute_dxheadsm_mm(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
dxheadsm_mm=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    dxheadsm_mm{1,i}=(trx(larva).xheadsm_mm(2:end)-trx(larva).xheadsm_mm(1:end-1))./trx(larva).dt;
end

units=parseunits('mm/s');
data=dxheadsm_mm;