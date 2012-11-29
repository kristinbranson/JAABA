%compute dxcentral_mm

function [data,units]=compute_dxcentral_mm(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
dxcentral_mm=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    dxcentral_mm{1,i}=(trx(larva).xcentral_mm(2:end)-trx(larva).xcentral_mm(1:end-1))./trx(larva).dt;
end

units=parseunits('mm/s');
data=dxcentral_mm;