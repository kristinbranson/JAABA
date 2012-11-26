%compute dxinflection_mm

function [data,units]=compute_dxinflection_mm(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
dxinflection_mm=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    dxinflection_mm{1,i}=(trx(larva).xinflection_mm(2:end)-trx(larva).xinflection_mm(1:end-1))./trx(larva).dt;
end

units=parseunits('mm/s');
data=dxinflection_mm;