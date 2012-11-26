%compute dyinflection_mm

function [data,units]=compute_dyinflection_mm(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
dyinflection_mm=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    dyinflection_mm{1,i}=(trx(larva).yinflection_mm(2:end)-trx(larva).yinflection_mm(1:end-1))./trx(larva).dt;
end

units=parseunits('mm/s');
data=dyinflection_mm;