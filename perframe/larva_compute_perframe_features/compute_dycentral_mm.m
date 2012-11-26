%compute dycentral_mm

function [data,units]=compute_dycentral_mm(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
dycentral_mm=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    dycentral_mm{1,i}=(trx(larva).ycentral_mm(2:end)-trx(larva).ycentral_mm(1:end-1))./trx(larva).dt;
end

units=parseunits('mm/s');
data=dycentral_mm;