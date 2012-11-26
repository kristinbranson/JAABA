%compute dyhead_mm

function [data,units]=compute_dyhead_mm(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
dyhead_mm=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    dyhead_mm{1,i}=(trx(larva).yhead_mm(2:end)-trx(larva).yhead_mm(1:end-1))./trx(larva).dt;
end

units=parseunits('mm/s');
data=dyhead_mm;