%compute dytailhead_mm

function [data,units]=compute_dytailhead_mm(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
dytailhead_mm=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    dytailhead_mm{1,i}=(trx(larva).ytailhead_mm(2:end)-trx(larva).ytailhead_mm(1:end-1))./trx(larva).dt;
end

units=parseunits('mm/s');
data=dytailhead_mm;