%compute dxtailhead_mm

function [data,units]=compute_dxtailhead_mm(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
dxtailhead_mm=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    dxtailhead_mm{1,i}=(trx(larva).xtailhead_mm(2:end)-trx(larva).xtailhead_mm(1:end-1))./trx(larva).dt;
end

units=parseunits('mm/s');
data=dxtailhead_mm;