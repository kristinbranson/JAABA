%compute dxhead_mm

function [data,units]=compute_dxhead_mm(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
dxhead_mm=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    dxhead_mm{1,i}=(trx(larva).xhead_mm(2:end)-trx(larva).xhead_mm(1:end-1))./trx(larva).dt;
end

units=parseunits('mm/s');
data=dxhead_mm;