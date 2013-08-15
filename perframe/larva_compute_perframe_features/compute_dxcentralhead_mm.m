%compute dxcentralhead_mm

function [data,units]=compute_dxcentralhead_mm(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
dxcentralhead_mm=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    dxcentralhead_mm{1,i}=(trx(larva).xcentralhead_mm(2:end)-trx(larva).xcentralhead_mm(1:end-1))./trx(larva).dt;
end

units=parseunits('mm/s');
data=dxcentralhead_mm;