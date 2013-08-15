%compute dycentralhead_mm

function [data,units]=compute_dycentralhead_mm(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
dycentralhead_mm=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    dycentralhead_mm{1,i}=(trx(larva).ycentralhead_mm(2:end)-trx(larva).ycentralhead_mm(1:end-1))./trx(larva).dt;
end

units=parseunits('mm/s');
data=dycentralhead_mm;