%compute dxtailcentral_mm

function [data,units]=compute_dxtailcentral_mm(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
dxtailcentral_mm=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    dxtailcentral_mm{1,i}=(trx(larva).xtailcentral_mm(2:end)-trx(larva).xtailcentral_mm(1:end-1))./trx(larva).dt;
end

units=parseunits('mm/s');
data=dxtailcentral_mm;