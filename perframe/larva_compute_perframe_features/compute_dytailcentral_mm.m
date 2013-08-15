%compute dytailcentral_mm

function [data,units]=compute_dytailcentral_mm(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
dytailcentral_mm=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    dytailcentral_mm{1,i}=(trx(larva).ytailcentral_mm(2:end)-trx(larva).ytailcentral_mm(1:end-1))./trx(larva).dt;
end

units=parseunits('mm/s');
data=dytailcentral_mm;