%compute dnormlength
function [data,units]=compute_dnormlength(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
dnormlength=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    dnormlength{1,i}=(trx(larva).normlength(2:end)-trx(larva).normlength(1:end-1))./trx(larva).dt;
end

units=parseunits('1/s');
data=dnormlength;
