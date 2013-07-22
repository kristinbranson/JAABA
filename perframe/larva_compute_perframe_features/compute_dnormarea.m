%compute dnormarea
function [data,units]=compute_dnormarea(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
dnormarea=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    dnormarea{1,i}=(trx(larva).normarea(2:end)-trx(larva).normarea(1:end-1))./trx(larva).dt;
end

units=parseunits('1/s');
data=dnormarea;
