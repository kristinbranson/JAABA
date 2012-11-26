%compute dtailcentralmag
function [data,units]=compute_dtailcentralmag(trx,n)


larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
dtailcentralmag=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    dtailcentralmag{1,i}=(trx(larva).tailcentralmag(2:end)-trx(larva).tailcentralmag(1:end-1))./trx(larva).dt;
end

units=parseunits('mm/s');
data=dtailcentralmag;
