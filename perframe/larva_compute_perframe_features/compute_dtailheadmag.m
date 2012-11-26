%compute dtailheadmag
function [data,units]=compute_dtailheadmag(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
dtailheadmag=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    dtailheadmag{1,i}=(trx(larva).tailheadmag(2:end)-trx(larva).tailheadmag(1:end-1))./trx(larva).dt;
end

units=parseunits('mm/s');
data=dtailheadmag;
