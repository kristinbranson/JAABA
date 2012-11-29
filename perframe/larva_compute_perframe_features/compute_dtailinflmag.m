%compute dtailinflmag
function [data,units]=compute_dtailinflmag(trx,n)

% 
larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
dtailinflmag=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    dtailinflmag{1,i}=(trx(larva).tailinflmag(2:end)-trx(larva).tailinflmag(1:end-1))./trx(larva).dt;
end

units=parseunits('mm/s');
data=dtailinflmag;
