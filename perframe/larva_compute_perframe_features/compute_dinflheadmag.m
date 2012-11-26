%compute dinflheadmag
function [data,units]=compute_dinflheadmag(trx,n)


larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
dinflheadmag=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    dinflheadmag{1,i}=(trx(larva).inflheadmag(2:end)-trx(larva).inflheadmag(1:end-1))./trx(larva).dt;
end

units=parseunits('mm/s');
data=dinflheadmag;
