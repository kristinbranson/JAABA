%compute dratioarealength

function [data,units]=compute_dratioarealength(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
dratioarealength=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    dratioarealength{1,i}=(trx(larva).ratioarealength(2:end)-trx(larva).ratioarealength(1:end-1))./trx(larva).dt;
end

units=parseunits('mm/s');
data=dratioarealength;