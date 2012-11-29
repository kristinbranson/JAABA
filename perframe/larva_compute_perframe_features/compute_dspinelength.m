%compute dspinelength
function [data,units]=compute_dspinelength(trx,n)

larvae = trx.exp2flies{n};
numlarvae = numel(larvae);
dspinelength=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    dspinelength{1,i}=(trx(larva).spinelength(2:end)-trx(larva).spinelength(1:end-1))./trx(larva).dt;
end


units=parseunits('mm/s');
data=dspinelength;
