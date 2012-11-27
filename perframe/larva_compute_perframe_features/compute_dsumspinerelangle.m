%compute dsumspinerelangle

function [data,units]=compute_dsumspinerelangle(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
dsumspinerelangle=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    dsumspinerelangle{1,i}=(bsxfun(@minus,trx(larva).sumspinerelangle(2:end),trx(larva).sumspinerelangle(1:end-1)))./trx(larva).dt;
end

units=parseunits('rad/s');
data=dsumspinerelangle;