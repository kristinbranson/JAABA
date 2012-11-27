%compute xtail_mm
function [data,units]=compute_xtail_mm(trx,n)
larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
xtail_mm=cell(1,numlarvae);

for i=1:numlarvae
    larva=larvae(i);
    xtail_mm{1,i}=trx(larva).xspine_mm(end,:);
end

data=xtail_mm;
units = parseunits('mm');
