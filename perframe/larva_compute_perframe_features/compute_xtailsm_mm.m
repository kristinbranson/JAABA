%compute xtailsm_mm
function [data,units]=compute_xtailsm_mm(trx,n)
larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
xtailsm_mm=cell(1,numlarvae);

for i=1:numlarvae
    larva=larvae(i);
    xtailsm_mm{1,i}=trx(larva).xmaxcurvstop_mm;
end

data=xtailsm_mm;
units = parseunits('mm');
