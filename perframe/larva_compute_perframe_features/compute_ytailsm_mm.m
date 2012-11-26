%compute ytailsm_mm
function [data,units]=compute_ytailsm_mm(trx,n)
larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
ytailsm_mm=cell(1,numlarvae);

for i=1:numlarvae
    larva=larvae(i);
    ytailsm_mm{1,i}=trx(larva).ymaxcurvstop_mm;
end

data=ytailsm_mm;
units = parseunits('mm');