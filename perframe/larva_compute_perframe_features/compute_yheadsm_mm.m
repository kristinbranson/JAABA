%compute yheadsm_mm
function [data,units]=compute_yheadsm_mm(trx,n)
larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
yheadsm_mm=cell(1,numlarvae);

for i=1:numlarvae
    larva=larvae(i);
    yheadsm_mm{1,i}=trx(larva).ymaxcurvbegin_mm;
end

data=yheadsm_mm;
units = parseunits('mm');