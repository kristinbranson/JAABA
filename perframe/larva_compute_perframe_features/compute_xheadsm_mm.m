%compute xheadsm_mm
function [data,units]=compute_xheadsm_mm(trx,n)
larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
xheadsm_mm=cell(1,numlarvae);

for i=1:numlarvae
    larva=larvae(i);
    xheadsm_mm{1,i}=trx(larva).xmaxcurvbegin_mm;
end

data=xheadsm_mm;
units = parseunits('mm');
