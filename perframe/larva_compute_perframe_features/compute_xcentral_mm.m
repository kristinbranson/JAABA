%compute xcentral_mm
function [data,units]=compute_xcentral_mm(trx,n)
larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
xcentral_mm=cell(1,numlarvae);

for i=1:numlarvae
    larva=larvae(i);
    xcentral_mm{1,i}=trx(larva).xspine_mm(6,:);
end

data=xcentral_mm;
units = parseunits('mm');
