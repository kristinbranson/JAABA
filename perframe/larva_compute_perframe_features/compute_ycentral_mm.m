%compute ycentral_mm
function [data,units]=compute_ycentral_mm(trx,n)
larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
ycentral_mm=cell(1,numlarvae);

for i=1:numlarvae
    larva=larvae(i);
    ycentral_mm{1,i}=trx(larva).yspine_mm(6,:);
end

data=ycentral_mm;
units = parseunits('mm');
