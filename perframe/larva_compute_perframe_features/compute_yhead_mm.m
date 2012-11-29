%compute yhead_mm
function [data,units]=compute_yhead_mm(trx,n)
larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
yhead_mm=cell(1,numlarvae);

for i=1:numlarvae
    larva=larvae(i);
    yhead_mm{1,i}=trx(larva).yspine_mm(1,:);
end

data=yhead_mm;
units = parseunits('mm');
