%compute xhead_mm
function [data,units]=compute_xhead_mm(trx,n)
larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
xhead_mm=cell(1,numlarvae);

for i=1:numlarvae
    larva=larvae(i);
    xhead_mm{1,i}=trx(larva).xspine_mm(1,:);
end

data=xhead_mm;
units = parseunits('mm');
