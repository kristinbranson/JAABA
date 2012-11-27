%compute ytail_mm
function [data,units]=compute_ytail_mm(trx,n)
larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
ytail_mm=cell(1,numlarvae);

for i=1:numlarvae
    larva=larvae(i);
    ytail_mm{1,i}=trx(i).yspine_mm(end,:);
end

data=ytail_mm;
units = parseunits('mm');
