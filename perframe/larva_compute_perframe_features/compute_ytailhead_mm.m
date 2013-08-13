%compute the ycoordinate of the head in the tails coordinate system

function [data,units]=compute_ytailhead_mm(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
ytailhead_mm=cell(1,numlarvae);

for i=1:numlarvae
    larva=larvae(i);
    ytailhead_mm{1,i}=trx(larva).yhead_mm-trx(larva).ytail_mm;
end

units=parseunits('mm');
data=ytailhead_mm;