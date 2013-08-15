%compute the ycoordinate of the center in the tails coordinate system

function [data,units]=compute_ytailcentral_mm(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
ytailcentral_mm=cell(1,numlarvae);

for i=1:numlarvae
    larva=larvae(i);
    ytailcentral_mm{1,i}=trx(larva).ycentral_mm-trx(larva).ytail_mm;
end

units=parseunits('mm');
data=ytailcentral_mm;
