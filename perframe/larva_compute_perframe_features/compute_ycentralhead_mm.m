%compute the ycoordinate of the head in the central coordinate system

function [data,units]=compute_ycentralhead_mm(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
ycentralhead_mm=cell(1,numlarvae);

for i=1:numlarvae
    larva=larvae(i);
    ycentralhead_mm{1,i}=trx(larva).yhead_mm-trx(larva).ycentral_mm;
end

units=parseunits('mm');
data=ycentralhead_mm;