%compute the xcoordinate of the head in the tails coordinate system

function [data,units]=compute_xtailhead_mm(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
xtailehad_mm=cell(1,numlarvae);

for i=1:numlarvae
    larva=larvae(i);
    xtailehad_mm{1,i}=trx(larva).xhead_mm-trx(larva).xtail_mm;
end

units=parseunits('mm');
data=xtailehad_mm;