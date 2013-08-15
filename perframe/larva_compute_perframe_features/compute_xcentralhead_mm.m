%compute the xcoordinate of the head in the central coordinate system

function [data,units]=compute_xcentralhead_mm(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
xcentralhead_mm=cell(1,numlarvae);

for i=1:numlarvae
    larva=larvae(i);
    xcentralhead_mm{1,i}=trx(larva).xhead_mm-trx(larva).xcentral_mm;
end

units=parseunits('mm');
data=xcentralhead_mm;