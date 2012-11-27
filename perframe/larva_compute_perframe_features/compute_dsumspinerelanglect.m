%compute dsumspinerelanglect

function [data,units]=compute_dsumspinerelanglect(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
dsumspinerelanglect=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    dsumspinerelanglect{1,i}=(bsxfun(@minus,trx(larva).sumspinerelanglect(2:end),trx(larva).sumspinerelanglect(1:end-1)))./trx(larva).dt;
end

units=parseunits('rad/s');
data=dsumspinerelanglect;