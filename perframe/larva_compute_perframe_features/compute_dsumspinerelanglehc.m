%compute dsumspinerelanglehc

function [data,units]=compute_dsumspinerelanglehc(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
dsumspinerelanglehc=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    dsumspinerelanglehc{1,i}=(bsxfun(@minus,trx(larva).sumspinerelanglehc(2:end),trx(larva).sumspinerelanglehc(1:end-1)))./trx(larva).dt;
end

units=parseunits('rad/s');
data=dsumspinerelanglehc;