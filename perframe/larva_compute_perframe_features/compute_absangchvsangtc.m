%compute abs of angchvsangtc

function [data,units]=compute_absangchvsangtc(trx,n)
larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
absangchvstc=cell(1,numlarvae);

for i=1:numlarvae
    larva=larvae(i);
    absangchvstc{1,i}=abs(trx(larva).angchvsangtc);
end

units=parseunits('rad');
data=absangchvstc;
