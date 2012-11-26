%compute length of the spine 

function [data,units]=compute_spinelength(trx,n)

larvae= trx.exp2flies{n};
numlarvae = numel(larvae);
spinereldist=cell(1,numlarvae);
spinelength=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    spinereldist{1,i}=bsxfun(@hypot,trx(larva).xspine_mm(1:end-1,:)-trx(larva).xspine_mm(2:end,:),trx(larva).yspine_mm(1:end-1,:)-trx(larva).yspine_mm(2:end,:));
    spinelength{1,i}=sum(spinereldist{1,i},1);
end


units=parseunits('mm');
data=spinelength;
