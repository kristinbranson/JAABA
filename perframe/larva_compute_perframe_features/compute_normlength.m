%compute normalization of the spinelength by the mean spinelength size in the whole video
function [data,units]=compute_normlength(trx,n)
larvae= trx.exp2flies{n};
numlarvae = numel(larvae);

meanlength=zeros(1,numlarvae);
numelements=zeros(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    meanlength(i)=nanmean(trx(larva).spinelength(1:end));
    numelements(i)=numel(trx(larva).spinelength(1:end));
end
normvalue=sum(meanlength.*numelements./(sum(numelements)));
normlength=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    normlength{1,i}=trx(larva).spinelength(1:end)/normvalue;
end

units=parseunits([]);
data=normlength;
