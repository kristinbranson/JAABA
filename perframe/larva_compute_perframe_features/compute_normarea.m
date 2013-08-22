%compute normalization of the area by the mean area size in the whole video
function [data,units]=compute_normarea(trx,n)
larvae= trx.exp2flies{n};
numlarvae = numel(larvae);

meanarea=zeros(1,numlarvae);
numelements=zeros(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    meanarea(i)=nanmean(trx(larva).area_mm(1:end));
    numelements(i)=numel(trx(larva).area_mm(1:end));

end
normvalue=sum(meanarea.*numelements./(sum(numelements)));
normarea=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    normarea{1,i}=trx(larva).area_mm(1:end)/normvalue;
end

units=parseunits([]);
data=normarea;
