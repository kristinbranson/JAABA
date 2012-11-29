%compute yinflection
function [data,units]=compute_yinflection_mm(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
yinflection_mm=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    sz = size(trx(larva).yspine_mm);
    yinflection_mm{1,i} = trx(larva).yspine_mm(sub2ind(sz,trx(larva).inflectionpointdistance,1:sz(2)));
% 
%     yinflection_mm{1,i}=zeros(1,size(trx(larva).yspine_mm,2));
%     for j=1:size(trx(larva).yspine_mm,2)
%     yinflection_mm{1,i}(j)=trx(larva).yspine_mm(trx(larva).inflectionpointdistance(j),j);
%     end
end

units=parseunits('mm');
data=yinflection_mm;
