%compute xinflection
function [data,units]=compute_xinflection_mm(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
xinflection_mm=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    sz = size(trx(larva).xspine_mm);
    xinflection_mm{1,i} = trx(larva).xspine_mm(sub2ind(sz,trx(larva).inflectionpointdistance,1:sz(2)));
%     xinflection_mm1{1,i}=zeros(1,size(trx(larva).xspine_mm,2));
%     for j=1:size(trx(larva).xspine_mm,2)
%     xinflection_mm1{1,i}(j)=trx(larva).xspine_mm(trx(larva).inflectionpointdistance(j),j);
%     end
%     if any(xinflection_mm1{i} ~= xinflection_mm{i}),
%       error('mismatch');
%     end
end

units=parseunits('mm');
data=xinflection_mm;
