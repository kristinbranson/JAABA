%compute dratiolengthvsth
function [data,units]=compute_dratiolengthvsth(trx,n)

larvae = trx.exp2flies{n};
numlarvae = numel(larvae);
dratiolengthvsth=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    dratiolengthvsth{1,i}=diff(trx(larva).ratiolengthvsth)./trx(larva).dt;
end

units=parseunits('unit/s');
data=dratiolengthvsth;