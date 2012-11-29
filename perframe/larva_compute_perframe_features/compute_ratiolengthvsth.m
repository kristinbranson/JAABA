%compute ratio lengthvsth
function [data,units]=compute_ratiolengthvsth(trx,n)

larvae = trx.exp2flies{n};
numlarvae = numel(larvae);
ratiolengthvsth=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    ratiolengthvsth{1,i}=trx(larva).spinelength./trx(larva).tailheadmag;
end

units=parseunits('unit');
data=ratiolengthvsth;
