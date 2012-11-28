%compute width as the ratio between area and spine length
function [data,units]=compute_ratioarealength(trx,n)
larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
ratioarealength=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    ratioarealength{1,i}=bsxfun(@rdivide,trx(larva).area_mm,trx(larva).spinelength);
    % KB: WAS
    %ratioarealength{1,i}=bsxfun(@rdivide,trx(larva).area,trx(larva).spinelength{1,i});
end


units=parseunits('mm');
data=ratioarealength;
