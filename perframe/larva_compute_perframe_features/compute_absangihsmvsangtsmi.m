%compute absang between vector tailsminfl and vector inflheadsm

function [data,units]=compute_absangihsmvsangtsmi(trx,n)


larvae = trx.exp2flies{n};
numlarvae = numel(larvae);
absangihsmvsangtsmi=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    absangihsmvsangtsmi{1,i}=real(acos(cos(trx(larva).tailsminflang).*cos(trx(larva).inflheadsmang)+sin(trx(larva).tailsminflang).*sin(trx(larva).inflheadsmang)));

end

units=parseunits('rad');
data=absangihsmvsangtsmi;