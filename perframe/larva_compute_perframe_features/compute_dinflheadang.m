%compute dinflheadang
function [data,units]=compute_dinflheadang(trx,n)


larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
%absdtailheadang=cell(1,numlarvae);
dinflheadang=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    dinflheadang{1,i}=modrange(diff(trx(larva).inflheadang),-pi,pi)./trx(larva).dt;
    % KB: this angle distance didn't make sense to me, changed to one that
    % made sense to me
    % dinflheadang{1,i}=(mod1(trx(larva).inflheadang(2:end)-trx(larva).inflheadang(1:end-1),pi)-pi)./trx(larva).dt;
end

units=parseunits('rad/s');
data=dinflheadang;
