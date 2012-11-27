%compute sumspinerelanglehc

function [data,units]=compute_sumspinerelanglehc(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
spineangles=cell(1,numlarvae);
spinerelangle=cell(1,numlarvae);
sumspinerelanglehc=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    nspinepts = size(trx(larva).yspine_mm,1);
    spineangles{1,i}=bsxfun(@atan2,trx(larva).yspine_mm(1:end-1,:)-trx(larva).yspine_mm(2:end,:),trx(larva).xspine_mm(1:end-1,:)-trx(larva).xspine_mm(2:end,:));
    %spinerelangle{1,i}=bsxfun(@atan2,sin(spineangles{1,i}(2:end,:))-sin(spineangles{1,i}(1:end-1,:)),cos(spineangles{1,i}(2:end,:))-cos(spineangles{1,i}(1:end-1,:)));
    % KB: this method of measuring the difference between angles didn't
    % make sense to me. 
    spinerelangle{1,i} = modrange(-diff(spineangles{i},1,1),-pi,pi);
    sumspinerelanglehc{1,i}=sum(spinerelangle{1,i}(1:floor(nspinepts/2)-1,:),1);
end

units=parseunits('rad');
data=sumspinerelanglehc;