%compute inflection point based on distance
% choose the spine point that is farthest from the line connecting the head
% and the tail
function [data,units]=compute_inflectionpointdistance(trx,n)
larvae=trx.exp2flies{n};
numlarvae=numel(larvae);

inflectionpointdistance=cell(1,numlarvae);
distances=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    % compute the distance from every spine point to the line from the tail
    % to the head of the larva
    m=(trx(larva).yhead_mm-trx(larva).ytail_mm)./(trx(larva).xhead_mm-trx(larva).xtail_mm);
    isvertical = isinf(m);
    n=(1/2)*((trx(larva).yhead_mm+trx(larva).ytail_mm)-m.*(trx(larva).xhead_mm+trx(larva).xtail_mm));
    m=repmat(m,11,1);
    n=repmat(n,11,1);
    distances{1,i}=abs(m.*trx(larva).xspine_mm-trx(larva).yspine_mm+n)./sqrt(m.^2+1);
    % above formula doesn't work for vertical lines
    distances{1,i}(:,isvertical) = abs(bsxfun(@minus,trx(larva).xspine_mm(:,isvertical),trx(larva).xspine_mm(1,isvertical)));
    % note that ties go to the spine point closest to the head
    [dcurr,inflectionpointdistance{1,i}]=max(distances{1,i});
    % sanity check
    isbad = isnan(dcurr) | isinf(dcurr);
    if any(isbad),
      warning('Found %d nan or inf values for inflection point distance, larva %d',nnz(isbad),larva);
    end
end

units=parseunits('units');
data=inflectionpointdistance;
