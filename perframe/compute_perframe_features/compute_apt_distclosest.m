function [data,dxyclosest_all,units] = compute_apt_distclosest(trx,n,apt_info,pt1,pt2)

flies = trx.exp2flies{n};
nflies = numel(flies);
dclosest_all = cell(1,nflies);
dxyclosest_all = cell(1,nflies);

n_parts = apt_info.params.n_classes;

parfor i1 = 1:nflies,
  fly1 = flies(i1);
%  fprintf('fly1 = %d\n',fly1);
  flies2 = flies(trx.roi(fly1)==trx.roi(flies));
  dclosest = nan(numel(flies2),trx(fly1).nframes);  
  dxy = nan(numel(flies2),2,trx(fly1).nframes);
  
  for i2 = 1:numel(flies2)
    fly2 = flies2(i2);
    if fly1 == fly2
      continue;
    end
    [dcur,dxycur] = dapt_pair(trx,fly1,fly2,pt1,pt2,n_parts);
    dclosest(i2,:) = dcur;
    dxy(i2,:,:) = dxycur;
  end
  [dclosest_all{i1},closesti] = min(dclosest,[],1,'omitnan');
  dxycur = nan(2,trx(fly1).nframes);
  for ndx = 1:trx(fly1).nframes
    if isnan(dclosest_all{i1}(ndx)), continue; end
    dxycur(:,ndx) = dxy(closesti(ndx),:,ndx);
  end
  dxyclosest_all{i1} = dxycur;
end

data = dclosest_all;
units = parseunits('unit');
end
