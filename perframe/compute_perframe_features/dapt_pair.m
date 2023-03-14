function [dist_apt,dxy_apt,i0,i1,j0,j1] = dapt_pair(trx,fly1,fly2,pt1,pt2,n_pts)
dist_apt = nan(1,trx(fly1).nframes);
dxy_apt = nan(2,trx(fly1).nframes);

% get start and end frames of overlap
t0 = max(trx(fly1).firstframe,trx(fly2).firstframe);
t1 = min(trx(fly1).endframe,trx(fly2).endframe);
i0 = nan;
i1 = nan; 
j0 = nan;
j1 = nan;
% no overlap
if t1 < t0, 
  return;
end
 
% indices for these frames
i0 = t0 + trx(fly1).off;
i1 = t1 + trx(fly1).off;
j0 = t0 + trx(fly2).off;
j1 = t1 + trx(fly2).off;

nframes = trx(fly1).nframes;
apt_data1 = trx(fly1).kpts;
apt_data1 = reshape(apt_data1,n_pts,[],nframes);

nframes2 = trx(fly2).nframes;
apt_data2 = trx(fly2).kpts;
apt_data2 = reshape(apt_data2,n_pts,[],nframes2);

a1 = apt_data1(pt1,:,:);
a2 = apt_data2(pt2,:,:);

if numel(pt2) > 1
  a1 = repmat(a1,[numel(pt2),1,1]);  
end
dd = a2(:,:,j0:j1)-a1(:,:,i0:i1);
d = sqrt(sum(dd.^2,2));
[d,ix] = min(d,[],1,'omitnan');
if size(dd,1)>1
  for ndx = j0:j1
    dxy_apt(:,i0+(ndx-j0)) = dd(ix(ndx-j0+1),:,ndx-j0+1);
  end
else
  dxy_apt(:,i0:i1) = reshape(dd,[size(dd,2),size(dd,3)]);
end
dist_apt(:,i0:i1) = d(1,:);
