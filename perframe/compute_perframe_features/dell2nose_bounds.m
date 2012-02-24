function [mindupper,dlower] = dell2nose_bounds(trx,fly1,flies)

nflies = numel(flies);
mindupper = inf(1,trx(fly1).nframes);
dlower = nan(nflies,trx(fly1).nframes);
for i2 = 1:nflies,
  fly2 = flies(i2);
  if fly1 == fly2,
    continue;
  end
  [dcenter2nose,off10,off11] = dcenter2nose_pair(trx,fly1,fly2);
  if off10 <= off11,
    mindupper(off10:off11) = min(mindupper(off10:off11),dcenter2nose(off10:off11) + 2*trx(fly1).a_mm(off10:off11));
    dlower(i2,off10:off11) = dcenter2nose(off10:off11) - 2*trx(fly1).a_mm(off10:off11);
  end
end
