function [mindupper,dlower] = dnose2ell_bounds(trx,fly1,flies)

nflies = numel(flies);
mindupper = inf(1,trx(fly1).nframes);
dlower = nan(nflies,trx(fly1).nframes);
for i2 = 1:nflies,
  fly2 = flies(i2);
  if fly1 == fly2,
    continue;
  end
  [dnose2center,off10,off11,off20,off21] = dnose2center_pair(trx,fly1,fly2);
  if off10 <= off11,
    mindupper(off10:off11) = min(mindupper(off10:off11),dnose2center(off10:off11) + 2*trx(fly2).a_mm(off20:off21));
    dlower(i2,off10:off11) = dnose2center(off10:off11) - 2*trx(fly2).a_mm(off20:off21);
  end
end
