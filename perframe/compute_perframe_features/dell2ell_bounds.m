%AR 3/18/2018

function [mindupper,dlower] = dell2ell_bounds(trx,fly1,flies)

nflies = numel(flies);
mindupper = inf(1,trx(fly1).nframes);
dlower = nan(nflies,trx(fly1).nframes);
for i2 = 1:nflies
    fly2 = flies(i2);
    if fly1 == fly2
        continue;
    end
    [dcenter,off10,off11,off20,off21] = dcenter_pair(trx,fly1,fly2);
    if off10 <= off11
        % upper bound is circle with x,y center of fly 1 and radius  = dcenter to fly2 + 1/2 bodylength of fly2
        % ideally dcenter would be dell2ctr but this is not calculated
        mindupper(off10:off11) = min(mindupper(off10:off11),dcenter(off10:off11) + 2*trx(fly2).a_mm(off20:off21));
        % lower bound is circle with x,y center of fly 1 and radius = dcenter to fly2 -1/2 bodylength of fly 2
        % ideally dcenter would be dell2ctr but this is not calculated 
        dlower(i2,off10:off11) = dcenter(off10:off11) - 2*trx(fly2).a_mm(off20:off21);
    end
end
