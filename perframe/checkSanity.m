function checkSanity(fastY,res_dumb,r,off,funcType,trans_type,extraStr)
  if nargin<6
    extraStr = '';
  end
  if any(isnan(fastY) ~= isnan(res_dumb)),
    fprintf('SANITY CHECK: %s, trans = %s, r = %d, off = %d, %s nan mismatch\n',funcType,trans_type,r,off,extraStr);
  else
    fprintf('SANITY CHECK: mean, trans = %s, r = %d, off = %d, %s max error = %f\n',funcType,trans_type,r,off,max(abs(fastY-res_dumb)),extraStr);
  end
