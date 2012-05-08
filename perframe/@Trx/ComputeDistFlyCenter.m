function dcenter = ComputeDistFlyCenter(obj,n,fly1)

nflies = obj.nfliespermovie(n);
flyidx1 = obj.getFlyIdx(n,fly1);

% initialize
x_mm1 = obj.GetPerFrameData('x_mm',n,fly1);
y_mm1 = obj.GetPerFrameData('y_mm',n,fly1);
nframes1 = obj.nframes(flyidx1);
firstframe1 = obj.firstframes(flyidx1);
endframe1 = obj.endframes(flyidx1);
dcenter = nan(nflies,nframes1);

for fly2 = 1:nflies,
  if fly2 == fly1, continue; end
  
  flyidx2 = obj.getFlyIdx(n,fly2);
  firstframe2 = obj.firstframes(flyidx2);
  endframe2 = obj.endframes(flyidx2);
  
  % get start and end frames of overlap
  t0 = max(firstframe1,firstframe2);
  t1 = min(endframe1,endframe2);
  
  % no overlap
  if t1 < t0, continue; end
  
  % indices for these frames
  offi = firstframe1-1;
  offj = firstframe2-1;
  i0 = t0 - offi;
  i1 = t1 - offi;
  j0 = t0 - offj;
  j1 = t1 - offj;
  
  x_mm2 = obj.GetPerFrameData('x_mm',n,fly2);
  y_mm2 = obj.GetPerFrameData('y_mm',n,fly2);
  
  % centroid distance
  dx = x_mm2(j0:j1)-x_mm1(i0:i1);
  dy = y_mm2(j0:j1)-y_mm1(i0:i1);
  z = sqrt(dx.^2 + dy.^2);
  dcenter(fly2,i0:i1) = z;
  
end

