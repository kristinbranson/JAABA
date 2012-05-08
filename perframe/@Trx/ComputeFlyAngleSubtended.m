function d = ComputeFlyAngleSubtended(obj,n,fly1)

nflies = obj.nfliespermovie(n);
flyidx1 = obj.getFlyIdx(n,fly1);

% initialize
x_mm1 = obj.GetPerFrameData('x_mm',n,fly1);
y_mm1 = obj.GetPerFrameData('y_mm',n,fly1);
a_mm1 = obj.GetPerFrameData('a_mm',n,fly1);
b_mm1 = obj.GetPerFrameData('b_mm',n,fly1);
theta_mm1 = obj.GetPerFrameData('theta_mm',n,fly1);
nframes1 = obj.nframes(flyidx1);
firstframe1 = obj.firstframes(flyidx1);
endframe1 = obj.endframes(flyidx1);
d = nan(nflies,nframes1);

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

  x_mm2 = obj.GetPerFrameData('x_mm',n,fly2);
  y_mm2 = obj.GetPerFrameData('y_mm',n,fly2);
  theta_mm2 = obj.GetPerFrameData('theta_mm',n,fly2);
  a_mm2 = obj.GetPerFrameData('a_mm',n,fly2);
  b_mm2 = obj.GetPerFrameData('b_mm',n,fly2);
  
  % distance from fly1's nose to fly2
  for t = t0:t1,
    i = t - offi;
    j = t - offj;
    
    d(fly2,i) = anglesubtended(...
      x_mm1(i),y_mm1(i),a_mm1(i),b_mm1(i),theta_mm1(i),...
      x_mm2(j),y_mm2(j),a_mm2(j),b_mm2(j),theta_mm2(j),...
      obj.fov);

  end
  
end

