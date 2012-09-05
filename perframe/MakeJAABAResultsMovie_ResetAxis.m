function ax = MakeJAABAResultsMovie_ResetAxis(x,y,i,nr,nc,pxwidthradius,pxheightradius,pxborder)

n = numel(x);

minx = x(i);
maxx = x(i);
miny = y(i);
maxy = y(i);
w = pxwidthradius-2*pxborder;
h = pxheightradius-2*pxborder;
for j = i+1:n,
  minx = min(minx,x(j));
  maxx = max(maxx,x(j));
  dx = ceil(maxx-minx)+1;
  if dx > w,
    break;
  end
  miny = min(miny,y(j));
  maxy = max(maxy,y(j));
  dy = ceil(maxy-miny)+1;
  if dy > h,
    break;
  end
end
mux = (minx+maxx)/2;
muy = (miny+maxy)/2;
ax = [mux+pxwidthradius*[-1,1],muy+pxheightradius*[-1,1]];

if ax(2)>nc && ax(1)<1,
  ax(1:2) = (nc+1)/2+pxwidthradius*[-1,1];
elseif ax(2)>nc,
  ax(2) = nc;
  ax(1) = nc-(2*pxwidthradius+1)+1;
elseif ax(1)<1,
  ax(1) = 1;
  ax(2) = 2*pxwidthradius+1;
end
if ax(2)>nc || ax(1)<1,
  ax(1:2) = (nc+1)/2+pxwidthradius*[-1,1];
end
if ax(4)>nr && ax(3)<1,
  ax(3:4) = (nr+1)/2+pxwidthradius*[-1,1];
elseif ax(4)>nr,
  ax(4) = nr;
  ax(3) = nr-(2*pxheightradius+1)+1;
elseif ax(3)<1,
  ax(3) = 1;
  ax(4) = 2*pxheightradius+1;
end
if ax(4)>nr || ax(3)<1,
  ax(3:4) = (nr+1)/2+pxheightradius*[-1,1];
end
