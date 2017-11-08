function hhof = plotHofFeatures(H,hax,maxv,colors,binedges,psize,offset)

if nargin<7,
  offset = [0 0];
end

hhof = [];
for xi = 1:size(H,2),
  cx = (psize/2 + 1 + (xi-1)*psize);
  for yi = 1:size(H,1),
    cy = (psize/2 + 1 + (yi-1)*psize);
    
    for bini = 1:size(H,3),
      tmp = linspace(binedges(bini),binedges(bini+1),20);
      xcurr = cx + [0,psize/2*cos(tmp),0];
      ycurr = cy + [0,psize/2*sin(tmp),0];
      hhof(yi,xi,bini) = patch(xcurr+offset(1),ycurr+offset(2),colors(bini,:),...
        'LineStyle','none','FaceAlpha',min(1,H(yi,xi,bini)/maxv),'Parent',hax);
    end
    
  end
end
