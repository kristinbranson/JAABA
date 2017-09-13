function  im = addHoGToImage(im,hogs,locs)
  
  iy = size(im,1);
  ix = size(im,2);
  for ndx = 1:numel(hogs)
    hx = size(hogs{ndx},2);
    hy = size(hogs{ndx},1);

    lx = round(locs{ndx}(1));
    ly = round(locs{ndx}(2));
    
    pyy = 1:hy;    pxx = 1:hx;
    pyy( (ly - hy/2 +pyy-1)<1)= [];
    pyy( (ly - hy/2+pyy-1)>iy) = [];
    pxx( (lx - hx/2+pxx-1)<1)= [];
    pxx( (lx - hy/2+pxx-1)>ix) = [];
    
    im( ly - hy/2 + pyy -1,...
            lx - hx/2 + pxx -1) = ...
         double(im( ly - hy/2 + pyy -1,...
            lx - hx/2 + pxx -1)) + ...
          hogs{ndx}(pyy,pxx);
          
  end
  
  


