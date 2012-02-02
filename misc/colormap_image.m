% Irgb = colormap_image(I,c)

function Irgb = colormap_image(I,c,clim)

if nargin < 2,
  c = colormap;
end
if nargin < 3,
  clim = nan(1,2);
end
if ndims(I) ~= 2,
  error('input image must be N x M');
end;
if isnan(clim(1)),
  a = min(I(:));
else
  a = clim(1);
end
if isnan(clim(2)),
  b = max(I(:));
else
  b = clim(2);
end
d = b-a;
n = size(c,1);
nr = size(I,1);
nc = size(I,2);
J = min(n,max(1,round((I - a)/d*(n-1)+1)));
Irgb = reshape(c(J(:),:),[nr,nc,3]);