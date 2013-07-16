function [binidx, nbins, featurenames, featureboundaries, featurecenters] = compute_spacetime_mask(meana, meanb)

boxwidth2 = round(meanb*6);
boxheight2 = round(meana*3);

nbinsr = 4;
nbinstheta = 8;

as = linspace(meana,boxheight2,nbinsr);
bs = linspace(meanb,boxwidth2,nbinsr);

[xgrid,ygrid] = meshgrid(-boxwidth2:boxwidth2,-boxheight2:boxheight2);
rbin = zeros(size(xgrid));
for i = 1:nbinsr,
  idx = rbin == 0 & xgrid.^2/bs(i)^2 + ygrid.^2/as(i)^2 <= 1;
  rbin(idx) = i;
end
thetalims = [-pi/2-pi/nbinstheta,3*pi/2-pi/nbinstheta];
thetagrid = modrange(atan2(ygrid,xgrid),thetalims(1),thetalims(2));
dtheta = 2*pi/nbinstheta;

thetabin = min(floor((thetagrid-thetalims(1))/dtheta) + 1,nbinstheta);
binidx = zeros(size(xgrid));
binidx(rbin>0) = sub2ind([nbinsr+1,nbinstheta],rbin(rbin>0),thetabin(rbin>0));
binidx(rbin==1) = 1;
[binis,~,binidx] = unique(binidx);
binidx(binidx==find(binis==0)) = 0;
binidx = binidx-1;
binis(binis==0) = [];
binidx = reshape(binidx,size(xgrid));
nbins = numel(binis);

featurenames = cell(1,nbins);
featurenames{1} = 'theta0_r0';
for i = 2:nbins,
  [r,t] = ind2sub([nbinsr+1,nbinstheta],binis(i));
  r = r-1;
  %featurenames{i} = sprintf('theta%d_r%d',round((thetalims(1)+dtheta*(t-.5))*180/pi),r);
  featurenames{i} = sprintf('theta%d_r%d',t,r);
end

for i = 1:nbins,
  col = min(find(sum(binidx==i)));
  row = min(find(binidx(:,col)==i));
  featureboundaries{i} = bwtraceboundary(binidx==i,[row col],'N');  
  x = mean(find(sum(binidx==i,1)));
  y = mean(find(sum(binidx==i,2)));
  featurecenters{i} = [x y];
end

%   subplot(1,4,2);
%   imagesc([-boxwidth2,boxwidth2],[-boxheight2,boxheight2], binidx);
%   axis image;
%   for i = 1:nbins,
%     x = mean(xgrid(binidx==i));
%     y = mean(ygrid(binidx==i));
%     text(x,y,featurenames{i},'Interpreter','none','HorizontalAlignment','center');
%   end