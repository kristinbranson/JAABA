function [Vx,Vy] = computeFlowBkgSup(im1curr,im2curr,params)
% function [Vx,Vy] = computeFlowBkgSup(im1curr,im2curr,params)
% Computes flow and then suppresses flow in the background.

%There is some flow towards the center that can't be estimated
% properly. The flow increases as we go away from the center.
% d_err is to account for that.
gparams = getParams;

sz = round(size(im1curr));
bwimg = zeros(sz(1),sz(2));
ctr = [ceil( (sz(1)+1)/2),ceil( (sz(2)+1)/2)];
bwimg(ctr(1),ctr(2))=1;
dimg = bwdist(bwimg,'euclidean');
[xx,yy]= meshgrid(1:sz(2),1:sz(1));
aimg = atan2(-(yy-ctr(1)),-(xx-ctr(2)));

dd_err = dimg/size(im1curr,1);
flow_thres = gparams.flow_thres;

uv = estimate_flow_interface(im1curr,im2curr,...
  'hs-brightness',{'max_warping_iters',gparams.warping_iters});

if ~params.stationary,
  Vx = uv(:,:,1);
  Vy = uv(:,:,2);
  return;
end

cdx = params.dx;
cdy = params.dy;
curt = params.theta;
rotd = [cos(curt) sin(curt); -sin(curt) cos(curt)]*[cdx;cdy];
cdx = -rotd(1); cdy = -rotd(2);

ctheta = params.dtheta;
rotflowu = dimg.*(cos(aimg+ctheta)-cos(aimg));
rotflowv = dimg.*(sin(aimg+ctheta)-sin(aimg));


dd1 = sqrt( (uv(:,:,1)-cdx-rotflowu).^2 + (uv(:,:,2)-cdy-rotflowv).^2);

crdx = params.rdx;  crdy = params.rdy;
rotd = [cos(curt) sin(curt); -sin(curt) cos(curt)]*[crdx;crdy];
fly_flow = rotd;

dd2 = sqrt( (uv(:,:,1)-fly_flow(1)).^2 + (uv(:,:,2)-fly_flow(2)).^2);
dd = min(dd1,dd2);
for ndx = 1:2
  tt = uv(:,:,ndx);
  tt( dd<(dd_err+flow_thres)) = 0;
  uv(:,:,ndx) = tt;
end
Vx = uv(:,:,1);
Vy = uv(:,:,2);
