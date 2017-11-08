function [Vx,Vy] = computeDeepFlowBkgSup(im1curr,im2curr,params)
% function [Vx,Vy] = computeFlowBkgSup(im1curr,im2curr,params)
% Computes flow and then suppresses flow in the background.

%There is some flow towards the center that can't be estimated
% properly. The flow increases as we go away from the center.
% d_err is to account for that.
global_params = getParams;
% flow_thres = 1;
% deepscale = 4;
stepsz = global_params.deepmatching_step/global_params.deepscale;

sz = round(size(im1curr)/stepsz);
bwimg = zeros(sz(1),sz(2));
ctr = [ceil( (sz(1)+1)/2),ceil( (sz(2)+1)/2)];
bwimg(ctr(1),ctr(2))=1;
dimg = bwdist(bwimg,'euclidean')*stepsz;
[xx,yy]= meshgrid(1:sz(2),1:sz(1));
aimg = atan2(-(yy-ctr(1)),-(xx-ctr(2)));
dd_err = dimg/size(im1curr,1)*stepsz;


[uv,Vs] = computeDeepFlow(im1curr,im2curr,global_params.deepscale);
uv = uv(2:stepsz:end,2:stepsz:end,:);
Vss = Vs(2:stepsz:end,2:stepsz:end);
selpx = Vss<global_params.deep_thres;
im1sm = imresize(im1curr,1/stepsz);
im2sm = imresize(im2curr,1/stepsz);

uvhs = estimate_flow_interface(im1sm,im2sm,...
  'hs-brightness',{'max_warping_iters',global_params.warping_iters});
uvhs = uvhs*stepsz;
% Multiply by stepsz to get the flow in the same range.

for ndx = 1:2
  tt = uv(:,:,ndx);
  tths = uvhs(:,:,ndx);
  tt(selpx) = tths(selpx);
  uv(:,:,ndx) = tt;
end


if ~params.stationary,
  uv = imresize(uv,stepsz);
  Vx = uv(:,:,1);
  Vy = uv(:,:,2);
  return;
end

% Account for background movement.
% flies movement by [dx,dy] is rotated by theta 
% because we rotate the patch after we grab the patch
% to make the fly face forward.

cdx = params.dx;
cdy = params.dy;
curt = params.theta;
rotd = [cos(curt) sin(curt); -sin(curt) cos(curt)]*[cdx;cdy];
cdx = -rotd(1); cdy = -rotd(2);

ctheta = params.dtheta;
rotflowu = dimg.*(cos(aimg+ctheta)-cos(aimg));
rotflowv = dimg.*(sin(aimg+ctheta)-sin(aimg));

dd1 = sqrt( (uv(:,:,1)-cdx-rotflowu).^2 + (uv(:,:,2)-cdy-rotflowv).^2);


% Account for motion within the fly.
% this occurs because we round off fly's location
% when we extract the patch.
crdx = params.rdx;  crdy = params.rdy;
rotd = [cos(curt) sin(curt); -sin(curt) cos(curt)]*[crdx;crdy];
fly_flow = rotd;

dd2 = sqrt( (uv(:,:,1)-fly_flow(1)).^2 + (uv(:,:,2)-fly_flow(2)).^2);
dd = min(dd1,dd2);
for ndx = 1:2
  tt = uv(:,:,ndx);
  tt( dd<(dd_err+global_params.flow_thres)) = 0;
  uv(:,:,ndx) = tt;
end

% Resize at the very end.
% Resizing it before computing background flow 
% doesn't suppress the background enough because of
% deep flows patchy output.
uv = imresize(uv,stepsz);

Vx = uv(:,:,1);
Vy = uv(:,:,2);
