ofile = 'epicflow/deepmatching/temp.out';
sz = size(im1curr);
A = dlmread(ofile);
Vx = zeros(sz); 
Vy = Vx; Vs = Vx;
idx = sub2ind(sz,round(A(:,2)-0.0006)+1,round(A(:,1)-0.0006)+1);
Vx(idx) = A(:,3)-A(:,1);
Vy(idx) = A(:,4)-A(:,2);
Vs(idx) = A(:,5);
uv = cat(3,Vx,Vy);
figure;  ax = [];
ns = 5;
ax(1) = subplot(1,ns,1);
imshowpair(im1curr,im2curr);
% hold on;
% quiver(Vx,Vy,0);

uv = uv(2:2:end,2:2:end,:);
Vss = Vs(2:2:end,2:2:end);
selpx = Vss<4.2;
% selpx = Vs<4;
tic;
uvold = estimate_flow_interface(im1curr,im2curr,'hs-brightness',...
  {'max_warping_iter',2});
toc;
uvoldss = uvold(2:2:end,2:2:end,:);
uv(selpx)=uvoldss(selpx);
uv = imresize(uv,2);
ax(2) = subplot(1,ns,2);
imshow(flowToColor(uv));
ax(3) = subplot(1,ns,3);
imagesc(Vs); axis image;
ax(4) = subplot(1,ns,4);
imshow(flowToColor(uvold)); 
ax(5) = subplot(1,ns,5);
imshow(flowToColor(cat(3,Vx,Vy))); 

linkaxes(ax);

