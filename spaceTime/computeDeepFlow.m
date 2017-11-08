function [uv, Vs, scale] = computeDeepFlow(im1curr,im2curr,scale)

if nargin <3,
  scale = 4;
end

isz = size(im1curr);
nghsz = 0.23;
% imwrite(repmat(uint8(im2curr),[1 1 3]),'epicflow/deepmatching/im2.png');
% imwrite(repmat(uint8(im1curr),[1 1 3]),'epicflow/deepmatching/im1.png');
% cmd = sprintf('%s %s %s %d %d %s %d %s','epicflow/deepmatching/deepmatching ',...
%   ' epicflow/deepmatching/im1.png epicflow/deepmatching/im2.png',...
%   '-resize ',isz(1)*scale,isz(2)*scale,...
%   ' -ngh_rad ', round(nghsz*isz(1)*scale),...
%   ' -out epicflow/deepmatching/temp.out');
% % fprintf(cmd);fprintf('\n');
% [a,b] = system(cmd);
% if ~a, 
%   error(b); 
% end
% %%
% 
% ofile = 'epicflow/deepmatching/temp.out';
% sz = size(im1curr);
% A = dlmread(ofile);

if ~isdeployed
  A = deepmex(single(im1curr),single(im2curr),...
    scale*isz(1),scale*isz(2),round(nghsz*scale*isz(1)));
else
  A = deepmex_cluster(single(im1curr),single(im2curr),...
    scale*isz(1),scale*isz(2),round(nghsz*scale*isz(1)));
end
A = A';
A(:,1:4) = A(:,1:4)/scale;
Vx = zeros(isz); 
Vy = Vx; Vs = Vx;
idx = sub2ind(isz,round(A(:,1)-0.0006)+1,round(A(:,2)-0.0006)+1);
Vx(idx) = A(:,4)-A(:,2);
Vy(idx) = A(:,3)-A(:,1);
Vs(idx) = A(:,5);
uv = cat(3,Vx,Vy);
