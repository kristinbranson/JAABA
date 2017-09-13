%% create and testing a ufmf movie.
% Last checked there was off-by-one error.

hh = ufmf_read_header('../mated10_20140714T131113/movie.ufmf');
outfile = '../temp/test.ufmf';
fig = figure('Visible','off');
fid = fopen(outfile,'w');

% Setup the ufmf.
if strcmp(hh.coding,'MONO8'),
  ncolors = 1;
end

gui = get(fig,'UserData');
[im,hh1,timestamp,bb,mu] = ufmf_read_frame(hh,1);

gui.bg.params.UpdatePeriod = inf;
gui.bg.params.NFrames = inf;
gui.bg.params.boxBasedCompression = true;
gui.bg.params.KeyframePeriod = inf;

nroi = 1;
roi = [1 1 hh.max_height hh.max_width];
gui.bg.lastkeyframetime = -inf;
tmp = struct('loc',cast([],'int64'),'timestamp',[]);
index = struct;
index.frame = tmp;
index.keyframe.mean = tmp;
gui.bg.index = repmat(index,[1,nroi]);
% (re)initialize the bg models
if ~isfield(gui.bg,'model') || length(gui.bg.model) ~= nroi,  
  gui.bg.model = struct('roi',cell(1,nroi),'mu',cell(1,nroi),'nframes',cell(1,nroi));
  for i = 1:nroi,
    gui.bg.model(i) = struct('roi',roi(i,:),...
      'mu',mu,...
      'nframes',0);
%      'mu',zeros([roi(i,4),roi(i,3),ncolors],'uint8'),...
  end
end
gui.bg.lastupdatetime = -inf;
gui.isBGModel = true;
set(fig,'UserData',gui);

writeUFMFHeader(fid,fig,1);

% Write few frames
writeUFMFKeyFrame(fig,timestamp,fid,1);

for ndx = 1:10
  [im,hh1,timestamp,bb] = ufmf_read_frame(hh,ndx);
  %bb = bb(:,[2 1 4 3]);
  mywriteUFMFFrame(fig,im,timestamp,fid,1,bb);
  
end

% wrap up the ufmf


wrapupUFMF(fid,fig,1);
fclose(fid);
fclose(hh.fid);
close(fig);
% read the created ufmf

hh1 = ufmf_read_header(outfile);
hh = ufmf_read_header('../mated10_20140714T131113/movie.ufmf');

%
figure;
for ndx = 1:10
  [im1,hh2,timestamp,bb] = ufmf_read_frame(hh1,ndx);
  subplot(1,3,1); imshow(im1);
  [im2,hh2,timestamp,bb] = ufmf_read_frame(hh,ndx);
  subplot(1,3,2); imshow(im2);
  dd = abs(double(im1)-double(im2));
  subplot(1,3,3); imshow(uint8(dd));
  title(sprintf('%d %d',max(dd(:)),nnz(dd(:))));
  pause;
end
fclose all
