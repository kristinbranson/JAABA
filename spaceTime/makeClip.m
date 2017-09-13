function makeClip(bdir,fly,fstart,fend)

moviename = fullfile(bdir,'movie.ufmf');
trackfilename = fullfile(bdir,'trx.mat');
[~,expname] = fileparts(bdir);
params = getParams;
npatches = params.npatches;
psize = params.psize;
nbins = params.nbins; 
patchsz = params.patchsz;
scale = params.scale;


tracks = load(trackfilename);
tracks = tracks.trx;

[readfcn,nframes,fid,headerinfo] = get_readframe_fcn(moviename);

im1 = [];
count = 1;
for fno = fstart:fend
  i1 = readfcn(fno);
  trackndx = fno - tracks(fly).firstframe + 1;
  locy = round(tracks(fly).y(trackndx));
  locx = round(tracks(fly).x(trackndx));
  im1(:,:,count) = extractPatch(i1,...
    locy,locx,tracks(fly).theta(trackndx),patchsz);
  count = count+1;
end

vid = VideoWriter(sprintf('clip_%s_Fly%d_From%d_To%d_%s.avi',expname,fly,fstart,fend,datestr(now,'yyyymmdd')));
open(vid);
for ndx = 1:size(im1,3)
  writeVideo(vid,uint8(imresize(im1(:,:,ndx),6)));
end
close(vid);
