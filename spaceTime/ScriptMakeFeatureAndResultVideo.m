%% Make an antennal grooming video.

% File locs and other parameters.
expdir = '/home/mayank/Dropbox/ForMayankFromAlice/grooming_GMR_30B01_AE_01_CsChr_RigB_20150903T161828/';
scorefilename = 'scores_antennal.mat';
trxfilename = 'registered_trx.mat';
% fly = 10;
% frames = 5551:5650;
fly = 1;
frames = 6276:6325;
outpath = '.';
outname = sprintf('grooming_result_fly%d_%dto%d.avi',fly,frames(1),frames(end));

%% smooth the trx file

trx = load(fullfile(expdir,trxfilename));
nsize = 30;
for flynum = 1:numel(trx.trx)
  trx.trx(flynum).x = smooth(trx.trx(flynum).x,nsize)';
  trx.trx(flynum).y = smooth(trx.trx(flynum).y,nsize)';
  
end
save(fullfile(expdir,'smoothed_trx.mat'),'-struct','trx');

%%
S = load(fullfile(expdir,scorefilename));
scores = S.allScores.scores{fly}(frames);
scnorm = S.allScores.scoreNorm;
scores = scores/scnorm;
scores(scores>1) = 1;
scores(scores<-1) = -1;

%% grab individual frames overlaid with features.
curtrxfilename = 'smoothed_trx.mat';

hofim = cell(1,numel(frames));
hogim = cell(1,numel(frames));
parfor count = 1:numel(frames)
  fnum = frames(count);
  curhofim = VisualizeFlowFeatures(expdir,fly,fnum,0,'trxfilename',curtrxfilename,'method','hs-sup');
  curhogim = VisualizeHogFeatures(expdir,fly,fnum,'trxfilename',curtrxfilename);
  hofim{count} = curhofim.cdata;
  hogim{count} = curhogim.cdata;
end


%% Make the video.

fig = figure;
params = getParams();

% Amount to crop should be a multiple of params.psize*params.scale 
% to avoid clean cropping of the overlaid hog/hof features.
ycropt = 1*params.psize*params.scale;
ycropb = ycropt;
xcropl = 0*params.psize*params.scale;
xcropr = 1*params.psize*params.scale;

% setting for the score ball and the text string.
ctr = [30,size(hofim{1},2)-xcropl-xcropr];
radius = 21;
txtstrloc = [ctr(1),ctr(2)+radius+5];
numPoints=100; 
theta=linspace(0,2*pi,numPoints);
rho=ones(1,numPoints)*radius;
[X,Y] = pol2cart(theta,rho);


vidobj = VideoWriter(fullfile(outpath,outname));
set(vidobj,'FrameRate',5,'Quality',95);
open(vidobj);
for ndx = 1:numel(frames)
  figure(fig);
  clf(fig);
  hax = axes('Position',[0,0,1,1]);
  set(fig,'Units','pixels');

  imshow(uint8([hofim{ndx}((ycropt+1):(end-ycropb),(xcropl+1):(end-xcropr),:) hogim{ndx}((ycropt+1):(end-ycropb),(xcropl+1):(end-xcropr),:)]));
  hold on;
  if scores(ndx)<0
    cc = [0 0 1];
    tstr = 'No Grooming';
  else
    cc = [1 0 0];
    tstr = 'Antennal Grooming';
  end
  
  patch(X+ctr(2),Y+ctr(1),cc,'LineStyle','none','FaceAlpha',abs(scores(ndx)));
  tobj = text(txtstrloc(2),txtstrloc(1),tstr,'Color',cc,'FontSize',16);
  
  % Uncommenting the below code will align the text such that 
  % the center of the textstring and the score ball are the same.
%   pos = get(tobj,'Position');
%   text_ext = get(tobj,'Extent');
%   tsz = text_ext(3);
%   pos(1) = ctr(2)-tsz/2;
%   set(tobj,'Position',pos);
%   title(sprintf('%d',ndx));

  truesize(fig);
  pause(0.1);
  curfr = getframe(hax);
  % sometimes the sizes of the first fetched frames doesnt match the rest.
  % not sure why.
  
  writeVideo(vidobj,curfr);
end
close(vidobj);
close(fig)
