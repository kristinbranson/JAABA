function MakeJAABAResultsMovie_PositiveClips(expdir,scorefn,varargin)

%%

textcolor = [0,0,0];
fontsize = 14;
trx_linewidth = 1;
target_linewidth = 2;
trxcolor = [0,0,0];
labelcolor = [.7,0,0];
trx_rad = 25;
compression = 'None';
allowedcompressions = {'Indeo3', 'Indeo5', 'Cinepak', 'MSVC', 'RLE', 'None','Uncompressed AVI','Motion JPEG AVI'};
useVideoWriter = exist('VideoWriter','file');
aviname = '';
fps = 10;
usemencoder = true;

[nzoomr,nzoomc,...
  figpos,...
  moviefilestr,trxfilestr,...
  hfig,...
  nframesout,...
  weight_id,weight_tmid,...
  weight_length,weight_sumscore,...
  weight_meanscore,...
  min_sumscore,...
  min_length,...
  pxwidthradius,...
  pxheightradius,...
  awidthradius,...
  aheightradius,...
  pxborder,...
  aborder,...
  nframesborder,...
  labelcolor,...
  useVideoWriter,...
  avifileTempDataFile,aviname,fps,...
  usemencoder,...
  behavior] ...
  = myparse(varargin,...
  'nzoomr',4,'nzoomc',5,...
  'figpos',[50,50,1500,1000],...
  'moviefilestr','movie.ufmf',...
  'trxfilestr','registered_trx.mat',...
  'hfig',1001,...
  'nframesout',500,...
  'weight_id',100,...
  'weight_tmid',1,...
  'weight_length',1,...
  'weight_sumscore',1,...
  'weight_meanscore',.01,...
  'min_sumscore',0,...
  'min_length',0,...
  'pxwidthradius',[],...
  'pxheightradius',[],...
  'awidthradius',3,...
  'aheightradius',3,...
  'pxborder',5,...
  'aborder',1,...
  'nframesborder',5,...
  'labelcolor',labelcolor,...
  'useVideoWriter',useVideoWriter,...
  'avifileTempDataFile','',...
  'aviname',aviname,...
  'fps',fps,...
  'usemencoder',usemencoder,...
  'behavior','');

TMID = 1;
LENGTH = 2;
SUMSCORE = 3;
MEANSCORE = 4;

w = nan(4,1);
w(TMID) = weight_tmid;
w(LENGTH) = weight_length;
w(SUMSCORE) = weight_sumscore;
w(MEANSCORE) = weight_meanscore;


if ~ischar(compression),
  compression = '';
end
if ~ispc && ~useVideoWriter,
  compression = 'None';
end
if ~isempty(compression) && ~any(strcmpi(compression,allowedcompressions)),
  fprintf('Unknown compressor %s\n',compression);
  compression = '';
end

if isempty(aviname),
  [~,experiment_name] = fileparts(expdir);
  if isempty(behavior),
    tmp = regexprep(scorefn,'\.mat$','');
    tmp = regexprep(tmp,'[sS]cores_?','');
    behavior = regexprep(tmp,'[sS]core_?','');
  end
  aviname = fullfile(expdir,sprintf('JAABA%sExamples_%s.avi',behavior,experiment_name));
end


%% open files
moviename = fullfile(expdir,moviefilestr);
[readframe,nframes,fid] = get_readframe_fcn(moviename);
im = readframe(1);
[nr,nc,ncolors] = size(im);

trxname = fullfile(expdir,trxfilestr);
[trx,trxname,loadsucceeded,timestamps] = load_tracks(trxname);
if ~loadsucceeded,
  error('Could not load trx from file %s',trxname);
end
nids = length(trx);
trxendframes = [trx.endframe];
trxfirstframes = [trx.firstframe];
T0 = min(trxfirstframes);
T1 = max(trxendframes);

scorefilename = fullfile(expdir,scorefn);
tmp = load(scorefilename);
scores = cell(1,nids);
labels = cell(1,nids);
for j = 1:nids,
  scores{j} = tmp.allScores.scores{j}(trxfirstframes(j):trxendframes(j));
  labels{j} = scores{j} > 0;
end

%% compute statistics of bouts

X = [];
ids = [];
allt0s = [];
allt1s = [];
alllengths = [];

for i = 1:nids,
  idx = scores{i} > 0;
  [t0s,t1s] = get_interval_ends(idx);
  t1s = t1s - 1;
  sumscores = [0,cumsum(scores{i})];
  allt0s = [allt0s,t0s]; %#ok<AGROW>
  allt1s = [allt1s,t1s]; %#ok<AGROW>
  ids = [ids,i+zeros(size(t0s))]; %#ok<AGROW>
  xcurr = nan(4,numel(t0s));
  xcurr(TMID,:) = (t0s+t1s)/2;
  xcurr(LENGTH,:) = t1s-t0s+1;
  alllengths = [alllengths,xcurr(LENGTH,:)]; %#ok<AGROW>
  xcurr(SUMSCORE,:) = sumscores(t1s+1)-sumscores(t0s);
  xcurr(MEANSCORE,:) = xcurr(SUMSCORE,:) ./ xcurr(LENGTH,:);
  X = [X,xcurr]; %#ok<AGROW>
end
badidx = X(SUMSCORE,:) < min_sumscore | X(LENGTH,:) < min_length;
X(:,badidx) = [];
allt0s(badidx) = [];
allt1s(badidx) = [];
ids(badidx) = [];
nbouts = size(X,2);

% get the rank for each of these stats
for i = 1:size(X,1),
  [~,order] = sort(X(i,:));
  [~,order] = sort(order);
  X(i,:) = order / nbouts; %#ok<AGROW>
end

ischosen = false(1,nbouts);
did = ones(1,nbouts);
dx = ones(size(X));
nframesfilled = zeros([nzoomr,nzoomc]);
boutis = cell(nzoomr,nzoomc);
while true,
  
  % choose the most different example
  d = weight_id*did + sum(bsxfun(@times,w,dx),1);
  d(ischosen) = -1;
  maxd = max(d);
  idxcurr = find(d == maxd);
  if numel(idxcurr) > 1,
    i = randsample(idxcurr,1);
  else
    i = idxcurr;
  end
  
  if ischosen(i),
    break;
  end
  
  % update distances
  did(ids == ids(i)) = 0;
  dx = min(dx,abs(bsxfun(@minus,X,X(:,i))));

  % try to put this somewhere
  [nfilledcurr,j] = min(nframesfilled(:));
  
  if nfilledcurr >= nframesout,
    break;
  end
  
  t0 = max(1,allt0s(i)-nframesborder);
  t1 = min(trx(ids(i)).nframes,allt1s(i)+nframesborder);
  n = t1-t0+1;
  
  [r,c] = ind2sub([nzoomr,nzoomc],j);
  boutis{r,c}(end+1) = i;
  nframesfilled(j) = nframesfilled(j) + n;
  ischosen(i) = true;
  
end

nframesout = min(nframesout,max(nframesfilled(:)));


tsout = nan(nzoomr,nzoomc,nframesout);
idsout = nan(nzoomr,nzoomc,nframesout);

for r = 1:nzoomr,
  for c = 1:nzoomc,
    off = 0;
    for i = 1:numel(boutis{r,c}),
      j = boutis{r,c}(i);
      t0 = max(1,allt0s(j)-nframesborder) + trx(ids(j)).firstframe - 1;
      t1 = min(trx(ids(j)).nframes,allt1s(j)+nframesborder) + trx(ids(j)).firstframe - 1;
      n = t1-t0+1;
      idsout(r,c,off+1:off+n) = ids(j);
      tsout(r,c,off+1:off+n) = t0:t1;
      off = off + n;
    end
  end
end

%% output window sizes


a = [trx.a];
meana = nanmedian(a)*4;
if isempty(pxwidthradius),
  pxwidthradius = awidthradius*meana;
end
if isempty(pxheightradius),
  pxheightradius = aheightradius*meana;
end
if isempty(pxborder),
  pxborder = aborder*meana;
end

axh = figpos(4)/nzoomr;
axw = figpos(3)/nzoomc;
pxheightradius1 = min((pxwidthradius)*axh / axw,(nr-1)/2);
pxwidthradius1 = min((pxheightradius)*axw / axh,(nc-1)/2);
pxheightradius = max(pxheightradius1,pxheightradius);
pxwidthradius = max(pxwidthradius1,pxwidthradius);
ax = MakeJAABAResultsMovie_ResetAxis(nc/2,nr/2,1,nr,nc,pxwidthradius,pxheightradius,pxborder);
axround = max(1,[min(nc,[floor(ax(1)),ceil(ax(2))]),...
  min(nr,[floor(ax(3)),ceil(ax(4))])]);
ax = repmat(ax(:),[1,nzoomr,nzoomc]);
axround = repmat(axround(:),[1,nzoomr,nzoomc]);

%% create the figure

figure(hfig);
clf;
set(hfig,'Position',figpos);

hax = createsubplots(nzoomr,nzoomc,0);
hax = reshape(hax,[nzoomr,nzoomc]);
hinfo = nan(nzoomr,nzoomc);
him = nan(nzoomr,nzoomc);
htrx = nan(nzoomr,nzoomc);
htarget = nan(nzoomr,nzoomc);
hlabel = nan(nzoomr,nzoomc);
for i = 1:nzoomr*nzoomc,
  him(i) = imagesc(0,'Parent',hax(i),[0,255]);
  axis(hax(i),'image');
  axis(hax(i),ax(:,i)');
  hold(hax(i),'on');
  htrx(i) = plot(hax(i),nan,nan,'.-','LineWidth',trx_linewidth,'color',trxcolor);
  hlabel(i) = plot(hax(i),nan,nan,'.','LineWidth',trx_linewidth,'color',labelcolor);
  htarget(i) = plot(hax(i),nan,nan,'-','LineWidth',target_linewidth,'color','k');
  hinfo(i) = text(ax(1,i),ax(4,i),'','Parent',hax(i),'Color',textcolor,...
    'FontUnits','pixels','FontSize',fontsize,...
    'HorizontalAlignment','left',...
    'VerticalAlignment','bottom');
  
  set(hax(i),'XTickLabel',{},'YTickLabel',{});
end
colormap gray;

%% 

isfirstframe = true;

for j = 1:nframesout,
  for r = 1:nzoomr,
    for c = 1:nzoomc,
      t = tsout(r,c,j);
      if isnan(t),
        set(him(r,c),'CData',0);
        set(htrx(r,c),'XData',nan,'YData',nan);
        set(hlabel(r,c),'XData',nan,'YData',nan);
        set(htarget(r,c),'XData',nan,'YData',nan);
        set(hinfo(r,c),'String','');
        continue;
      end
      fly = idsout(r,c,j);
      i = t+trx(fly).off;
      im = readframe(t);
    
      i0 = max(1,min(trx(fly).nframes,i-trx_rad));
      i1 = max(1,min(trx(fly).nframes,i+trx_rad));
      tmpx = trx(fly).x(i0:i1);
      tmpy = trx(fly).y(i0:i1);
    
      set(htrx(r,c),'XData',tmpx,...
        'YData',tmpy);
      set(hlabel(r,c),'XData',tmpx(labels{fly}(i0:i1)),...
        'YData',tmpy(labels{fly}(i0:i1)));
    
      updatefly(htarget(r,c),trx(fly).x(i),trx(fly).y(i),...
        trx(fly).theta(i),trx(fly).a(i),trx(fly).b(i));
      if labels{fly}(i),
        set(htarget(r,c),'Color',labelcolor);
      else
        set(htarget(r,c),'Color',trxcolor);
      end
      
      outofbounds = trx(fly).x(i)-pxborder < ax(1,r,c) || ...
        trx(fly).x(i)+pxborder > ax(2,r,c) || ...
        trx(fly).y(i)-pxborder < ax(3,r,c) || ...
        trx(fly).y(i)+pxborder > ax(4,r,c);

      if j == 1 || tsout(r,c,j)-1 ~= tsout(r,c,j-1) || outofbounds,
        ax(:,r,c) = MakeJAABAResultsMovie_ResetAxis(trx(fly).x,trx(fly).y,i,nr,nc,pxwidthradius,pxheightradius,pxborder);
        axis(hax(r,c),ax(:,r,c)');
        set(hinfo(r,c),'Position',ax([1,4]));
        axround(:,r,c) = max(1,[min(nc,[floor(ax(1,r,c)),ceil(ax(2,r,c))]),...
          min(nr,[floor(ax(3,r,c)),ceil(ax(4,r,c))])]);
        set(hinfo(r,c),'Position',ax([1,4],r,c));
      end
      set(him(r,c),'CData',im(axround(3,r,c):axround(4,r,c),axround(1,r,c):axround(2,r,c)),...
        'XData',axround([1,2],r,c),'YData',axround([3,4],r,c));
      set(hinfo(r,c),'String',sprintf('Target %d, t = %fs',fly,timestamps(t)));
    end
  end
  
  
  if isfirstframe,
    if useVideoWriter,
      if strcmpi(compression,'None') || strcmpi(compression,'Uncompressed AVI'),
        profile = 'Uncompressed AVI';
      else
        profile = 'Motion JPEG AVI';
      end
      aviobj = VideoWriter(aviname,profile); %#ok<TNMLP>
      set(aviobj,'FrameRate',fps);
      if ~strcmpi(profile,'Uncompressed AVI'),
        set(aviobj,'Quality',100);
      end
      open(aviobj);
    else
      if isempty(avifileTempDataFile),
        aviobj = avifile(aviname,'fps',fps,'quality',100,'compression',compression);  %#ok<REMFF1>
      else
        aviobj = myavifile(aviname,'fps',fps,'quality',100,'compression',compression,...
          'TempDataFile',avifileTempDataFile);
        fprintf('Temporary data file for avi writing: %s\n',aviobj.TempDataFile);
      end
    end
      
    isfirstframe = false;
  end

  fr = getframe(hfig);
  height = size(fr.cdata,1);
  width = size(fr.cdata,2);
  if useVideoWriter,
    writeVideo(aviobj,fr);
  else
    aviobj = addframe(aviobj,fr);
  end
  
  
  pos = get(hfig,'Position');
  if any(pos(3:4)~=figpos(3:4)),
    pos(3:4) = figpos(3:4);
    set(hfig,'Position',pos);
  end
  
  drawnow;
end


%% clean up
  
if useVideoWriter,
  close(aviobj);
else
  aviobj = close(aviobj); %#ok<NASGU>
end

fclose(fid);

%% compress using mencoder

if isunix && usemencoder,
  [path,base,ext] = fileparts(aviname);
  tmpfile = fullfile(path,[base,'_tmp',ext]);
  newheight = 4*ceil(height/4);
  newwidth = 4*ceil(width/4);
  % subtitles are upside down, so encode with subtitles and flip, then flip
  % again
  cmd = sprintf('mencoder %s -o %s -ovc xvid -xvidencopts fixed_quant=4 -vf scale=%d:%d',...
    aviname,tmpfile,newwidth,newheight);
  status = system(cmd);
  if status ~= 0,
    error('Failed to compress avi to %s',tmpfile);
  end
  unix(sprintf('mv %s %s',tmpfile,aviname));
end

