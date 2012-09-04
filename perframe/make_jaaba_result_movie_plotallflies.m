% [succeeded,aviname] = make_ctrax_result_movie('param',value,...)
% Options:
% 'moviename': name of raw movie to annotate
% 'trxname': name of mat file containing trajectories, or struct of
% trajectories
% 'aviname': name of movie to output to
% 'colors': colors to plot each fly
% 'nzoomr': number of rows of zoom fly boxes
% 'nzoomc': number of columns of zoom fly boxes
% 'boxradius': radius of the zoom fly box in pixels
% 'taillength': number of frames of trajectory to plot behind each fly
% 'zoomflies': flies to zoom in on
% 'fps': frames per second of output movie
% 'maxnframes': number of frames to output
% 'firstframe': first frame to output
% 'compression': compressor to use when outputting (pc only, will be set to
% 'none' for linux). 
% 'figpos': position of figure
% if any parameters are not given, the user will be prompted for these
function [succeeded,aviname,figpos,height,width] = make_jaaba_result_movie_plotallflies(expdir,classifierparamsfiles,varargin)

succeeded = false;
defaults.boxradius = 1.5;
defaults.taillength = 30;
defaults.fps = 20;
defaults.zoomflies = [];
defaults.nzoomr = 5;
defaults.nzoomc = 3;
defaults.compression = 'None';
allowedcompressions = {'Indeo3', 'Indeo5', 'Cinepak', 'MSVC', 'RLE', 'None','Uncompressed AVI','Motion JPEG AVI'};
useVideoWriter = exist('VideoWriter','file');
mencoderoptions = '';
mencoder_maxnframes = inf;
[moviefilestr,trxfilestr,aviname,colors,behavior_cm,...
  zoomflies,nzoomr,nzoomc,boxradius,...
  taillength,fps,maxnframes,firstframes,compression,figpos,movietitle,...
  useVideoWriter,mencoderoptions,mencoder_maxnframes,...
  avifileTempDataFile,titletext,dynamicflyselection,...
  doshowsex,doplotwings,...
  classifierparamsfiles,...
  behaviors,...
  scorefns,...
  score_norm_prctile,...
  behaviorweight,...
  printnone] = ...
  myparse(varargin,'moviefilestr','movie.ufmf','trxfilestr','registered_trx.mat',...
  'aviname','','colors',[],'behavior_cm',@jet,...
  'zoomflies',[],'nzoomr',defaults.nzoomr,'nzoomc',defaults.nzoomc,...
  'boxradius',nan,'taillength',defaults.taillength,'fps',defaults.fps,...
  'maxnframes',inf,'firstframes',1,'compression',defaults.compression,...
  'figpos',[50,50,1500,1000],'movietitle','','useVideoWriter',useVideoWriter,...
  'mencoderoptions',mencoderoptions,'mencoder_maxnframes',mencoder_maxnframes,...
  'avifileTempDataFile','',...
  'titletext',true,...
  'dynamicflyselection',true,...
  'doshowsex',true,...
  'doplotwings',true,...
  'classifierparamsfiles',{},...
  'behaviors',{},...
  'scorefns',{},...
  'score_norm_prctile',80,...
  'behaviorweight',[],...
  'printnone',false);

if ~ischar(compression),
  compression = '';
end
if ~isempty(compression) && ~any(strcmpi(compression,allowedcompressions)),
  fprintf('Unknown compressor %s\n',compression);
  compression = '';
end

if isempty(behaviors) || isempty(scorefns),
  for i = 1:numel(classifierparamsfiles),
    classifierparams = ReadClassifierParamsFile(classifierparamsfile);
    for j = 1:numel(classifierparams),
      if iscell(classifierparams(j).behaviors.names),
        behaviors{end+1} = sprintf('%s_',classifierparams(j).behaviors.names{:}); %#ok<AGROW>
        behaviors{end} = behaviors{end}(1:end-1);
      else
        behaviors{end+1} = classifierparams(j).behaviors.names; %#ok<AGROW>
      end
      behaviors{end} = regexprep(behaviors{i},'(_|^)([a-z])','${upper($2)}');
      scorefns{end+1} = classifierparams(j).file.scorefilename; %#ok<AGROW>
    end
  end
end
nbehaviors = numel(behaviors);

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
endframes = min(T1,firstframes+maxnframes-1);

scores = cell(1,nids);
for j = 1:nids,
  scores{j} = nan(nbehaviors,trx(j).nframes);
end
for i = 1:nbehaviors,
  scorefilename = fullfile(expdir,scorefns{i});
  tmp = load(scorefilename);
  for j = 1:nids,
    scores{j}(i,:) = tmp.allScores.scores{j}(trxfirstframes(j):trxendframes(j));
  end
end


%% output avi name
haveaviname = false;
if ischar(aviname) && ~isempty(aviname)
  [~,ext] = splitext(aviname);
  if strcmpi(ext,'.avi'),
    haveaviname = true;
  end
end
if ~haveaviname,
  fprintf('Choose avi file to output annotated version of %s\n',moviename);
  [base,~] = splitext(movienameonly);
  aviname = [moviepath,'ctraxresults_',base,'.avi'];
  helpmsg = {};
  helpmsg{1} = 'Choose avi file to write annotated movie to.';
  helpmsg{2} = sprintf('Raw movie input: %s',moviename);
  helpmsg{3} = sprintf('Trx file name input: %s',trxname);
  [avinameonly,avipath] = uiputfilehelp('*.mat',sprintf('Choose output avi for %s',movienameonly),aviname,'helpmsg',helpmsg);
  if ~ischar(avinameonly),
    return;
  end
  aviname = [avipath,avinameonly];
end

%% more parameter parsing

doshowsex = doshowsex && isfield(trx,'sex') && ~any(cellfun(@isempty,{trx.sex}));
if doshowsex,
  sexes = {};
  for i = 1:numel(trx),
    sexes = union(sexes,unique(trx(i).sex));
  end
  sexes = upper(sexes);
  sexes = sort(sexes);
  sexmarkers = {'none','*','x','o','+','d','s','p','h'};
  if numel(sexes) > numel(sexmarkers),
    sexmarkers = repmat(sexmarkers,[1,ceil(numel(sexes)/numel(sexmarkers))]);
  end
  sexmarkers = sexmarkers(1:numel(sexes));
end

doplotwings = doplotwings && all(isfield(trx,{'xwingl','ywingl','xwingr','ywingr'}));

if isnan(boxradius),
  meana = 4*nanmean([trx.a]);
  boxradius = ceil(defaults.boxradius*meana);
end

%% which flies to zoom in on

if ~isempty(zoomflies),
  nzoom = numel(zoomflies);
  if max(zoomflies) > nids || min(zoomflies) < 1,
    error('Illegal values for zoomflies');
  end
  nzoomc = ceil(nzoom/nzoomr);
else
  nzoom = nzoomr*nzoomc;
end

% choose some random flies to zoom in on
nframesoverlap = zeros(1,numel(trx));
for i = 1:numel(firstframes),
  nframesoverlap = nframesoverlap + ...
      max(0,min(endframes(i),trxendframes)-max(firstframes(i),trxfirstframes) + 1);
end
fliesmaybeplot = find(nframesoverlap > 0);

if isempty(zoomflies),
  if numel(fliesmaybeplot) < nzoom,
    zoomflies = [fliesmaybeplot,nan(1,nzoom-length(fliesmaybeplot))];
    fprintf('Not enough flies to plot\n');
  else
    
    if dynamicflyselection,
      
      % choose flies to plot in each frame
      fliesplotperframe = nan(nzoom,nframes);
      
      % flies currently chosen
      fliesplotcurr = nan(nzoom,1);
      allfliesplotted = [];
      
      % loop through all frames
      
      for i = 1:numel(firstframes),
        for f = firstframes(i):endframes(i),
          
          % flies that are currently alive
          isalive = trxfirstframes <= f & trxendframes >= f;
          fliesalive = find(isalive);
          
          % remove newly dead flies
          openzoomboxes = isnan(fliesplotcurr) | ~ismember(fliesplotcurr,fliesalive);
          fliesplotcurr(openzoomboxes) = nan;
          
          % how many flies do we need to choose
          nflieschoose = nnz(openzoomboxes);
          
          % flies we can choose to add
          fliesleft = setdiff(fliesalive,fliesplotcurr);
          
          nfliesshort = max(0,nflieschoose-numel(fliesleft));
          if nfliesshort > 0,
            newfliesplot = fliesleft;
          else
            [~,order] = sort(-[trx(fliesleft).endframe]);
            newfliesplot = fliesleft(order(1:nflieschoose));
          end
          allfliesplotted = [allfliesplotted,newfliesplot]; %#ok<AGROW>
          fliesplotcurr(openzoomboxes) = [newfliesplot,nan(1,nfliesshort)];
          fliesplotperframe(:,f) = fliesplotcurr;
          
        end
      end

      zoomflies = fliesplotperframe(:,firstframes(1));
      
    else
      fliesmaybeplot = fliesmaybeplot(randperm(length(fliesmaybeplot)));
      [~,flieswithmostframes] = sort(-nframesoverlap(fliesmaybeplot));
      zoomflies = sort(fliesmaybeplot(flieswithmostframes(1:nzoom)));
    end
  end
elseif nzoom > length(zoomflies),
  zoomflies = [zoomflies,nan(1,nzoom-length(zoomflies))];
end
zoomflies = reshape(zoomflies,[nzoomr,nzoomc]);
rowszoom = floor(nr/nzoomr);

%% choose colors for each behavior

nadd = nbehaviors - size(colors,1);
if nadd > 0,
  tmp = behavior_cm(256);
  colors = cat(1,colors,tmp(round(linspace(1,256,nadd)),:));
end

%% get the score normalizations

score_norms = nan(nbehaviors,1);
for i = 1:nbehaviors,
  tmp = [];
  for j = 1:nids,
    tmp = [tmp,scores{j}(i,:)]; %#ok<AGROW>
  end
  score_norms(i) = prctile(abs(tmp),score_norm_prctile);
end

for fly = 1:nids,
  for behaviori = 1:nbehaviors,
    scores{fly}(behaviori,:) = min(1,max(-1,scores{fly}(behaviori,:)/score_norms(behaviori)));
  end
end

%% choose weights for each behavior

if isempty(behaviorweight),
  behaviorweight = nan(nbehaviors,1);
  for behaviori = 1:nbehaviors,
    tmp = [];
    for j = 1:nids,
      tmp = [tmp,scores{j}(i,:)]; %#ok<AGROW>
    end
    behaviorweight(behaviori) = nnz(tmp>0);
  end
  behaviorweight = behaviorweight / sum(behaviorweight);
end
if numel(behaviorweight) == 1,
  behaviorweight = ones(nbehaviors,1);
end
behaviorweight = behaviorweight(:);

%% initialize plots

figure(1);
clf;
hold on;
hax = gca;
set(hax,'position',[0,0,1,1]);
axis off;
isdisplay = ispc || ~strcmpi(get(1,'XDisplay'),'nodisplay');

% corners of zoom boxes in plotted image coords
x0 = nc+(0:nzoomc-1)*rowszoom+1;
y0 = (0:nzoomr-1)*rowszoom+1;
x1 = x0 + rowszoom - 1;
y1 = y0 + rowszoom - 1;

% relative frame offset
nframesoff = trxfirstframes - 1;

% pre-allocate
himzoom = zeros(nzoomr,nzoomc);
htail = zeros(2*taillength+1,nids);
htri = zeros(1,nids);
hwing = zeros(1,nids);
hsexmarker = zeros(1,nids);
scalefactor = rowszoom / (2*boxradius+1);
hzoom = zeros(nzoomr,nzoomc);
hzoomwing = zeros(nzoomr,nzoomc);
htextzoom = zeros(nzoomr,nzoomc);
hzoombox = zeros(nzoomr,nzoomc);

% initialize plots
him = image([1,nc],[1,nr],zeros([nr,nc,3]));
axis image;
axis([.5,x1(end)+.5,.5,y1(end)+.5]);
axis off;

if titletext && isdisplay,
  htext = text(.5,.5,'frame info','Parent',hax,'BackgroundColor','k','Color','g','VerticalAlignment','bottom','interpreter','none');
end

for i = 1:nzoomr,
  for j = 1:nzoomc,
    himzoom(i,j) = image([x0(j),x1(j)],[y0(i),y1(i)],repmat(uint8(123),[boxradius*2+1,boxradius*2+1,3]));
  end
end

for fly = 1:nids,
  htail(:,fly) = plot(nan(2,2*taillength+1),nan(2,2*taillength+1),'k-');
  if doshowsex,
    hsexmarker(fly) = plot(nan,nan,'k.');
  end
  htri(fly) = plot(nan,nan,'k-');
  if doplotwings,
    hwing(fly) = plot(nan,nan,'k.-');
  end
end

for i = 1:nzoomr,
  for j = 1:nzoomc,
    hzoom(i,j) = plot(nan,nan,'k-');
    hzoomwing(i,j) = plot(nan,nan,'k.-');
    xtmp = x0(j) + [0,1]*(x1(j)-x0(j)) + [1,-1];
    ytmp = y0(i) + [0,1]*(y1(i)-y0(i)) + [1,-1];
    hzoombox(i,j) = plot(xtmp([1,2,2,1,1]),ytmp([1,1,2,2,1]),'k-','LineWidth',2,...
      'Visible','on');
    htextzoom(i,j) = text((x0(j)+x1(j))/2,.95*y0(i)+.05*y1(i),'',...
      'color','k','horizontalalignment','center',...
      'verticalalignment','bottom','fontweight','bold');
  end
end

mencoder_nframes = 0;

tailx = nan(2*taillength+1,nids);
taily = nan(2*taillength+1,nids);
tailcolor = zeros(2*taillength+1,3,nids);
tailisbehavior = false(2*taillength+1,nids);


%%

for segi = 1:numel(firstframes),
  firstframe = firstframes(segi);
  endframe = endframes(segi);

  frame = firstframe - 1;
  idx = frame - nframesoff;
  tailend = 1;
  for fly = 1:nids,
    tailx(:,fly) = padgrab(trx(fly).x,nan,1,1,idx(fly)-taillength,idx(fly)+taillength);
    taily(:,fly) = padgrab(trx(fly).y,nan,1,1,idx(fly)-taillength,idx(fly)+taillength);
    tmpscores = padgrab(scores{fly},0,1,nbehaviors,idx(fly)-taillength,idx(fly)+taillength);
    tailisbehavior(:,fly) = any(tmpscores > 0,1);
    w = max(0,bsxfun(@times,tmpscores,behaviorweight));
    w = permute(w,[1,3,2]);
    z = sum(w,1);
    z(z==0) = 1;
    w = bsxfun(@rdivide,w,z);
    tailcolortmp = sum(bsxfun(@times,w,colors),1);
    tailcolor(:,:,fly) = permute(tailcolortmp,[3,2,1]);
    for i = 1:taillength*2,
      set(htail(i,fly),'XData',tailx([i,i+1],fly),'YData',taily([i,i+1],fly),'Color',tailcolor(i,:,fly));
      if tailisbehavior(i,fly),
        set(htail(i,fly),'LineWidth',2);
      else
        set(htail(i,fly),'LineWidth',1);
      end
    end
  end
  
  for frame = firstframe:endframe,
    if mod(frame - firstframe,20) == 0,
      fprintf('frame %d\n',frame);
    end
      
    % relative frame
    idx = frame - nframesoff;
    
    isalive = frame >= trxfirstframes & ...
      frame <= trxendframes;
    
    % draw the unzoomed image
    im = uint8(readframe(frame));
    if ncolors == 1,
      im = repmat(im,[1,1,3]);
    end
    set(him,'cdata',im);
    
    % draw frame number text box
    % text doesn't show up in no display mode
    if titletext && isdisplay,
      framestr = sprintf('Frame %d, t = %.2f s',frame,timestamps(frame)-timestamps(1));
      if ~isempty(movietitle),
        framestr = {framestr,movietitle}; %#ok<AGROW>
      end
      set(htext,'String',framestr);
    end
    
    % draw the zoomed image
    if dynamicflyselection && exist('fliesplotperframe','var'),
      zoomflies = reshape(fliesplotperframe(:,frame),[nzoomr,nzoomc]);
    end
    for i = 1:nzoomr,
      for j = 1:nzoomc,
        fly = zoomflies(i,j);
        
        % fly not visible?
        if isnan(fly) || ~isalive(fly),
          set(himzoom(i,j),'cdata',repmat(uint8(123),[boxradius*2+1,boxradius*2+1,3]));
          continue;
        end
        
        % grab a box around (x,y)
        x = round(trx(fly).x(idx(fly)));
        y = round(trx(fly).y(idx(fly)));
        boxradx1 = min(boxradius,x-1);
        boxradx2 = min(boxradius,size(im,2)-x);
        boxrady1 = min(boxradius,y-1);
        boxrady2 = min(boxradius,size(im,1)-y);
        box = uint8(zeros(2*boxradius+1));
        box(boxradius+1-boxrady1:boxradius+1+boxrady2,...
          boxradius+1-boxradx1:boxradius+1+boxradx2) = ...
          im(y-boxrady1:y+boxrady2,x-boxradx1:x+boxradx2);
        set(himzoom(i,j),'cdata',repmat(box,[1,1,3]));
        
      end
    end
    
    idxfuture = frame + taillength - nframesoff;
    colorfuture = zeros(nids,3);
    isalivefuture = frame + taillength >= trxfirstframes & ...
      frame + taillength <= trxendframes;

    tailx = cat(1,tailx(2:end,:),nan(1,nids));
    taily = cat(1,taily(2:end,:),nan(1,nids));
    tailisbehavior = cat(1,tailisbehavior(2:end,:),false(1,nids));
    for fly = 1:nids,
      if ~isalivefuture(fly), continue; end
      w = max(0,scores{fly}(:,idxfuture(fly))).*behaviorweight;
      tailisbehavior(end,fly) = any(scores{fly}(:,idxfuture(fly))>0,1);
      z = sum(w);
      if z == 0,
        colorfuture(fly,:) = [0,0,0];
      else
        w = w /z;
        %colortmp(fly,:) = sum(bsxfun(@times,w,colors),1)/z*max(scores{fly}(:,idx(fly)));
        colorfuture(fly,:) = sum(bsxfun(@times,w,colors),1);
      end
      tailx(end,fly) = trx(fly).x(idxfuture(fly));
      taily(end,fly) = trx(fly).y(idxfuture(fly));
    end
    
    tailcolor = cat(1,tailcolor(2:end,:,:),permute(colorfuture,[3,2,1]));
    colortmp = permute(tailcolor(taillength+1,:,:),[3,2,1]);
    isbehavior = tailisbehavior(taillength+1,:);
    
    % sanity check
    if debug,
      tmpcolortmp = zeros(nids,3);
      tmptailx = nan(2*taillength+1,nids);
      tmptaily = nan(2*taillength+1,nids);
      tmpisbehavior = false(1,nids);
      for fly = 1:nids,
        if ~isalive(fly), continue; end
        w = max(0,scores{fly}(:,idx(fly))).*behaviorweight;
        z = sum(w);
        if z == 0,
          tmpcolortmp(fly,:) = [0,0,0];
        else
          w = w /z;
          %colortmp(fly,:) = sum(bsxfun(@times,w,colors),1)/z*max(scores{fly}(:,idx(fly)));
          tmpcolortmp(fly,:) = sum(bsxfun(@times,w,colors),1);
        end
        tmptailx(:,fly) = padgrab(trx(fly).x,nan,1,1,idx(fly)-taillength,idx(fly)+taillength);
        tmptaily(:,fly) = padgrab(trx(fly).y,nan,1,1,idx(fly)-taillength,idx(fly)+taillength);
        
        tmpisbehavior(fly) = any(scores{fly}(:,idx(fly)) > 0);
        
      end
      if ~all(isnan(tailx(:)) & isnan(tmptailx(:)) | (tailx(:) == tmptailx(:)))
        error('tailx does not pass sanity check');
      end
      if ~all(isnan(taily(:)) & isnan(tmptaily(:)) | (taily(:) == tmptaily(:)))
        error('taily does not pass sanity check');
      end
      if ~all(colortmp(:) == tmpcolortmp(:)),
        error('colortmp does not pass sanity check');
      end
      if ~all(isbehavior(:) == tmpisbehavior(:)),
        error('isbehavior does not pass sanity check');
      end
    end
    
    % plot the zoomed out position
    for fly = 1:nids,
      if isalive(fly),
        set(htail(tailend,fly),'XData',tailx(end-1:end,fly),...
          'YData',taily(end-1:end,fly),...
          'Color',tailcolor(end-1,:,fly));
        if tailisbehavior(end-1,fly),
          set(htail(tailend,fly),'LineWidth',2);
        else
          set(htail(tailend,fly),'LineWidth',1);
        end
%         i0 = max(1,idx(fly)-taillength);
%         set(htail(fly),'xdata',trx(fly).x(i0:idx(fly)),...
%           'ydata',trx(fly).y(i0:idx(fly)),'color',colortmp(fly,:));
        if doshowsex,
          sexi = find(strcmpi(trx(fly).sex{idx(fly)},sexes),1);
          sexmarker = sexmarkers{sexi};
          set(hsexmarker(fly),'xdata',trx(fly).x(idx(fly)),...
            'ydata',trx(fly).y(idx(fly)),...
            'color',colortmp(fly,:),...
            'marker',sexmarker,...
            'markerfacecolor',colortmp(fly,:));
        end
        updatefly(htri(fly),trx(fly),idx(fly));
        set(htri(fly),'Color',colortmp(fly,:));
        if isbehavior(fly),
          set(htri(fly),'Linewidth',2);
        else
          set(htri(fly),'Linewidth',1);
        end
        if doplotwings,
          xwing = [trx(fly).xwingl(idx(fly)),trx(fly).x(idx(fly)),trx(fly).xwingr(idx(fly))];
          ywing = [trx(fly).ywingl(idx(fly)),trx(fly).y(idx(fly)),trx(fly).ywingr(idx(fly))];
          set(hwing(fly),'XData',xwing,'YData',ywing,'Color',colortmp(fly,:));
        end
      else
        set(htail(:,fly),'xdata',nan,'ydata',nan);
        set(htri(fly),'xdata',[],'ydata',[]);
        if doshowsex,
          set(hsexmarker(fly),'xdata',[],'ydata',[]);
        end
        if doplotwings,
          set(hwing(fly),'XData',[],'YData',[]);
        end
      end
    end
    tailend = tailend + 1;
    if tailend >= 2*taillength+1,
      tailend = 1;
    end
    
    % plot the zoomed views
    for i = 1:nzoomr,
      for j = 1:nzoomc,
        fly = zoomflies(i,j);
        if ~isnan(fly) && isalive(fly),
          x = trx(fly).x(idx(fly));
          y = trx(fly).y(idx(fly));
          x = boxradius + (x - round(x))+.5;
          y = boxradius + (y - round(y))+.5;
          x = x * scalefactor;
          y = y * scalefactor;
          x = x + x0(j) - 1;
          y = y + y0(i) - 1;
          a = trx(fly).a(idx(fly))*scalefactor;
          b = trx(fly).b(idx(fly))*scalefactor;
          theta = trx(fly).theta(idx(fly));
          
          labelis = find(scores{fly}(:,idx(fly)) > 0);
          tmpscores = scores{fly}(labelis,idx(fly));
          s = cell(1,numel(labelis));
          for iii = 1:numel(labelis),
            s{iii} = sprintf('\\color[rgb]{%f %f %f}%s',tmpscores(iii)*colors(labelis(iii),:),behaviors{labelis(iii)});
          end
          if printnone && isempty(s),
            s = {'None'};
          end
          
%           if doshowsex,
%             s = sprintf('%d, %s',fly,trx(fly).sex{idx(fly)});
%           else
%             s = sprintf('%d',fly);
%           end
          if doplotwings,
            xwingl = trx(fly).xwingl(idx(fly)) - round(trx(fly).x(idx(fly))) + boxradius + .5;
            ywingl = trx(fly).ywingl(idx(fly)) - round(trx(fly).y(idx(fly))) + boxradius + .5;
            xwingl = xwingl * scalefactor;
            ywingl = ywingl * scalefactor;
            xwingl = xwingl + x0(j) - 1;
            ywingl = ywingl + y0(i) - 1;
            xwingr = trx(fly).xwingr(idx(fly)) - round(trx(fly).x(idx(fly))) + boxradius + .5;
            ywingr = trx(fly).ywingr(idx(fly)) - round(trx(fly).y(idx(fly))) + boxradius + .5;
            xwingr = xwingr * scalefactor;
            ywingr = ywingr * scalefactor;
            xwingr = xwingr + x0(j) - 1;
            ywingr = ywingr + y0(i) - 1;
            xwing = [xwingl,x,xwingr];
            ywing = [ywingl,y,ywingr];
          end

          updatefly(hzoom(i,j),x,y,theta,a,b);
          if doplotwings,
            set(hzoomwing(i,j),'XData',xwing,'YData',ywing,'Color',colortmp(fly,:));
          end
          set(htextzoom(i,j),'string',s,'color',colortmp(fly,:));
          if isbehavior(fly),
            set(hzoombox(i,j),'Visible','on','color',colortmp(fly,:));
            set(hzoom(i,j),'Linewidth',2);
          else
            set(hzoombox(i,j),'Visible','off');
            set(hzoom(i,j),'Linewidth',1);
          end
            
%           if doshowsex,
%             sexi = find(strcmpi(trx(fly).sex{idx(fly)},sexes),1);
%             sexmarker = sexmarkers{sexi};
%             set(htextzoom(i,j),'xdata',x,...
%               'ydata',y,...
%               'color',colortmp(fly,:),...
%               'marker',sexmarker,...
%               'markerfacecolor',colortmp(fly,:));
%           end
          set(hzoom(i,j),'color',colortmp(fly,:));
        else
          set(hzoom(i,j),'xdata',[],'ydata',[]);
          if doplotwings,
            set(hzoomwing(i,j),'XData',[],'YData',[]);
          end
          set(htextzoom(i,j),'string','');
%           if doshowsex,
%             set(htextzoom(i,j),'xdata',[],'ydata',[]);
%           end
        end
      end
    end
    
    if frame == firstframes(1),
      if ~isempty(figpos),
        set(1,'Position',figpos);
      else
        input('Resize figure 1 to the desired size, hit enter when done.');
        figpos = get(1,'Position');
      end
      %set(1,'visible','off');
      if ~debug,
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
      end
    end
    
    if ~debug,
    if frame == firstframes(1),
      fr = getframe_invisible(hax);
      [height,width,~] = size(fr);
%       height = ceil(height/4)*4;
%       width = ceil(width/4)*4;
%       fr = getframe_invisible(hax,[height,width]);
    else
      fr = getframe_invisible(hax,[height,width]);
    end
    if useVideoWriter,
      writeVideo(aviobj,fr);
    else
      aviobj = addframe(aviobj,fr);
    end
    set(1,'Position',figpos);
    else
      drawnow;
      if any(isbehavior(:)),
        input('');
      end
    end
    
  end
  
end
  
%% clean up

if useVideoWriter,
  close(aviobj);
else
  aviobj = close(aviobj); %#ok<NASGU>
end
if fid > 0,
  fclose(fid);
end

succeeded = true;