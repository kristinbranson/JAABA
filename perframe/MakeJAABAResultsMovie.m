% MakeJAABAResultsMovie(expdir,classifierparamsfile,'nframesperseg',500,'nintervals',10,'aviname','JAABAResultsMovie_FCF_cantons_1500002_None_Rig1Plate15BowlA_20120519T160453.avi')

function MakeJAABAResultsMovie(expdir,classifierparamsfile,varargin)

%% parse parameters

frac_bottom = .1;
border_bottom = .05;
border_middle = .02;
trx_linewidth = 1;
target_linewidth = 2;
bar_width = .8;
trx_rad = 25;
pxborder = [];
aborder = 1;
trxcolor = [0,0,0];
figpos = [50,50,1500,1000];
barbackcolor = [.5,.5,.5];
textcolor = [0,0,0];
fontsize = 30;
aback = .25;
printnone = true;
mintrajlength = 100;
nintervals = 1;
compression = 'None';
allowedcompressions = {'Indeo3', 'Indeo5', 'Cinepak', 'MSVC', 'RLE', 'None','Uncompressed AVI','Motion JPEG AVI'};
useVideoWriter = exist('VideoWriter','file');
aviname = '';
fps = 10;
usemencoder = true;

[framestarts,nframesperseg,fracstarts,fracperseg,targets,hfig,...
  pxwidthradius,pxheightradius,awidthradius,aheightradius,frac_bottom,...
  trx_linewidth,target_linewidth,...
  behavior_colors,behavior_cm,...
  border_bottom,border_middle,bar_width,...
  trx_rad,score_norm_prctile,...
  pxborder,aborder,...
  trxcolor,...
  figpos,...
  barbackcolor,...
  behaviorweight,...
  textcolor,fontsize,...
  pxback,aback,...
  printnone,...
  mintrajlength,...
  nintervals,...
  useVideoWriter,...
  avifileTempDataFile,aviname,fps,...
  usemencoder,titlepagetext,...
  titlepausetime] = ...
  myparse(varargin,'framestarts',[],'nframesperseg',[],...
  'fracstarts',[],'fracperseg',1/60,...
  'targets',[],...
  'hfig',1,...
  'pxwidthradius',[],...
  'pxheightradius',[],...
  'awidthradius',3,...
  'aheightradius',3,...
  'frac_bottom',frac_bottom,...
  'trx_linewidth',trx_linewidth,...
  'target_linewidth',target_linewidth,...
  'behavior_colors',zeros(0,3),...
  'behavior_cm',@(x) jet(x)*.7,...
  'border_bottom',border_bottom,...
  'border_middle',border_middle,...
  'bar_width',bar_width,...
  'trx_rad',trx_rad,...
  'score_norm_prctile',80,...
  'pxborder',pxborder,...
  'aborder',aborder,...
  'trxcolor',trxcolor,...
  'figpos',figpos,...
  'barbackcolor',barbackcolor,...
  'behaviorweight',[],...
  'textcolor',textcolor,...
  'fontsize',fontsize,...
  'pxback',[],...
  'aback',aback,...
  'printnone',printnone,...
  'mintrajlength',mintrajlength,...
  'nintervals',nintervals,...
  'useVideoWriter',useVideoWriter,...
  'avifileTempDataFile','',...
  'aviname',aviname,...
  'fps',fps,...
  'usemencoder',usemencoder,...
  'titlepagetext',{},...
  'titlepausetime',2);

classifierparams = ReadClassifierParamsFile(classifierparamsfile);
nbehaviors = numel(classifierparams);
behaviors = cell(1,nbehaviors);
scorefns = cell(1,nbehaviors);
for i = 1:nbehaviors,
  if iscell(classifierparams(i).behaviors.names),
    behaviors{i} = sprintf('%s_',classifierparams(i).behaviors.names{:});
    behaviors{i} = behaviors{i}(1:end-1);
  else
    behaviors{i} = classifierparams(i).behaviors.names;
  end
  behaviors{i} = regexprep(behaviors{i},'(_|^)([a-z])','${upper($2)}');
  scorefn = classifierparams(i).file.scorefilename;
  scorefns{i} = regexprep(scorefn,'\.mat$','');
end

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
  aviname = fullfile(expdir,sprintf('JAABAResultsMovie_%s.avi',experiment_name));
end

%% create trx structure

trx = Trx('trxfilestr',classifierparams(1).file.trxfilename,...
  'moviefilestr',classifierparams(1).file.moviefilename,...
  'perframedir',classifierparams(1).file.perframedir);
if isfield(classifierparams,'perframe') && isfield(classifierparams(1).perframe,'params'),
  trx.SetPerFrameParams(classifierparams(1).perframe.params);
end  
if isfield(classifierparams,'perframe') && isfield(classifierparams(1).perframe,'landmark_params'),
  trx.SetLandmarkParams(classifierparams(1).perframe.landmark_params);
end  

trx.AddExpDir(expdir,'openmovie',false);

%% for video reading

moviename = fullfile(expdir,trx.moviefilestr);
[readframe,~,fid,headerinfo] = get_readframe_fcn(moviename);

%% choose frames, target to show

% how many intervals
nintervals = max([nintervals,numel(targets),numel(nframesperseg),...
  numel(framestarts),numel(fracperseg),numel(fracstarts)]);

% choose frame starts evenly spaced throughout the movie
if isempty(framestarts) && isempty(fracstarts),
  framestarts = round(linspace(max(trx.firstframe),min(trx.endframe),nintervals));
end

if ~isempty(nframesperseg) && numel(nframesperseg) < nintervals,
  nframesperseg = [nframesperseg,nframesperseg(end)+zeros(1,nintervals-numel(nframesperseg))];
end

% make sizes match
if isempty(targets),
  for segi = 1:nintervals,
    if numel(nframesperseg) >= segi,
      mintrajlength_curr = max([mintrajlength,nframesperseg(segi)]);
      idx = trx.nframes>mintrajlength_curr;
      if numel(framestarts) >= segi,
        idx = idx & trx.firstframe<=framestarts(segi) & ...
          trx.endframe >= framestarts(segi);
      end
      idx = setdiff(find(idx),targets(1:segi-1));
      targets(segi) = randsample(idx,1);
    end
  end
elseif numel(targets) < nintervals,
  targets = [targets,targets(end)+zeros(1,nintervals-numel(targets))];
end
if numel(fracperseg) < nintervals,
  fracperseg = [fracperseg,fracperseg(end)+zeros(1,nintervals-numel(fracperseg))];
end
if numel(fracstarts) < nintervals && ~isempty(fracstarts),
  fracstarts = [fracstarts,fracstarts(end)+zeros(1,nintervals-numel(fracstarts))];
end

frameends = nan(1,nintervals);
for targeti = 1:nintervals,
  target = targets(targeti);
  firstframe = trx.firstframes(target);
  endframe = trx.endframes(target);
  nframes_curr = endframe - firstframe + 1;

  if numel(nframesperseg) < targeti || isnan(nframesperseg(targeti)),
    nframesperseg(targeti) = round(nframes_curr*fracperseg(targeti));
  end
  if numel(framestarts) < targeti || isnan(framestarts(targeti)),
    framestarts(targeti) = firstframe + round(fracstarts(targeti)*nframes_curr);
  end
  frameends(targeti) = max(firstframe,min(endframe,framestarts(targeti)+nframesperseg(targeti)-1));
end
%nframesperseg = frameends-framestarts+1;

%% get the score normalizations

score_norms = nan(nbehaviors,1);
for i = 1:nbehaviors,
  tmp = trx(targets).(scorefns{i});
  if iscell(tmp),
    tmp = [tmp{:}];
  end
  score_norms(i) = prctile(abs(tmp),score_norm_prctile);
end

%% choose weights for each behavior

if isempty(behaviorweight),
  behaviorweight = nan(nbehaviors,1);
  for behaviori = 1:nbehaviors,
    tmp = trx(targets).(scorefns{behaviori});
    if iscell(tmp),
      tmp = [tmp{:}];
    end
    behaviorweight(behaviori) = nnz(tmp>0);
  end
  behaviorweight = behaviorweight / sum(behaviorweight);
end
if numel(behaviorweight) == 1,
  behaviorweight = ones(nbehaviors,1);
end
behaviorweight = behaviorweight(:);

%% choose colors for each behavior

nadd = nbehaviors - size(behavior_colors,1);
if nadd > 0,
  tmp = behavior_cm(256);
  behavior_colors = cat(1,behavior_colors,tmp(round(linspace(1,256,nadd)),:));
end

%% initialize plot

if ishandle(hfig),
  set(0,'CurrentFigure',hfig);
  clf(hfig);
else
  hfig = figure(hfig);
end
set(hfig,'Units','pixels','Position',figpos);

mainpos = [0,frac_bottom+border_bottom+border_middle,1,1-frac_bottom-border_bottom-border_middle];
haxmain = axes('Position',mainpos);
haxbottom = axes('Position',[0,border_bottom,1,frac_bottom]);

t = framestarts(1);
im = readframe(t);
ncolors = size(im,3);
if ncolors == 1,
  him = imagesc(im2double(im),'Parent',haxmain,[0,1]);
  colormap(haxmain,'gray');
else
  him = image(im2double(im),'Parent',haxmain);
end
hold(haxmain,'on');
%axis(haxmain,'image');
htrx = plot(haxmain,nan,nan,'.-','LineWidth',trx_linewidth,'color',trxcolor);
htarget = plot(haxmain,nan,nan,'-','LineWidth',target_linewidth,'color','k');
set(haxmain,'XTick',[],'YTick',[]);
i = t+trx(target).off;

a = trx(targets).a;
if iscell(a),
  a = [a{:}];
end
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
if isempty(pxback),
  pxback = meana*aback;
end

set(haxmain,'Units','pixels');
mainpos = get(haxmain,'Position');
set(haxmain,'Units','normalized');
pxheightradius1 = min((pxwidthradius)*mainpos(4) / mainpos(3),(headerinfo.nr-1)/2);
pxwidthradius1 = min((pxheightradius)*mainpos(3) / mainpos(4),(headerinfo.nc-1)/2);
pxheightradius = max(pxheightradius1,pxheightradius);
pxwidthradius = max(pxwidthradius1,pxwidthradius);
ax = ResetAxis(headerinfo.nc/2,headerinfo.nr/2,1,headerinfo.nr,headerinfo.nc,pxwidthradius,pxheightradius,pxborder);

hinfo = text(ax(1),ax(4),'','Parent',haxmain,'Color',textcolor,...
  'FontUnits','pixels','FontSize',fontsize,...
  'HorizontalAlignment','left',...
  'VerticalAlignment','bottom');

axis(haxmain,'image');
axis(haxmain,ax);

hbar = nan(1,nbehaviors);
hbarback = nan(1,nbehaviors);
for i = 1:nbehaviors,
  hbarback(i) = patch(i+bar_width/2*[-1,1,1,-1,-1],zeros(1,5)+.1,barbackcolor,'EdgeColor','none','Parent',haxbottom);
  if i == 1,
    hold(haxbottom,'on');
  end
  hbar(i) = patch(i+bar_width/2*[-1,1,1,-1,-1],zeros(1,5)+.1,behavior_colors(i,:),'EdgeColor','none','Parent',haxbottom);
end
htext = text(nan,nan,'','Parent',haxmain,'Color',textcolor,...
  'FontUnits','pixels','FontSize',fontsize);
plot(haxbottom,[0,nbehaviors+1],[0,0],'--','color',[.5,.5,.5]);
set(haxbottom,'XTick',1:nbehaviors,'XTickLabel',behaviors,'YTick',[],...
  'XLim',[0,nbehaviors+1],'YLim',[-1,1],'XColor','w','Color','k','TickLength',[0,0]);

set(hfig,'Color','k');

%% update plots

isfirstframe = true;

if ~isempty(titlepagetext),
  titlepagetext2 = sprintf(', %s',titlepagetext{:});
end


for segi = 1:numel(framestarts),
  
  % data for this target
  target = targets(segi);

  nframes_curr = trx(target).nframes;
  x = trx(target).x;
  y = trx(target).y;
  a = trx(target).a;
  b = trx(target).b;
  sex = trx(target).sex;
  theta = trx(target).theta;
  % try to make sector somewhat continuous
  fil = [0.000263865082737   0.106450771973592   0.786570725887342   0.106450771973592   0.000263865082737];
  fil = fil / sum(fil(:));
  smooththeta = imfilter(theta,fil,'same','replicate');  
  sector = floor(modrange((smooththeta-pi/8),0,2*pi)/pi*4)+1;
  % close holes
  if numel(sector) > 2,
    for i = 1:numel(sector),
      if ((i==1) || (i==numel(sector)) || (sector(i-1)==sector(i+1))) && ...
          (i==1 || sector(i)~=sector(i-1)) && ...
          (i==numel(sector) || sector(i)~=sector(i+1)) && ...
          (i==1 || abs(modrange(theta(i)-theta(i-1),-pi,pi))<pi/12) && ...
          (i==numel(sector) || abs(modrange(theta(i)-theta(i+1),-pi,pi))<pi/12),
        if i == 1,
          sector(i) = sector(i+1);
        else
          sector(i) = sector(i-1);
        end
      end
    end
  end

  scores = nan(nbehaviors,nframes_curr);
  for behaviori = 1:nbehaviors,
    scores(behaviori,:) = min(1,max(-1,trx(target).(scorefns{behaviori})/score_norms(behaviori)));
  end
  
  for t = framestarts(segi):frameends(segi),
    im = readframe(t);
    set(him,'CData',im2double(im));
    i = t+trx(target).off;
    
    i0 = max(1,min(nframes_curr,i-trx_rad(1)));
    i1 = max(1,min(nframes_curr,i+trx_rad(1)));
    
    set(htrx,'XData',x(i0:i1),...
      'YData',y(i0:i1));
    
    updatefly(htarget,x(i),y(i),...
      theta(i),a(i),b(i));
    
    for behaviori = 1:nbehaviors,
      set(hbar(behaviori),'YData',[0,0,1,1,0]*scores(behaviori,i));
      set(hbarback(behaviori),'YData',[0,0,1,1,0]*sign(scores(behaviori,i)));
    end
    w = max(0,scores(:,i)).*behaviorweight;
    z = sum(w);
    if z == 0,
      colortmp = [0,0,0];
    else
      colortmp = sum(bsxfun(@times,w,behavior_colors),1)/z*max(scores(:,i));
    end
    set(htarget,'Color',colortmp);
    
    tailpos = [x(i)-(pxback+2*a(i))*cos(theta(i)),y(i)-(pxback+2*a(i))*sin(theta(i))];      
    sectori = sector(i);
    % floor(modrange((smooththeta(i)-pi/8),0,2*pi)/pi*4)+1;
    
    if ismember(sectori,[1,7,8]),
      halign = 'right';
    elseif ismember(sectori,[2,6]),
      halign = 'center';
    elseif ismember(sectori,[3,4,5]),
      halign = 'left';
    else
      error('Sanity check: sectori > 8 or < 1');
    end
    if ismember(sectori,[1,2,3]),
      valign = 'bottom';
    elseif ismember(sectori,[4,8]),
      valign = 'middle';
    elseif ismember(sectori,[5,6,7]),
      valign = 'top';
    else
      error('Sanity check: sectori > 8 or < 1');
    end
    
    labelis = find(scores(:,i) > 0);
    s = cell(1,numel(labelis));
    for iii = 1:numel(labelis),
      s{iii} = sprintf('\\color[rgb]{%f %f %f}%s',behavior_colors(labelis(iii),:),behaviors{labelis(iii)});
    end
    if printnone && isempty(s),
      s = {'None'};
    end
    
    set(htext,'String',s,...
      'Position',tailpos,...
      'HorizontalAlignment',halign,...
      'VerticalAlignment',valign);

    ss = sprintf('Target %d (%s), Frame %d',target,sex{i},t);
    if ~isempty(titlepagetext),
      ss = [ss,titlepagetext2];
    end

    set(hinfo,'String',ss)
    
    outofbounds = x(i)-pxborder < ax(1) || ...
      x(i)+pxborder > ax(2) || ...
      y(i)-pxborder < ax(3) || ...
      y(i)+pxborder > ax(4);
    if t == framestarts(segi) || outofbounds,
      ax = ResetAxis(x,y,i,headerinfo.nr,headerinfo.nc,pxwidthradius,pxheightradius,pxborder);
      axis(haxmain,ax);
      set(hinfo,'Position',ax([1,4]));
    end
    drawnow;
    
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
      
      if ~isempty(titlepagetext),
        set(htext,'Visible','off');
        htitle = text(mean(ax(1:2)),mean(ax(3:4)),...
          titlepagetext,'HorizontalAlignment','center',...
          'VerticalAlignment','middle',...
          'Color',[.5,0,0],'Parent',haxmain,...
          'FontSize',40,'FontWeight','bold');
        fr = getframe(hfig);
        for pausei = 1:fps*titlepausetime,
          if useVideoWriter,
            writeVideo(aviobj,fr);
          else
            aviobj = addframe(aviobj,fr);
          end
        end
        set(htext,'Visible','on');
        delete(htitle);
      end
      
      isfirstframe = false;
    end

    % pause for a few frames at the start of a new segment
    if t == framestarts(segi) && segi > 1,
      set(hinfo,'Color',[.5,0,0],'FontWeight','bold');
      fr = getframe(hfig);
      for pausei = 1:fps,
        if useVideoWriter,
          writeVideo(aviobj,fr);
        else
          aviobj = addframe(aviobj,fr);
        end
      end
      set(hinfo,'Color',textcolor,'FontWeight','normal');
    end
    
    fr = getframe(hfig);
    height = size(fr.cdata,1);
    width = size(fr.cdata,2);
    if useVideoWriter,
      writeVideo(aviobj,fr);
    else
      aviobj = addframe(aviobj,fr);
    end
    
    
    set(hfig,'Position',figpos);
    
  end
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


function ax = ResetAxis(x,y,i,nr,nc,pxwidthradius,pxheightradius,pxborder)

n = numel(x);

minx = x(i);
maxx = x(i);
miny = y(i);
maxy = y(i);
w = pxwidthradius-2*pxborder;
h = pxheightradius-2*pxborder;
for j = i+1:n,
  minx = min(minx,x(j));
  maxx = max(maxx,x(j));
  dx = ceil(maxx-minx)+1;
  if dx > w,
    break;
  end
  miny = min(miny,y(j));
  maxy = max(maxy,y(j));
  dy = ceil(maxy-miny)+1;
  if dy > h,
    break;
  end
end
mux = (minx+maxx)/2;
muy = (miny+maxy)/2;
ax = [mux+pxwidthradius*[-1,1],muy+pxheightradius*[-1,1]];

if ax(2)>nc && ax(1)<1,
  ax(1:2) = (nc+1)/2+pxwidthradius*[-1,1];
elseif ax(2)>nc,
  ax(2) = nc;
  ax(1) = nc-(2*pxwidthradius+1)+1;
elseif ax(1)<1,
  ax(1) = 1;
  ax(2) = 2*pxwidthradius+1;
end
if ax(2)>nc || ax(1)<1,
  ax(1:2) = (nc+1)/2+pxwidthradius*[-1,1];
end
if ax(4)>nr && ax(3)<1,
  ax(3:4) = (nr+1)/2+pxwidthradius*[-1,1];
elseif ax(4)>nr,
  ax(4) = nr;
  ax(3) = nr-(2*pxheightradius+1)+1;
elseif ax(3)<1,
  ax(3) = 1;
  ax(4) = 2*pxheightradius+1;
end
if ax(4)>nr || ax(3)<1,
  ax(3:4) = (nr+1)/2+pxheightradius*[-1,1];
end

