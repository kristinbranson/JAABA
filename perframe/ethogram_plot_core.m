function [boutmat,line_names] = ethogram_plot_core(hfig,labels,T0,T1,expdirs,behaviornames,jabfiles,varargin)
% [boutmat,line_names] = ethogram_plot_core(hfig,labels,expdirs,T0,T1,behaviornames,jabfiles,p1,v1,...)
%
% labels: nexps x nsamps x nbehs label array. 0s and 1s are plotted as
% none and behavior-present, respectively. -1's may also be supplied, which
% are plotted as negative-behavior (faded color).
%
% Optional PVs:
%
% Input arg jabfiles are used only to enable click-on-bout-to-open-in-JAABA.

[trialdividers,pulldownselect,behaviorcolorslines] = ...
  myparse(varargin,...
  'trialdividers',false,...
  'pulldownselect',false,...
  'behaviorcolorslines',false);
  
figure(hfig);
clf;
ax = axes;

nexps = numel(expdirs);
nT = T1-T0+1;
nbehaviors = numel(behaviornames);
assert(size(labels,1)==nexps);
assert(size(labels,2)==nT);
assert(size(labels,3)==nbehaviors);
assert(numel(jabfiles)==nbehaviors,...
  'Input arg ''jabfiles'' must be a cellstr specifying the containing jab for each behavior');

if behaviorcolorslines || nbehaviors>7
  behaviorcolors = lines(nbehaviors);
else
  BEH_COLORS = [
         0         0    1.0000; ...
         0    0.5000         0; ...
    1.0000         0         0; ...
         0    0.7500    0.7500; ...
    0.7500         0    0.7500; ...
    0.7500    0.6000         0; ...
    0.2500    0.2500    0.2500];
  behaviorcolors = BEH_COLORS;
end
nobehaviorcolors = behaviorcolors + 0.67*(1-behaviorcolors); % gray transform
yrad = 1/nbehaviors;
xl = [T0-10,T1+10];

h = cell(nexps,nbehaviors); % patch handles; we keep track of these but ends up unused
hNo = cell(nexps,nbehaviors);
hlegend = nan(1,nbehaviors);
boutmat = cell(nexps,nbehaviors); 
isfirst = true;
for j = 1:nbehaviors,
  behclr = behaviorcolors(j,:);
  nobehclr = nobehaviorcolors(j,:);
  for i = 1:nexps,
    
    if trialdividers && j==1
      % plot dividers for each exp
      plot(xl,repmat(i-yrad/2,1,2),'k:');
    end
    if isfirst,
      isfirst = false;
      hold on;
    end
    
    [starts,ends] = get_interval_ends(labels(i,:,j)>0);
    boutmat{i,j} = [starts(:) ends(:)];
       
    % pos behaviors
    if numel(starts)>0
      h{i,j} = hlpPatch(starts,ends,i,j,yrad,T0,behclr,{'FaceAlpha',0.4});
      if isnan(hlegend(j)),
        hlegend(j) = h{i,j}(1);
      end
    end
    
    % negative behaviors
    [starts,ends] = get_interval_ends(labels(i,:,j)<0);
    if numel(starts)>0
      hNo{i,j} = hlpPatch(starts,ends,i,j,yrad,T0,nobehclr,{'FaceAlpha',0.4});
    end
    
  end
end

line_names = cell(1,nexps);
for i = 1:nexps,
  [~,line_names{i}] = myfileparts(expdirs{i});
end
legend(hlegend(~isnan(hlegend)),behaviornames(~isnan(hlegend)));
yoffset = (nbehaviors-1)/nbehaviors/2;
set(ax,'YTick',yoffset+(1:nexps),'YTickLabel',line_names,'XLim',xl,'YLim',[1,nexps+1],'YDir','reverse','Box','off','TickDir','out');
if isprop(ax,'TickLabelInterpreter')
  ax.TickLabelInterpreter = 'none';
end

% set up axis guidata; data used in callbacks
gdata = struct();
gdata.T0 = T0;
gdata.T1 = T1;
gdata.labels = labels;
gdata.jabfiles = jabfiles;
gdata.expdirs = expdirs;
gdata.Nexp = numel(expdirs);
gdata.behaviornames = behaviornames;
gdata.Nbeh = numel(behaviornames);
gdata.behaviorcolors = behaviorcolors;
gdata.hAutoMarks = [];

% gdata.EXP_SELECTED_COLOR = [1 1 0];
% gdata.EXP_SELECTED_ALPHA = 1;
gdata.EXPZOOMPLUSMINUS = 1; % when zooming on exp, include +/- this many neighboring exps
% gdata.tfExpSelInProgress = false;
% gdata.tfExpIsUserSelected = false(numel(expdirs),1);
% gdata.hUserSelectionPatches = nan(numel(expdirs),1);
% gdata.customGroupVizInfo = cell(0,1); % customGroupVizInfo{i}.name % name of grp
% customGroupVizInfo{i}.tf  % indicator vec for ith viz group
% customGroupVizInfo{i}.h   % plot handles for ith viz group
% Note: due to arbitrary removal, customGroupVizInfo{i} may be [] for
% arbitrary i.
% gdata.CUSTOMGROUPMODN = 7; % this the "mod" for lines()

% gdata.axGroupStats = []; % second/overlay axis for group stats
guidata(ax,gdata);

if pulldownselect
  set(hfig,'units','pixels');
  set(hfig,'position',[100 100 1000 800],'toolbar','figure'); % adding uicontrols to figure below
  
  tmpunits = get(ax,'Units');
  set(ax,'Units','normalized');
  axpos = get(ax,'Position');
  set(ax,'Units',tmpunits);
  
  selections = [{'<all exps>'};expdirs(:)];
  uicontrol('style','popupmenu','String',selections,'fontsize',7,'units','normalized',...
    'position',[axpos(1) 0.955 axpos(3) 0.02],'callback',@(zS,zE)cbkExpZoom(zS,ax));
  uicontrol('style','text','String','Experiment zoom: ','fontweight','bold','horizontalalignment','right',...
    'units','normalized','position',[axpos(1)-0.2 0.95 0.2 0.02],'backgroundcolor',get(hfig,'color'));
end
  
set(ax,'ButtonDownFcn',@(src,evt)cbkAxBDF(src));

function hPatch = hlpPatch(starts,ends,i,j,yrad,T0,clr,patchArgs)
nbouts = numel(starts);
ends = ends-0.5;
starts = starts-0.5;
x = [starts;starts;ends;ends;starts]+T0-1;
y = i+yrad*[0;1;1;0;0]+yrad*(j-1);
xx = x;
yy = repmat(y,[1 nbouts]);
hPatch = patch(xx,yy,-1*ones(size(xx)),clr,...
  'MarkerFaceColor',clr,'HitTest','off','EdgeColor','none',patchArgs{:});


function cbkExpZoom(src,ax)
gdata = guidata(ax);

val = get(src,'value');
str = get(src,'String');
str = str{val};

expnum = find(strcmp(str,gdata.expdirs));
if isempty(expnum)
  ylim(ax,'auto');
else
  ylim(ax,[expnum-gdata.EXPZOOMPLUSMINUS expnum+gdata.EXPZOOMPLUSMINUS+1]);
end

function cbkAxBDF(ax)

gdata = guidata(ax);
ty = get(ancestor(ax,'figure'),'SelectionType');
switch ty
  case 'extend' % Shift-click
    yl = ylim(ax);
    if diff(yl)<=1.5
      % fully zoomed; zoom out
      ylim(ax,'auto');
    elseif diff(yl)<=gdata.EXPZOOMPLUSMINUS*2+1
      % medium zoom; go to full zoom (single exp)
      pt = get(ax,'currentpoint');
      expnum = floor(pt(1,2));
      ylim(ax,[expnum expnum+1]);
    else
      % go to medium zoom
      pt = get(ax,'currentpoint');
      expnum = floor(pt(1,2));
      ylim(ax,[expnum-gdata.EXPZOOMPLUSMINUS expnum+gdata.EXPZOOMPLUSMINUS+1]);
    end
  case 'alt' % Ctrl-click
%     pt = get(ax,'currentpoint');
%     expnum = floor(pt(1,2));
%     gdata = guidata(ax);
%     if expnum > gdata.Nexp % Can occur if ctrl-click slightly 'below' axis
%       warning('ethogram_plot:noSuchExp','Clicked point beyond range of experiments.');
%     else
%       gdata = lclToggleExpSelection(gdata,expnum);
%     end
%     gdata.tfExpSelInProgress = true;
%     guidata(ax,gdata);
  otherwise
    patchClickCallback(ax);
end

function patchClickCallback(h)
% Check for double click

% AL 20150113: This is not working

persistent chk

dothis = true;
if isempty(chk)
  chk = 1;
  pause(0.5); %Add a delay to distinguish single click from a double click
  if chk == 1
    chk = [];
    return;
  end
else
  dothis = false;
  chk = [];
end

if ~dothis, return; end
% Below executes only on double click

gdata = guidata(h);

% if gdata.scoresnotjabs
%   warning('ethogram_plot_core:nojabs','No jab files specified; clicking on bouts is disabled.');
%   return;
% end

pt = get(h,'currentpoint');
frameno = round(pt(1,1));
if frameno<0 || frameno>size(gdata.labels,2),
  return;
end
expnum = floor(pt(1,2));

% determine behavior number/index so we know which jab to open
behnum = find(gdata.labels(expnum,frameno,:));
if numel(behnum)==0, behnum = 1:numel(gdata.behaviornames); end
if numel(behnum)>1
  [sel,ok] = listdlg('ListSize',[650 220],'Name','Select a Behavior to see in JAABA',...
    'ListString',gdata.behaviornames(behnum),'SelectionMode','single','OKString','Select',...
    'InitialValue',1);
  if ok == 0, return; end
  behnum = behnum(sel);
end
assert(any(behnum==(1:gdata.Nbeh)));
selectedJab = gdata.jabfiles{behnum};

JLabelInterf.openJabExpFrame(selectedJab,gdata.expdirs{expnum},frameno);
