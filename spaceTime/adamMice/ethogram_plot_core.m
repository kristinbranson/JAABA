function [boutmat,line_names] = ethogram_plot_core(hfig,labels,T0,T1,expdirs,behaviornames,jabfiles,varargin)
% [boutmat,line_names] = ethogram_plot_core(hfig,labels,expdirs,T0,T1,behaviornames,jabfiles,p1,v1,...)
%
% labels: nexps x nsamps x nbehs label array. 0s and 1s are plotted as
% none and behavior-present, respectively. -1's may also be supplied, which
% are plotted as negative-behavior (faded color).
%
% boutmat: nexps x nbehaviors x 2 cell array. Rows/Cols of boutmat are
% labeled by line_names/behaviornames resp. boutmat{:}{:}{1} contain
% positive-behaviors, boutmat{:}{:}{2} contain no-behaviors.
% line_names: vector cellstr corresponding to expdirs.
%
% Optional PVs:
% - doautomarks. scalar logical.
% - automarkdata. Used if doautomarks is true. Data struct from ExpPP.
% - scoresnotjabs.
%
% Input arg jabfiles are used only to enable
% click-on-bout-to-open-in-JAABA.

[doautomarks,automarkdata,expppstatprefix,expppcontrolargs,scoresnotjabs,trialdividers,pulldownselect,behaviorcolorslines] = ...
  myparse(varargin,...
  'doautomarks',false,...
  'automarkdata',[],...
  'expppstatprefix','auto',...
  'expppcontrolargs',{},...
  'scoresnotjabs',false,...
  'trialdividers',false,...
  'pulldownselect',true,...
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
boutmat = cell(nexps,nbehaviors,2); % boutmat{iExp}{iBeh}{1=positive beh, 2=no-beh}
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
    boutmat{i,j,1} = [starts(:) ends(:)];
    if numel(starts)>0
      h{i,j} = hlpPatch(starts,ends,i,j,yrad,T0,behclr,doautomarks,{'FaceAlpha',0.4});
      if isnan(hlegend(j)),
        hlegend(j) = h{i,j}(1);
      end
    end
    
    % negative behaviors
    [starts,ends] = get_interval_ends(labels(i,:,j)<0);
    boutmat{i,j,2} = [starts(:) ends(:)];
    if numel(starts)>0
      hNo{i,j} = hlpPatch(starts,ends,i,j,yrad,T0,nobehclr,doautomarks,{'FaceAlpha',0.4});
    end
  end
end

line_names = cell(1,nexps);
for i = 1:nexps,
  [~,fname,ext] = fileparts(expdirs{i});
  line_names{i} = [fname ext];
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
gdata.scoresnotjabs = scoresnotjabs;
gdata.expdirs = expdirs;
gdata.Nexp = numel(expdirs);
gdata.behaviornames = behaviornames;
gdata.Nbeh = numel(behaviornames);
gdata.behaviorcolors = behaviorcolors;
gdata.hAutoMarks = [];
gdata.automarkFldPfix = expppstatprefix;

gdata.EXP_SELECTED_COLOR = [1 1 0];
gdata.EXP_SELECTED_ALPHA = 1;
gdata.EXPZOOMPLUSMINUS = 1; % when zooming on exp, include +/- this many neighboring exps
gdata.tfExpSelInProgress = false;
gdata.tfExpIsUserSelected = false(numel(expdirs),1);
gdata.hUserSelectionPatches = nan(numel(expdirs),1);
gdata.customGroupVizInfo = cell(0,1); % customGroupVizInfo{i}.name % name of grp
% customGroupVizInfo{i}.tf  % indicator vec for ith viz group
% customGroupVizInfo{i}.h   % plot handles for ith viz group
% Note: due to arbitrary removal, customGroupVizInfo{i} may be [] for
% arbitrary i.
gdata.CUSTOMGROUPMODN = 7; % this the "mod" for lines()

gdata.axGroupStats = []; % second/overlay axis for group stats
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
  
if doautomarks  
  initSelectedExpPatches(ax);
  
  set(ax,'ButtonDownFcn',@(src,evt)cbkAxBDF(ax));
  set(hfig,'ButtonDownFcn',@(src,evt)cbkFigBDF(ax));
  set(hfig,'WindowButtonMotionFcn',@(src,evt)cbkFigWBMF(ax));
  set(hfig,'WindowButtonUpFcn',@(src,evt)cbkFigWBUF(ax));
%   set(hfig,'ResizeFcn',@(src,evt)cbkFigResized(ax));
  
  callbacks.drawAutoMarks = @(zData)cbkDrawAutoMarks(ax,zData);
  callbacks.drawAutoMarkGroupStats = @(zStats,zGrpFlds)cbkDrawGroupStats(ax,zStats,zGrpFlds);
  callbacks.getSelectedExps = @()cbkGetSelectedExps(ax);
  callbacks.addCustomGroupViz = @(zGrpname,zTf)cbkAddCustomGroupViz(ax,zGrpname,zTf);
  callbacks.rmCustomGroupViz = @(grpname)cbkRmCustomGroupViz(ax,grpname);
  callbacks.setSelectedExps = @(tf)cbkSetSelectedExps(ax,tf);
  callbacks.exportFig = @()cbkExportFig(ax);
  
  assert(numel(automarkdata)==gdata.Nexp);

  ExpPPControl(automarkdata,callbacks,{'computeonlabels' true 'statframelimits' [-inf T1]},hfig,expppcontrolargs{:});
else
  set(ax,'ButtonDownFcn',@(src,evt)patchClickCallback(src));
end

function hPatch = hlpPatch(starts,ends,i,j,yrad,T0,clr,doautomarks,patchArgs)
nbouts = numel(starts);
ends = ends-0.5;
starts = starts-0.5;
x = [starts;starts;ends;ends;starts]+T0-1;
if doautomarks
  y = i+yrad*[1;0.5;0.5;1;1]+yrad*(j-1);
else
  y = i+yrad*[0;1;1;0;0]+yrad*(j-1);
end
xx = x;
yy = repmat(y,[1 nbouts]);
hPatch = patch(xx,yy,-1*ones(size(xx)),clr,...
  'MarkerFaceColor',clr,'HitTest','off','EdgeColor','none',patchArgs{:});
      
function cbkExportFig(ax)
hFig = ancestor(ax,'figure');
ethogram_export(hFig);

function cbkDrawAutoMarks(ax,dauto)
% cbkDrawAutoMarks(ax,dauto)
% Visualize automarks.
% 
% - dauto is computed externally (eg in ExpPPControl), but
% knowledge of how to visualize appropriately within the ethogramplot is
% put here.
% - Various plot metadata is stored in the axis guidata. Axis guidata is
% also used to store persistent state for automarks (eg patch handles). 
% - ethogram_plot_core may ultimately be more naturally structured as a
% class for state-sharing, depending on ultimate usage patterns.
% - Currently this duplicates a bit of ethogram_plot_core, maybe suboptimal 
% but wanted to (try to) enable cleaner on-the-fly adjustment of 
% postprocessing params.

gdata = guidata(ax);
delete(gdata.hAutoMarks(ishandle(gdata.hAutoMarks)));
gdata.hAutoMarks = [];

dautoexpfulls = {dauto.expfull}';
nexps = numel(gdata.expdirs);
nbehaviors = numel(gdata.behaviornames);
yrad = 1/nbehaviors;
for j = 1:nbehaviors
  currbeh = gdata.behaviornames{j};
  tfbeh = regexpmatch(currbeh,ExpPP.BASICBEHAVIORS,'caseinsens',true);
  if any(tfbeh)
    assert(nnz(tfbeh)==1);
    currbeh = ExpPP.BASICBEHAVIORS{tfbeh};
    
    % Get fields to plot for this behavior
    fldt0 = sprintf('%s_%s_t0s',gdata.automarkFldPfix,currbeh);
    fldt1 = sprintf('%s_%s_t1s',gdata.automarkFldPfix,currbeh);
    assert(all(isfield(dauto,{fldt0;fldt1})),'Data structure missing field ''%s'' or ''%s''.',fldt0,fldt1);
    fldgs00 = sprintf('%s_GS00_%s_0',gdata.automarkFldPfix,currbeh);
    fldgs01 = sprintf('%s_GS01_%s_0',gdata.automarkFldPfix,currbeh);
    fldgsss = sprintf('%s_GSSS_%s_0',gdata.automarkFldPfix,currbeh);
    fldSuccess = sprintf('%s_Grab_success',gdata.automarkFldPfix);
    fldSuccessType = sprintf('%s_Grab_successtype',gdata.automarkFldPfix);
    
    hBoutPatches = nan(nexps,1);
    hGS00Marks = nan(nexps,1);
    hGS01Marks = nan(nexps,1);
    hGSSSMarks = nan(nexps,1);
    hSuccessPatches = nan(nexps,1);    
    for i = 1:nexps
      currexp = gdata.expdirs{i};
      tfexp = strcmp(currexp,dautoexpfulls);
      if any(tfexp)
        assert(nnz(tfexp)==1);
    
        y = i+yrad*(j-1);

        % plot filled in bouts
        starts = dauto(tfexp).(fldt0);
        ends = dauto(tfexp).(fldt1);
        ends = ends-0.5;
        starts = starts-0.5;
        nbouts = numel(starts);
        if nbouts>0
          xx = [starts(:)';starts(:)';ends(:)';ends(:)';starts(:)'];
          yy = repmat(y+yrad*[0;0.5;0.5;0;0],[1 nbouts]);
          hBoutPatches(i) = patch(xx,yy,-1*ones(size(xx)),gdata.behaviorcolors(j,:),...
            'EdgeColor','none','MarkerFaceColor',gdata.behaviorcolors(j,:),'FaceAlpha',0.4,...
            'HitTest','off','parent',ax);
        end        
    
%         % plot t0s 
%         t0vec = dauto(tfexp).(fldt0);
%         if ~isempty(t0vec)
%           tmp = plot(ax,t0vec(:),y,'o','Color',gdata.behaviorcolors(j,:),'HitTest','off');
%           gdata.hAutoMarks = [gdata.hAutoMarks tmp(:)']; 
%         end
%         % plot flds01
%         for f = flds0(:)', f = f{1}; %#ok<FXSET>
%           xx = dauto(tfexp).(f);
%           gdata.hAutoMarks(end+1) = plot(ax,xx,y,'o','Color',gdata.behaviorcolors(j,:),...
%             'MarkerFaceColor',gdata.behaviorcolors(j,:),'HitTest','off'); 
%         end        

        % plot gs00/gs01/gsss 
        hGS00Marks(i) = hlpGSMarks(ax,dauto(tfexp),fldgs00,y+0.125*yrad,'x',gdata.behaviorcolors(j,:),'udGS00','MarkerSize',8);
        hGS01Marks(i) = hlpGSMarks(ax,dauto(tfexp),fldgs01,y+0.25*yrad,'^',gdata.behaviorcolors(j,:),'udGS01','MarkerSize',5,'MarkerEdgeColor','none');
        hGSSSMarks(i) = hlpGSMarks(ax,dauto(tfexp),fldgsss,y+0.375*yrad,'s',gdata.behaviorcolors(j,:),'udGSSS','MarkerSize',5,'MarkerEdgeColor','none');
        
        % shading for grab success-state
        SUCCESS_HIGHLIGHT_ALPHA = 0.1;
        if strcmpi(currbeh,'grab') && dauto(tfexp).(fldSuccess)
          successtype = dauto(tfexp).(fldSuccessType);
          successcolor = ExpPP.TRIALSUCCESSCOLOR.(successtype);
          xx = [gdata.T0 gdata.T0 gdata.T1 gdata.T1];
          yy = [i i+1 i+1 i];
          hSuccessPatches(i) = patch(xx,yy,-2.5*ones(size(xx)),successcolor,...
            'EdgeColor','none','FaceAlpha',SUCCESS_HIGHLIGHT_ALPHA,...
            'HitTest','off','Parent',ax);
        end
      end
    end
    gdata.hAutoMarks = [gdata.hAutoMarks;hBoutPatches(:);hGS00Marks(:);hGS01Marks(:);hGSSSMarks(:);hSuccessPatches(:)];
  end
end

guidata(ax,gdata);

function hLine = hlpGSMarks(ax,dEl,statfld,yy,mrkr,clr,udata,varargin)
hLine = nan;
if isfield(dEl,statfld)
  xx = dEl.(statfld);
  if ~isnan(xx)
    hLine = plot(ax,xx,yy,mrkr,'Color',clr,...
      'MarkerFaceColor',clr,'HitTest','off','userdata',udata,varargin{:});
  end
end

function cbkDrawGroupStats(ax,stats,grpFlds)
% - stats: stats structure
% - grpFlds : cellstr, fields of dauto that specify indicator vecs for
% custom-group-stats plotting. Can be empty in which case no
% custom-group-stats are shown.

gdata = guidata(ax);

% if it doesn't exist, create overlay axis for stats (purpose is to contain
% second legend)
if isempty(gdata.axGroupStats)
  gdata.axGroupStats = copyobj(ax,ancestor(ax,'figure')); % TODO: This only occurs once but takes a loooong time 
  set(gdata.axGroupStats,'Color','none','XTick',[],'YTick',[],'Box','off','userData','udGroupStatsAx','hittest','off');  
  linkaxes([ax gdata.axGroupStats]);
  set(ax,'ActivePositionProperty','position');
  set(gdata.axGroupStats,'ActivePositionProperty','position');
  set(gdata.axGroupStats,'Position',get(ax,'Position'));
end
delete(get(gdata.axGroupStats,'Children')); % delete all lines

% plot stats 
STATFLDS = {'lift' sprintf('%s_Lift_0',gdata.automarkFldPfix); ...   
            'grab' sprintf('%s_Grab_0',gdata.automarkFldPfix)};
STATGRPMARKERS = {'-' ':' '-.' '--'};
assert(iscellstr(grpFlds));
if numel(grpFlds)>4
  warning('ethogram_plot:tooManyGrps','Can only show 4 groups at once.');
  grpFlds = grpFlds(1:4);
end
STATGRPS(:,1) = grpFlds(:);
Ngrps = size(STATGRPS,1);
STATGRPS(:,2) = STATGRPMARKERS(1:Ngrps);

tfLegUpdated = false;
for iBeh = 1:size(STATFLDS,1)
  beh = STATFLDS{iBeh,1};
  tf = regexpmatch(gdata.behaviornames,beh,'caseinsens',true);
  if any(tf)
    assert(nnz(tf)==1,'Multiple ''%s'' behaviors specified.',beh);
    fld = STATFLDS{iBeh,2};
    clr = gdata.behaviorcolors(tf,:);

    hStatLines = zeros(0,1); % stats lines for current beh, all groups
    for iGrp = 1:Ngrps
      grpname = STATGRPS{iGrp,1};
      marker = STATGRPS{iGrp,2};
      statsGrp = stats.(fld).(grpname);
      fprintf(1,'Stats for ''%s'', grouping ''%s''. %d exps. N=%d, mdn=%.3f, sd=%.3f.\n',...
        fld,grpname,gdata.Nexp,statsGrp.N,statsGrp.median,statsGrp.std);
      hStatLines(iGrp,1) = plot(gdata.axGroupStats,[statsGrp.median statsGrp.median],[0 gdata.Nexp+1],...
        marker,'Color',clr,'HitTest','off'); 
    end
    
    set(gdata.axGroupStats,'Color','None');
    
    if ~tfLegUpdated
      if Ngrps==0
        legend(gdata.axGroupStats,'off');
      else
        legend(hStatLines,STATGRPS(:,1),'Location','SouthEast','interpreter','none');
      end
      tfLegUpdated = true;  
    end
  end
end

guidata(ax,gdata);

function cbkAddCustomGroupViz(ax,grpname,tf)

gdata = guidata(ax);
assert(isvarname(grpname),'Custom group name ''%s'' is not a valid varname.',grpname);
validateattributes(tf,{'logical'},{'vector' 'numel' gdata.Nexp});

% make sure group doesn't already exist
cgInfo = gdata.customGroupVizInfo;
tfEmpty = cellfun(@isempty,cgInfo);
cgNames = cellfun(@(x)x.name,cgInfo(~tfEmpty),'uni',0);
if ismember(grpname,cgNames)
  error('ethogram_plot:groupExists','Custom group ''%s'' already exists.',grpname);
end

% create a new customGroupVizInfo element
newCGVI = struct();
newCGVI.name = grpname;
newCGVI.tf = tf;

% find a place for the new element
if any(tfEmpty)
  idx = find(tfEmpty,1,'first');
else
  idx = numel(tfEmpty)+1;
end
idxMod = mod(idx-1,gdata.CUSTOMGROUPMODN)+1; % idxMod in [1,CUSTOMGROUPMODN]

% draw it
%assert(idx<=gdata.MAXCUSTOMGROUPS,'Cannot support more than %d custom groups.',gdata.MAXCUSTOMGROUPS);
colors = lines(gdata.CUSTOMGROUPMODN);
newCGVI.h = lclDrawCustomGroup(ax,idx,tf,colors(idxMod,:));

% save it
cgInfo{idx,1} = newCGVI;
gdata.customGroupVizInfo = cgInfo;
guidata(ax,gdata);

function cbkRmCustomGroupViz(ax,grpname)
assert(ischar(grpname) && ~isempty(grpname));

gdata = guidata(ax);
cgInfo = gdata.customGroupVizInfo;
Ninfo = numel(cgInfo);
idxRm = zeros(0,1);
for i = 1:Ninfo
  if ~isempty(cgInfo{i}) && strcmp(grpname,cgInfo{i}.name)
    idxRm(end+1,1) = i; %#ok<AGROW>
  end
end

if isempty(idxRm)
  warning('ethogram_plot:nonexistentGroup','Custom group ''%s'' doesn''t exist.',grpname);
elseif isscalar(idxRm)
  delete(cgInfo{idxRm}.h);
  cgInfo{idxRm} = [];
  gdata.customGroupVizInfo = cgInfo;  
else
  assert(false,'Group names should be unique.');
end
guidata(ax,gdata);
  
function h = lclDrawCustomGroup(ax,idx,tf,clr)
% idx: 1-based index for custom group
% tf: indicator vec, one element for each experiment
% clr: 1x3 color

gdata = guidata(ax);
assert(numel(tf)==gdata.Nexp);

%assert(idx<=gdata.MAXCUSTOMGROUPS,'Max number of groups exceeded.');
dx = (gdata.T1-gdata.T0)/250;
idxDiv = floor((idx-1)/gdata.CUSTOMGROUPMODN); % [1,CUSTOMGROUPMODN] -> 0; next CUSTOMGROUPMODN -> 1, etc.
idxMod = mod(idx-1,gdata.CUSTOMGROUPMODN)+1; % idxMod in [1,CUSTOMGROUPMODN]
x = gdata.T0+1+dx*idxDiv;
h = zeros(0,1);
for i = 1:gdata.Nexp
  if tf(i)
    y = i+(idxMod-0.5)/gdata.CUSTOMGROUPMODN;
    h(end+1,1) = plot(ax,x,y,'s','Color',clr,'MarkerFaceColor',clr,'userdata','udCustomGroupMarker'); %#ok<AGROW>
  end
end

function initSelectedExpPatches(ax)
gdata = guidata(ax);
x0 = gdata.T0;
x1 = gdata.T0+10;
xx = [x0 x0 x1 x1];
for i = 1:gdata.Nexp
  yy = [i i+1 i+1 i];
  gdata.hUserSelectionPatches(i) = patch(xx,yy,-2*ones(size(xx)),gdata.EXP_SELECTED_COLOR,...
    'EdgeColor','none','FaceAlpha',gdata.EXP_SELECTED_ALPHA,...
    'HitTest','off','Parent',ax,'Visible','off','userdata','udUserSelectionPatch');
end
guidata(ax,gdata);

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
    pt = get(ax,'currentpoint');
    expnum = floor(pt(1,2));
    gdata = guidata(ax);
    if expnum > gdata.Nexp % Can occur if ctrl-click slightly 'below' axis
      warning('ethogram_plot:noSuchExp','Clicked point beyond range of experiments.');
    else
      gdata = lclToggleExpSelection(gdata,expnum);
    end
    gdata.tfExpSelInProgress = true;
    guidata(ax,gdata);
  otherwise
    patchClickCallback(ax);
end

function cbkFigBDF(ax)
gdata = guidata(ax);
ty = get(ancestor(ax,'Figure'),'SelectionType');
switch ty
  case 'alt' % Ctrl-click
    cbkSetSelectedExps(ax,false(gdata.Nexp,1));
end

function cbkFigWBMF(ax)
gdata = guidata(ax);
if gdata.tfExpSelInProgress
  pt = get(ax,'CurrentPoint');
  expnum = floor(pt(1,2));
  if 0 < expnum && expnum <= gdata.Nexp && ~gdata.tfExpIsUserSelected(expnum)
    set(gdata.hUserSelectionPatches(expnum),'Visible','on');
    gdata.tfExpIsUserSelected(expnum) = true;
  end
end
guidata(ax,gdata);  

function cbkFigWBUF(ax)
gdata = guidata(ax);
if gdata.tfExpSelInProgress
  gdata.tfExpSelInProgress = false;
end
guidata(ax,gdata);

% function cbkFigResized(ax)
% gdata = guidata(ax);
% % delete overlay axis for stats, it doesn't resize in sync with the main
% % axis, leading to incorrect display
% delete(gdata.axGroupStats(ishandle(gdata.axGroupStats)));
% gdata.axGroupStats = [];
% guidata(ax,gdata);

function tf = cbkGetSelectedExps(ax)
gdata = guidata(ax);
tf = gdata.tfExpIsUserSelected;

function cbkSetSelectedExps(ax,tf)
gdata = guidata(ax);
validateattributes(tf,{'logical'},{'vector' 'numel' gdata.Nexp});
gdata.tfExpIsUserSelected = tf;
set(gdata.hUserSelectionPatches(~tf),'Visible','off');
set(gdata.hUserSelectionPatches(tf),'Visible','on');
guidata(ax,gdata);

function gdata = lclToggleExpSelection(gdata,expnum)
assert(ismember(expnum,1:gdata.Nexp));
if gdata.tfExpIsUserSelected(expnum)
  % Was selected, now deselect
  set(gdata.hUserSelectionPatches(expnum),'Visible','off');
  gdata.tfExpIsUserSelected(expnum) = false;
else
  % Was not selected, now select
  set(gdata.hUserSelectionPatches(expnum),'Visible','on');
  gdata.tfExpIsUserSelected(expnum) = true;
end

function patchClickCallback(h)
% Check for double click
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

if gdata.scoresnotjabs
  warning('ethogram_plot_core:nojabs','No jab files specified; clicking on bouts is disabled.');
  return;
end

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
