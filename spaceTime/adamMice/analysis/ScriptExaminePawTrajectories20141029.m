%% set up paths

addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/misc;
addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/filehandling;
rootdir = '/tier2/hantman';

% datafiles = {'m119nocnoraw.csv'
%   'm119cnoraw.csv'};
ppfiles = {'data/PostProcessed_M119_20141110.mat'};
%   'data/PostProcessed_M122_20141210.csv'
%   'data/PostProcessed_M130_20141210.csv'};
trxfile = 'data/FixedTrackingResults_M119_20141212.mat';
%   'data/FixedTrackingResults_M122_20141215.mat'
%   'data/FixedTrackingResults_M130_20141214.mat'};

timestamp = now;
savefigdir = sprintf('Trajectories%s',datestr(timestamp,'yyyymmdd'));
if ~exist(savefigdir,'dir'),
  mkdir(savefigdir);
end

%% read in data files

trxinfo = load(trxfile);
expdirs_all = cell(1,numel(trxinfo.moviefiles_all));
for i = 1:numel(trxinfo.moviefiles_all),
  [expdirs_all{i}] = fileparts(trxinfo.moviefiles_all{i});
end

if ~ispc && ~isempty(regexp(expdirs_all{1},'^[A-Z]:','once')),
  expdirs_all = strrep(expdirs_all,'\','/');
  expdirs_all = regexprep(expdirs_all,'^[A-Z]:',rootdir);
end
assert(all(cellfun(@exist,expdirs_all)==7));

trxdata = GetExpTypeFromDir(expdirs_all);

expi = 1;
expdir = trxdata(expi).expdir;
[readframe] = get_readframe_fcn(fullfile(expdir,'movie_comb.avi'));
im = readframe(1);
imsz = size(im);

for i = 1:numel(trxdata),
  trxdata(i).x1 = trxinfo.p_all{i}(:,1);
  trxdata(i).x2 = trxinfo.p_all{i}(:,2);
  trxdata(i).y1 = trxinfo.p_all{i}(:,3);
  trxdata(i).y2 = trxinfo.p_all{i}(:,4);
  td = load(fullfile(trxdata(i).expdir,'trx.mat'));
  fns = fieldnames(td.trx(1).arena);
  for j = 1:numel(fns),
    if ~isempty(regexp(fns{j},'front','once')),
      trxdata(i).(fns{j}) = td.trx(1).arena.(fns{j})+[imsz(2)/2,0];
    else      
      trxdata(i).(fns{j}) = td.trx(1).arena.(fns{j});
    end
  end
  trxdata(i).x1norm = (trxdata(i).x1-trxdata(i).perch(1))/(trxdata(i).food(1)-trxdata(i).perch(1));
  trxdata(i).y1norm = (trxdata(i).y1-trxdata(i).perch(2))/(trxdata(i).food(2)-trxdata(i).perch(2));
  trxdata(i).x2norm = trxdata(i).x2;
  trxdata(i).y2norm = (trxdata(i).y2-trxdata(i).foodfront(2))/(trxdata(i).mouthfront(2)-trxdata(i).foodfront(2));
end

rawdata = [];
for i = 1:numel(ppfiles),
  
  [~,~,ext] = fileparts(ppfiles{i});
  switch ext,
    case '.csv'
      rawdatacurr = ReadRawDataFile(ppfiles{i});
      rawdata = structappend(rawdata,rawdatacurr);
    case '.mat'
      rawdatacurr = load(ppfiles{i});
      rawdata = structappend(rawdata,rawdatacurr.data);
  end
  
end

% line these up
[ism,idx] = ismember({trxdata.id},{rawdata.id}); % AL20141205: warning, changed definition of .id produced by ExpPP due to collisions
if ~all(ism),
  fprintf('The following ids for trajectory info do not have corresponding ids for the behavior data:\n');
  fprintf('%s\n',trxdata(~ism).id);
end
rawdata = rawdata(idx(ism));
trxdata = trxdata(ism);
assert(all(strcmp({trxdata.id},{rawdata.id}))); % AL20141205 etc

%% get behavior names and order

% fnsspecial = {'auto_Grab_success','auto_Grab_successtype'};
% fns = fieldnames(rawdata);
% datafns = setdiff(fns(~cellfun(@isempty,regexp(fns,'^auto'))),fnsspecial);
% nstats = numel(datafns);
% 
% % some behaviors don't have counts
% tmp = regexp(datafns,'^auto_([^_]+)_0$','tokens','once');
% behaviors = unique([tmp{:}]);

behaviors = {'Lift','Handopen','Grab','Sup','Atmouth','Chew'};
behaviorfns = {'auto_GS00_Lift_0','auto_GS00_Handopen_0','auto_GS00_Grab_0','auto_GSSS_Sup_0',...
  'auto_GSSS_Atmouth_0','auto_GSSS_Chew_0'};
nbehaviors = numel(behaviors);

%% make a video showing trajectory

lw = 30;
nframespad = 200;
ms = 30;
fs = 40;
dosave = true;

colors = flipud(jet(128));
colors = colors(round(linspace(1,size(colors,1),nbehaviors+1)),:);

%expi = 3;
expi = 221;
expdir = trxdata(expi).expdir;
[readframe,nframes] = get_readframe_fcn(fullfile(expdir,'movie_comb.avi'));
assert(nframes==numel(trxdata(expi).x1));
ts = nan(1,nbehaviors+1);
ts(1) = 1;
for i = 1:nbehaviors,
  ts(i+1) = rawdata(expi).(behaviorfns{i});
end
assert(~any(isnan(ts)));
assert(numel(unique(ts))==nbehaviors+1);

im = readframe(1);
hfig = 1;
figure(hfig);
set(hfig,'Units','pixels','Position',[10,10,size(im,2)/.96*2,size(im,1)/.96*2],'Color','w');
clf;
hax = axes('Position',[.02,.02,.96,.96]);

him = imagesc(im,[0,255]); axis image;
hold on;
axis off;

fil = normpdf(-5:5,0,1);
fil = fil / sum(fil);
npad = floor(numel(fil)/2);
xsmooth1 = conv(trxdata(expi).x1',fil,'valid');
xsmooth1 = [zeros(1,npad),xsmooth1,zeros(1,npad)];
ysmooth1 = conv(trxdata(expi).y1',fil,'valid');
ysmooth1 = [zeros(1,npad),ysmooth1,zeros(1,npad)];
xsmooth2 = conv(trxdata(expi).x2',fil,'valid');
xsmooth2 = [zeros(1,npad),xsmooth2,zeros(1,npad)];
ysmooth2 = conv(trxdata(expi).y2',fil,'valid');
ysmooth2 = [zeros(1,npad),ysmooth2,zeros(1,npad)];
  
colori = 1;
hprev1 = patch(nan,nan,[0,0,0],'FaceColor','none',...
  'EdgeAlpha','interp','EdgeColor','interp','LineWidth',lw,'FaceVertexCData',[0,0,0],'FaceVertexAlphaData',1);
hprev2 = patch(nan,nan,[0,0,0],'FaceColor','none',...
  'EdgeAlpha','interp','EdgeColor','interp','LineWidth',lw,'FaceVertexCData',[0,0,0],'FaceVertexAlphaData',1);
alphasprev = 1;
colorsprev = colors(colori,:);
tstart = npad+2;


htrx1 = plot(nan,nan,'wo','MarkerSize',ms,'LineWidth',lw);
htrx2 = plot(nan,nan,'wo','MarkerSize',ms,'LineWidth',lw);

htext = text(imsz(2)/2,imsz(1)-20,'','Color','w','HorizontalAlignment','center',...
  'FontWeight','bold','FontSize',fs);

if dosave
  savefile = fullfile(savefigdir,'ExampleTracking.avi');
  if exist(savefile,'file'),
    delete(savefile);
  end
  
  vidobj = VideoWriter(savefile);
  open(vidobj);
end
  
lastt = nan;
for t = tstart:ts(end)+nframespad,
  i = find(ts==t);
  if ~isempty(i),
    lastt = t;
    colori = i;
    set(htext,'String',behaviors{colori-1});
  end
  if t - lastt == 40,
    set(htext,'String','');
  end
  im = readframe(t);
  set(him,'CData',im);
%   for i = 1:numel(hprev1),
%     set(hprev1(i),'EdgeAlpha',get(hprev1(i),'EdgeAlpha')*.975);
%     set(hprev2(i),'EdgeAlpha',get(hprev2(i),'EdgeAlpha')*.975);
%   end
  if t > tstart,
    alphasprev = [alphasprev*.98,1];
    colorsprev = [colorsprev;colors(colori,:)]; %#ok<AGROW>
    UpdateInterpColorLine(hprev1,'x',xsmooth1(tstart:t),...
      'y',ysmooth1(tstart:t),'colors',colorsprev,'alphas',alphasprev);
    UpdateInterpColorLine(hprev2,'x',xsmooth2(tstart:t),...
      'y',ysmooth2(tstart:t),'colors',colorsprev,'alphas',alphasprev);
%     hprev1(end+1) = patch(xsmooth1(t-1:t),ysmooth1(t-1:t),[0,0,0],'FaceColor','none',...
%       'EdgeAlpha',1,'EdgeColor',colors(colori,:),'LineWidth',lw);
%     hprev2(end+1) = patch(xsmooth2(t-1:t),ysmooth2(t-1:t),[0,0,0],'FaceColor','none',...
%       'EdgeAlpha',1,'EdgeColor',colors(colori,:),'LineWidth',lw);
%     hprev1(end+1) = plot(trxdata(expi).x1(t-1:t),trxdata(expi).y1(t-1:t),'-','Color',colors(colori,:),'LineWidth',lw);
%     hprev2(end+1) = plot(trxdata(expi).x2(t-1:t),trxdata(expi).y2(t-1:t),'-','Color',colors(colori,:),'LineWidth',lw);
  end
  set(htrx1,'XData',xsmooth1(t),'YData',ysmooth1(t),'Color',colors(colori,:));
  set(htrx2,'XData',xsmooth2(t),'YData',ysmooth2(t),'Color',colors(colori,:));
  
  drawnow;
  
  if dosave,
    fr = getframe(hax);
    writeVideo(vidobj,fr);
  end
end

if dosave,
  close(vidobj);
end

% 
% for i = nbehaviors:-1:1,
%   plot(trxdata(expi).x1(ts(i):ts(i+1)),trxdata(expi).y1(ts(i):ts(i+1)),'-','Color',colors(i,:),'LineWidth',lw);
%   plot(trxdata(expi).x2(ts(i):ts(i+1)),trxdata(expi).y2(ts(i):ts(i+1)),'-','Color',colors(i,:),'LineWidth',lw);
%   
%   plot(5+[0,10],10*i+[0,0],'-','LineWidth',lw,'Color',colors(i,:));
%   if i < nbehaviors,
%     plot(15+[0,10],10*i+[0,0],'-','LineWidth',lw,'Color',colors(i+1,:));
%   end
%   text(30,10*i+1,behaviors{i},'Color','w','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',14);
% end
% axis off;
% set(hfig,'InvertHardCopy','off');

%% plot one trajectory

colors = prism(nbehaviors);

%expi = 3;
expi = 221;
expdir = trxdata(expi).expdir;
[readframe,nframes] = get_readframe_fcn(fullfile(expdir,'movie_comb.avi'));
assert(nframes==numel(trxdata(expi).x1));
ts = nan(1,nbehaviors+1);
ts(1) = 1;
for i = 1:nbehaviors,
  ts(i+1) = rawdata(expi).(behaviorfns{i});
end

im = readframe(1);
hfig = 1;
figure(hfig);
clf;
set(hfig,'Units','pixels','Position',[10,10,1705,657]);

imagesc(im,[0,255]); axis image;
hold on;
lw = 3;
for i = nbehaviors:-1:1,
  plot(trxdata(expi).x1(ts(i):ts(i+1)),trxdata(expi).y1(ts(i):ts(i+1)),'-','Color',colors(i,:),'LineWidth',lw);
  plot(trxdata(expi).x2(ts(i):ts(i+1)),trxdata(expi).y2(ts(i):ts(i+1)),'-','Color',colors(i,:),'LineWidth',lw);
  
  plot(5+[0,10],10*i+[0,0],'-','LineWidth',lw,'Color',colors(i,:));
  if i < nbehaviors,
    plot(15+[0,10],10*i+[0,0],'-','LineWidth',lw,'Color',colors(i+1,:));
  end
  text(30,10*i+1,behaviors{i},'Color','w','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',14);
end
axis off;
set(hfig,'InvertHardCopy','off');
savefig(fullfile(savefigdir,'SampleTrajectory.png'),hfig,'png');

%% compare cno vs normal


colors = prism(nbehaviors);

[days,~,dayidx] = unique({trxdata.day});
ndays = numel(days);

hfig = 2;
figure(hfig);
clf;
nr = 3;
nc = ceil(ndays/3);
%nr = ceil(sqrt(ndays));
%nc = ceil(ndays/nr);
% nr = ndays;
% nc = 1;
hax = createsubplots(nr,nc,[.01,.01]);
delete(hax(ndays+1:end));
hax = hax(1:ndays);
lw = 1;


dts = nan(0,nbehaviors-1);
for expi = 1:numel(rawdata),  
  if ~rawdata(expi).auto_Grab_success || ...
      ~strcmp(rawdata(expi).auto_Grab_successtype,'successful_1grab'),
    continue;
  end
  
  ts = nan(1,nbehaviors);
  for i = 1:nbehaviors,
    tcurr = rawdata(expi).(behaviorfns{i});
    ts(i) = tcurr;
  end
  dts(end+1,:) = diff(ts);
end
dts(dts<=0) = nan;
meandts = nanmean(dts,1);
stddts = nanstd(dts,1,1);
maxdts = meandts+stddts*5;
ismaxdt = false(1,nbehaviors);
ismaxdt(strcmpi(behaviors,'Atmouth')) = true;

for dayi = 1:ndays,
  
  idxcurr = find(dayidx==dayi);
  hold(hax(dayi),'on');
  iscno = unique([trxdata(idxcurr).iscno]);
  assert(numel(iscno)==1);
  for expi = idxcurr,
    
    if ~rawdata(expi).auto_Grab_success || ~strcmp(rawdata(expi).auto_Grab_successtype,'successful_1grab'),
      continue;
    end
    
    ts = nan(1,nbehaviors+1);
    ts(1) = 1;
    for i = 1:nbehaviors,
      tcurr = rawdata(expi).(behaviorfns{i});
      ts(i+1) = tcurr;
    end
    dtscurr = diff(ts(2:end));
    if any(dtscurr(ismaxdt) > maxdts(ismaxdt)),
      continue;
    end
    
    for i = nbehaviors:-1:1,
      
      plot(hax(dayi),trxdata(expi).x1(ts(i):ts(i+1)),trxdata(expi).y1(ts(i):ts(i+1)),'-','Color',colors(i,:),'LineWidth',lw);
      plot(hax(dayi),trxdata(expi).x2(ts(i):ts(i+1)),trxdata(expi).y2(ts(i):ts(i+1)),'-','Color',colors(i,:),'LineWidth',lw);

    end

  end
  axis(hax(dayi),'image','ij');
  if iscno,
    text(size(im,2)-5,5,sprintf('%s, CNO',days{dayi}),'HorizontalAlignment','right','VerticalAlignment','top','Color','w','Parent',hax(dayi));
  else
    text(size(im,2)-5,5,days{dayi},'HorizontalAlignment','right','VerticalAlignment','top','Color','w','Parent',hax(dayi));
  end
end

set(hfig,'InvertHardCopy','off');

lw = 2;
for i = nbehaviors:-1:1,
  plot(hax(1),5+[0,10],15*i+[0,0],'-','LineWidth',lw,'Color',colors(i,:));
  if i < nbehaviors,
    plot(hax(1),15+[0,10],15*i+[0,0],'-','LineWidth',lw,'Color',colors(i+1,:));
  end
  text(30,15*i+1,behaviors{i},'Color','w','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',12,'Parent',hax(1));
end

set(hax,'XLim',[.5,size(im,2)+.5],'YLim',[.5,size(im,1)+.5],'Color','k','XTick',[],'YTick',[]);
savefig(fullfile(savefigdir,'AllTrajectoriesByDay.png'),hfig,'png');
savefig(fullfile(savefigdir,'AllTrajectoriesByDay.pdf'),hfig,'pdf');

%% plot all trajectories after first lift

hfig = 3;
figure(hfig);
clf;
% nr = ceil(sqrt(ndays));
% nc = ceil(ndays/nr);
nr = 3;
nc = ceil(ndays/nr);
% nr = ndays;
% nc = 1;
hax = createsubplots(nr,nc,.01);
delete(hax(ndays+1:end));
hax = hax(1:ndays);
lw = 1;
colors = jet(100);

for dayi = 1:ndays,
  
  idxcurr = find(dayidx==dayi);
  hold(hax(dayi),'on');
  iscno = unique([trxdata(idxcurr).iscno]);
  assert(numel(iscno)==1);
  
  colori = round(linspace(1,100,numel(idxcurr)));

  h = nan(numel(idxcurr),2);
  for expii = 1:numel(idxcurr),
    expi = idxcurr(expii);
        
    ts = nan(1,2);
    ts(1) = rawdata(expi).auto_GS00_Lift_0;
    ts(2) = rawdata(expi).auto_GS00_Handopen_0;
    if any(isnan(ts)),
      continue;
    end
    
    i = 1;
    h(expii,1) = plot(hax(dayi),trxdata(expi).x1(ts(i):ts(i+1)),trxdata(expi).y1(ts(i):ts(i+1)),'-','Color',colors(colori(expii),:),'LineWidth',lw);
    h(expii,2) = plot(hax(dayi),trxdata(expi).x2(ts(i):ts(i+1)),trxdata(expi).y2(ts(i):ts(i+1)),'-','Color',colors(colori(expii),:),'LineWidth',lw);

  end
  axis(hax(dayi),'image','ij');
  set(hax(dayi),'Children',h(randperm(numel(h))));
  if iscno,
    text(size(im,2)-5,5,sprintf('%s, CNO',days{dayi}),'HorizontalAlignment','right','VerticalAlignment','top','Color','w','Parent',hax(dayi));
  else
    text(size(im,2)-5,5,days{dayi},'HorizontalAlignment','right','VerticalAlignment','top','Color','w','Parent',hax(dayi));
  end
end

set(hfig,'InvertHardCopy','off');
% 
% lw = 2;
% for i = nbehaviors:-1:1,
%   plot(hax(1),5+[0,10],15*i+[0,0],'-','LineWidth',lw,'Color',colors(i,:));
%   if i < nbehaviors,
%     plot(hax(1),15+[0,10],15*i+[0,0],'-','LineWidth',lw,'Color',colors(i+1,:));
%   end
%   text(30,15*i+1,behaviors{i},'Color','w','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',12,'Parent',hax(1));
% end

set(hax,'XLim',[.5,size(im,2)+.5],'YLim',[.5,size(im,1)+.5],'Color','k','XTick',[],'YTick',[]);

savefig(fullfile(savefigdir,'AllReachTrajectoriesByDay.png'),hfig,'png');
savefig(fullfile(savefigdir,'AllReachTrajectoriesByDay.pdf'),hfig,'pdf');

%% plot all trajectories after first lift colored by time

hfig = 4;
figure(hfig);
clf;
nr = 3;
nc = ceil(ndays/nr);
% nr = ndays;
% nc = 1;
hax = createsubplots(nr,nc,.01);
delete(hax(ndays+1:end));
hax = hax(1:ndays);
lw = 1;

dts = ([rawdata.auto_GS00_Handopen_0] - [rawdata.auto_GS00_Lift_0])+1;
maxdt = round(prctile(dts(~isnan(dts)),95));

colors = jet(maxdt);
    


for dayi = 1:ndays,
  
  idxcurr = find(dayidx==dayi);
  hold(hax(dayi),'on');
  iscno = unique([trxdata(idxcurr).iscno]);
  assert(numel(iscno)==1);
  
  h = nan(numel(idxcurr),2);
  for expii = 1:numel(idxcurr),
    expi = idxcurr(expii);
        
    ts = nan(1,2);
    ts(1) = rawdata(expi).auto_GS00_Lift_0;
    ts(2) = rawdata(expi).auto_GS00_Handopen_0;
    if any(isnan(ts)),
      continue;
    end
    ts(2) = min(ts(2),ts(1)+maxdt-1);
        
    i = 1;
    h(expii,1) = PlotInterpColorLine(trxdata(expi).x1(ts(i):ts(i+1)),trxdata(expi).y1(ts(i):ts(i+1)),colors(min(maxdt,1:diff(ts)+1),:),'Parent',hax(dayi));
    h(expii,2) = PlotInterpColorLine(trxdata(expi).x2(ts(i):ts(i+1)),trxdata(expi).y2(ts(i):ts(i+1)),colors(min(maxdt,1:diff(ts)+1),:),'Parent',hax(dayi));

  end
  set(hax(dayi),'Children',h(randperm(numel(h))));

  axis(hax(dayi),'image','ij');
  if iscno,
    text(size(im,2)-5,5,sprintf('%s, CNO',days{dayi}),'HorizontalAlignment','right','VerticalAlignment','top','Color','w','Parent',hax(dayi));
  else
    text(size(im,2)-5,5,days{dayi},'HorizontalAlignment','right','VerticalAlignment','top','Color','w','Parent',hax(dayi));
  end
end

set(hfig,'InvertHardCopy','off');
% 
% lw = 2;
% for i = nbehaviors:-1:1,
%   plot(hax(1),5+[0,10],15*i+[0,0],'-','LineWidth',lw,'Color',colors(i,:));
%   if i < nbehaviors,
%     plot(hax(1),15+[0,10],15*i+[0,0],'-','LineWidth',lw,'Color',colors(i+1,:));
%   end
%   text(30,15*i+1,behaviors{i},'Color','w','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',12,'Parent',hax(1));
% end

set(hax,'XLim',[.5,size(im,2)+.5],'YLim',[.5,size(im,1)+.5],'Color','k','XTick',[],'YTick',[]);
set(hax,'CLim',[0,maxdt-1]);
colorbar('Peer',hax(1),'Location','West','XColor','w','YColor','w');
savefig(fullfile(savefigdir,'AllReachTrajectoriesColorByTime.png'),hfig,'png');

%% plot time after lift vs distance from food

hfig = 5;
figure(hfig);
clf;
% nr = ceil(sqrt(ndays));
% nc = ceil(ndays/nr);
nr = ndays;
nc = 2;
hax = createsubplots(nr,nc,[.05,.025;.025,.01]);
lw = 1;

hax = reshape(hax,[nr,nc]);

mu1 = nan(ndays,maxdt);
sig1 = nan(ndays,maxdt);
mu2 = nan(ndays,maxdt);
sig2 = nan(ndays,maxdt);
maxd = 0;
for dayi = 1:ndays,
  
  idxcurr = find(dayidx==dayi);
  hold(hax(dayi,1),'on');
  hold(hax(dayi,2),'on');
  iscno = unique([trxdata(idxcurr).iscno]);
  assert(numel(iscno)==1);
 
  h = nan(numel(idxcurr),2);
  colors = jet(numel(idxcurr));
  d1 = nan(numel(idxcurr),maxdt);
  d2 = nan(numel(idxcurr),maxdt);
  for expii = 1:numel(idxcurr),
    expi = idxcurr(expii);
        
    ts = nan(1,2);
    ts(1) = rawdata(expi).auto_GS00_Lift_0;
    ts(2) = rawdata(expi).auto_GS00_Handopen_0;
    if any(isnan(ts)),
      continue;
    end
    ts(2) = min(ts(2),ts(1)+maxdt-1);
    
    i = 1;
    d1(expii,1:diff(ts)+1) = sqrt((trxdata(expi).x1(ts(i):ts(i+1))-trxdata(expi).food(1)).^2+...
      (trxdata(expi).y1(ts(i):ts(i+1))-trxdata(expi).food(2)).^2);
    d2(expii,1:diff(ts)+1) = sqrt((trxdata(expi).x2(ts(i):ts(i+1))-trxdata(expi).foodfront(1)).^2+...
      (trxdata(expi).y2(ts(i):ts(i+1))-trxdata(expi).foodfront(2)).^2);
    
    h(expii,1) = plot(hax(dayi,1),1:diff(ts)+1,d1(expii,1:diff(ts)+1),'Color',colors(expii,:));
    h(expii,2) = plot(hax(dayi,2),1:diff(ts)+1,d2(expii,1:diff(ts)+1),'Color',colors(expii,:));

  end
  order = randperm(size(h,1));
  set(hax(dayi,1),'Children',h(order,1));
  set(hax(dayi,2),'Children',h(order,2));
  mu1(dayi,:) = nanmean(d1,1);
  sig1(dayi,:) = nanstd(d1,1,1);
  mu2(dayi,:) = nanmean(d2,1);
  sig2(dayi,:) = nanstd(d2,1,1);
  maxd = max([maxd,max(d1(:)),max(d2(:))]);
    
  if iscno,
    text(maxdt*.99,maxd*.99,sprintf('%s, CNO',days{dayi}),'HorizontalAlignment','right','VerticalAlignment','top','Color','w','Parent',hax(dayi,1));
    text(maxdt*.99,maxd*.99,sprintf('%s, CNO',days{dayi}),'HorizontalAlignment','right','VerticalAlignment','top','Color','w','Parent',hax(dayi,2));
  else
    text(maxdt*.99,maxd*.99,days{dayi},'HorizontalAlignment','right','VerticalAlignment','top','Color','w','Parent',hax(dayi,1));
    text(maxdt*.99,maxd*.99,days{dayi},'HorizontalAlignment','right','VerticalAlignment','top','Color','w','Parent',hax(dayi,2));
  end
  
  if dayi < ndays,
    set(hax(dayi,:),'XTickLabel',{});
  end
  set(hax(dayi,2),'YTickLabel',{});
end

set(hfig,'InvertHardCopy','off');
set(hax,'Color','k','XLim',[0,maxdt],'YLim',[0,maxd]);
linkaxes(hax)
xlabel(hax(end,1),'Time after first lift (frames)');
ylabel(hax(end,1),{'Dist. to food (px)','side view'});
xlabel(hax(end,2),'Time after first lift (frames)');
ylabel(hax(end,2),'front view');
savefig(fullfile(savefigdir,'Dist2FoodAllReaches.png'),hfig,'png');
savefig(fullfile(savefigdir,'Dist2FoodAllReaches.pdf'),hfig,'pdf');

hfig = 6;
figure(hfig);
clf;
hax = createsubplots(1,2,[.05,.1]);
hold(hax(1),'on');
hold(hax(2),'on');
h = nan(1,2);
for dayi = 1:ndays,
  
  idxcurr = find(dayidx==dayi);
  iscno = unique([trxdata(idxcurr).iscno]);
  assert(numel(iscno)==1);
  if iscno,
    color = [1,0,0];
    i = 1;
  else
    color = [1,1,1];
    i = 2;
  end
  patch([1:maxdt,maxdt:-1:1],...
    [mu1(dayi,:)+sig1(dayi,:),fliplr(mu1(dayi,:)-sig1(dayi,:))],color,'LineStyle','none','FaceAlpha',.25,...
    'Parent',hax(1));
  patch([1:maxdt,maxdt:-1:1],...
    [mu2(dayi,:)+sig2(dayi,:),fliplr(mu2(dayi,:)-sig2(dayi,:))],color,'LineStyle','none','FaceAlpha',.25,...
    'Parent',hax(2));
  h(i) = plot(hax(1),1:maxdt,mu1(dayi,:),'-','Color',color,'LineWidth',2);
  plot(hax(2),1:maxdt,mu2(dayi,:),'-','Color',color,'LineWidth',2);
end
set(hax,'Color','k','XLim',[0,maxdt],'YLim',[0,maxd]);
hleg = legend(h,{'CNO','Control'});
set(hleg,'TextColor','w')

set(hfig,'InvertHardCopy','off');
xlabel(hax(1),'Time after first lift (frames)');
xlabel(hax(2),'Time after first lift (frames)');
ylabel(hax(1),'Distance to food, side view (px)');
ylabel(hax(2),'Distance to food, front view (px)');
axisalmosttight([],hax(1));
axisalmosttight([],hax(2));
linkaxes(hax);

savefig(fullfile(savefigdir,'MeanDist2FoodAllReaches.png'),hfig,'png');

%% plot all trajectories after first lift colored by speed

hfig = 7;
figure(hfig);
clf;
%nr = ceil(sqrt(ndays));
nr = 3;
nc = ceil(ndays/nr);
% nr = ndays;
% nc = 1;
hax = createsubplots(nr,nc,.01);
delete(hax(ndays+1:end));
hax = hax(1:ndays);
lw = 1;
colors = jet(256);

gfil = normpdf(-5:5,0,1);
gfil = gfil / sum(gfil);
fil = conv([-1,0,1],gfil);
npad = floor(numel(fil)/2);
for i = 1:numel(trxdata),
  
  dx = conv(trxdata(i).x1',fil,'valid');
  dx = [zeros(1,npad),dx,zeros(1,npad)];
  dy = conv(trxdata(i).y1',fil,'valid');
  dy = [zeros(1,npad),dy,zeros(1,npad)];
  speed = sqrt(dx.^2+dy.^2);
  trxdata(i).speed1 = speed;

  dx = conv(trxdata(i).x2',fil,'valid');
  dx = [zeros(1,npad),dx,zeros(1,npad)];
  dy = conv(trxdata(i).y2',fil,'valid');
  dy = [zeros(1,npad),dy,zeros(1,npad)];
  speed = sqrt(dx.^2+dy.^2);
  trxdata(i).speed2 = speed;
  
end

maxspeed = prctile([trxdata.speed1,trxdata.speed2],99);

set(hfig,'Renderer','zbuffer');
for dayi = 1:ndays,
  
  idxcurr = find(dayidx==dayi);
  hold(hax(dayi),'on');
  iscno = unique([trxdata(idxcurr).iscno]);
  assert(numel(iscno)==1);
  
  h = nan(numel(idxcurr),2);
  for expii = 1:numel(idxcurr),
    expi = idxcurr(expii);
    
    ts = nan(1,2);

    ts(1) = rawdata(expi).auto_GS00_Lift_0;
    ts(2) = rawdata(expi).auto_GS00_Handopen_0;
    if any(isnan(ts)),
      continue;
    end
    ts(2) = min(ts(2),ts(1)+maxdt-1);
    
    colori1 = min(size(colors,1),floor(1+trxdata(expi).speed1(ts(1):ts(2))/maxspeed*(size(colors,1)-1)));
    colori2 = min(size(colors,1),floor(1+trxdata(expi).speed2(ts(1):ts(2))/maxspeed*(size(colors,1)-1)));
    
    i = 1;
    h(expii,1) = PlotInterpColorLine(trxdata(expi).x1(ts(i):ts(i+1)),trxdata(expi).y1(ts(i):ts(i+1)),colors(colori1,:),'Parent',hax(dayi));
    h(expii,2) = PlotInterpColorLine(trxdata(expi).x2(ts(i):ts(i+1)),trxdata(expi).y2(ts(i):ts(i+1)),colors(colori2,:),'Parent',hax(dayi));

  end
  set(hax(dayi),'Children',h(randperm(numel(h))));

  axis(hax(dayi),'image','ij');
  if iscno,
    text(size(im,2)-5,5,sprintf('%s, CNO',days{dayi}),'HorizontalAlignment','right','VerticalAlignment','top','Color','w','Parent',hax(dayi));
  else
    text(size(im,2)-5,5,days{dayi},'HorizontalAlignment','right','VerticalAlignment','top','Color','w','Parent',hax(dayi));
  end
end

set(hfig,'InvertHardCopy','off');
% 
% lw = 2;
% for i = nbehaviors:-1:1,
%   plot(hax(1),5+[0,10],15*i+[0,0],'-','LineWidth',lw,'Color',colors(i,:));
%   if i < nbehaviors,
%     plot(hax(1),15+[0,10],15*i+[0,0],'-','LineWidth',lw,'Color',colors(i+1,:));
%   end
%   text(30,15*i+1,behaviors{i},'Color','w','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',12,'Parent',hax(1));
% end

set(hax,'XLim',[.5,size(im,2)+.5],'YLim',[.5,size(im,1)+.5],'Color','k','XTick',[],'YTick',[]);
set(hax,'CLim',[0,maxspeed]);
colorbar('Peer',hax(1),'Location','East','XColor','w','YColor','w');
savefig(fullfile(savefigdir,'AllReachTrajectoriesColorBySpeed.png'),hfig,'png');

%% plot time after lift vs speed

hfig = 8;
figure(hfig);
clf;
%nr = ceil(sqrt(ndays));
%nc = ceil(ndays/nr);
nr = ndays;
nc = 2;
hax = createsubplots(nr,nc,[.05,.025;.025,.01]);
hax = reshape(hax,[nr,nc]);
lw = 1;

mu = nan(ndays,maxdt);
sig = nan(ndays,maxdt);
maxspeed = 0;
for dayi = 1:ndays,
  
  idxcurr = find(dayidx==dayi);
  hold(hax(dayi,1),'on');
  hold(hax(dayi,2),'on');
  iscno = unique([trxdata(idxcurr).iscno]);
  assert(numel(iscno)==1);
 
  h = nan(numel(idxcurr),2);
  colors = jet(numel(idxcurr));
  speed1 = nan(numel(idxcurr),maxdt);
  speed2 = nan(numel(idxcurr),maxdt);
  for expii = 1:numel(idxcurr),
    expi = idxcurr(expii);
        
    ts = nan(1,2);
    ts(1) = rawdata(expi).auto_GS00_Lift_0;
    ts(2) = rawdata(expi).auto_GS00_Handopen_0;
    if any(isnan(ts)),
      continue;
    end
    ts(2) = min(ts(2),ts(1)+maxdt-1);
        
    i = 1;
    speed1(expii,1:diff(ts)+1) = trxdata(expi).speed1(ts(1):ts(2));
    speed2(expii,1:diff(ts)+1) = trxdata(expi).speed2(ts(1):ts(2));
    maxspeed = max([maxspeed,max(speed1(:)),max(speed2(:))]);
    
    h(expii,1) = plot(hax(dayi,1),1:maxdt,speed1(expii,:),'Color',colors(expii,:));
    h(expii,2) = plot(hax(dayi,2),1:maxdt,speed2(expii,:),'Color',colors(expii,:));

  end

  order = randperm(size(h,1));
  set(hax(dayi,1),'Children',h(order,1));
  set(hax(dayi,2),'Children',h(order,2));
  mu1(dayi,:) = nanmean(speed1,1);
  sig1(dayi,:) = nanstd(speed1,1,1);
  mu2(dayi,:) = nanmean(speed2,1);
  sig2(dayi,:) = nanstd(speed2,1,1);
  
  if iscno,
    text(maxdt*.99,maxspeed*.99,sprintf('%s, CNO',days{dayi}),'HorizontalAlignment','right','VerticalAlignment','top','Color','w','Parent',hax(dayi,1));
    text(maxdt*.99,maxspeed*.99,sprintf('%s, CNO',days{dayi}),'HorizontalAlignment','right','VerticalAlignment','top','Color','w','Parent',hax(dayi,2));
  else
    text(maxdt*.99,maxspeed*.99,days{dayi},'HorizontalAlignment','right','VerticalAlignment','top','Color','w','Parent',hax(dayi,1));
    text(maxdt*.99,maxspeed*.99,days{dayi},'HorizontalAlignment','right','VerticalAlignment','top','Color','w','Parent',hax(dayi,2));
  end

  
  if dayi < ndays,
    set(hax(dayi,:),'XTickLabel',{});
  end
  set(hax(dayi,2),'YTickLabel',{});

  
end

set(hfig,'InvertHardCopy','off');
set(hax,'Color','k','XLim',[0,maxdt],'YLim',[0,maxspeed]);

xlabel(hax(end,1),'Time after first lift (frames)');
ylabel(hax(end,1),{'Speed (px/fr)','side view'});
xlabel(hax(end,2),'Time after first lift (frames)');
ylabel(hax(end,2),'front view');


savefig(fullfile(savefigdir,'SpeedAllReaches.png'),hfig,'png');
savefig(fullfile(savefigdir,'SpeedAllReaches.pdf'),hfig,'pdf');

hfig = 9;
figure(hfig);
clf;
hax = createsubplots(1,2,[.05,.1]);
hold(hax(1),'on');
hold(hax(2),'on');
h = nan(1,2);
for dayi = 1:ndays,
  
  idxcurr = find(dayidx==dayi);
  iscno = unique([trxdata(idxcurr).iscno]);
  assert(numel(iscno)==1);
  if iscno,
    color = [1,0,0];
    i = 1;
  else
    color = [1,1,1];
    i = 2;
  end
  patch([1:maxdt,maxdt:-1:1],...
    [mu1(dayi,:)+sig1(dayi,:),fliplr(mu1(dayi,:)-sig1(dayi,:))],color,'LineStyle','none','FaceAlpha',.25,'Parent',hax(1));
  patch([1:maxdt,maxdt:-1:1],...
    [mu2(dayi,:)+sig2(dayi,:),fliplr(mu2(dayi,:)-sig2(dayi,:))],color,'LineStyle','none','FaceAlpha',.25,'Parent',hax(2));
  h(i) = plot(hax(1),1:maxdt,mu1(dayi,:),'-','Color',color,'LineWidth',2);
  plot(hax(2),1:maxdt,mu2(dayi,:),'-','Color',color,'LineWidth',2);
end
set(hax,'Color','k','XLim',[0,maxdt],'YLim',[0,maxspeed]);
hleg = legend(h,{'CNO','Control'});
set(hleg,'TextColor','w')

set(hfig,'InvertHardCopy','off');
xlabel(hax(1),'Time after first lift (frames)');
xlabel(hax(2),'Time after first lift (frames)');
ylabel(hax(1),'Speed (px/fr), side view');
ylabel(hax(2),'Speed (px/fr), front view');
savefig(fullfile(savefigdir,'MeanSpeedAllReaches.png'),hfig,'png');


%% plot aligned trajectories
% 
% repexpi = nan(1,ndays);
% for dayi = 1:ndays,
% 
%   idxcurr = find(dayidx==dayi);
%   iscno = unique([trxdata(idxcurr).iscno]);
%   assert(numel(iscno)==1);
%  
%   x1 = nan(numel(idxcurr),maxdt);
%   x2 = nan(numel(idxcurr),maxdt);
%   y1 = nan(numel(idxcurr),maxdt);
%   y2 = nan(numel(idxcurr),maxdt);
%   meandt = 0;
%   for expii = 1:numel(idxcurr),
%     expi = idxcurr(expii);
%     
%     ts = nan(1,3);
%     for i = 1:3,
%       tcurr = rawdata(expj).(behaviorfns{i});
%       ts(i) = tcurr;
%     end
%     if any(isnan(ts)),
%       continue;
%     end
%     ts(2) = min(ts(2),ts(1)+maxdt-1);
%         
%     i = 1;
%     x1(expii,1:diff(ts)+1) = trxdata(expi).x1(ts(1):ts(2));
%     x2(expii,1:diff(ts)+1) = trxdata(expi).x2(ts(1):ts(2));
%     y1(expii,1:diff(ts)+1) = trxdata(expi).y1(ts(1):ts(2));
%     y2(expii,1:diff(ts)+1) = trxdata(expi).y2(ts(1):ts(2));
%     meandt = meandt + ts(2)-ts(1)+1;
%   end
%   meandt = meandt/numel(idxcurr);
% 
%   meanx1 = nanmean(x1,1);
%   meanx2 = nanmean(x2,1);
%   meany1 = nanmean(y1,1);
%   meany2 = nanmean(y2,1);
%   D1 = bsxfun(@minus,x1,meanx1).^2+bsxfun(@minus,y1,meany1).^2;
%   meanD1 = nanmean(D1,1);
%   nD1 = bsxfun(@rdivide,D1,meanD1);
%   [mind,expii] = min(nanmean(nD1,2));
%   repexpi(dayi) = idxcurr(expii);
% end
% 

colors = [1,0,1
  0,1,1
  0,1,1
  0,0,1
  1,0,1];

nbehsplot = 2;

lw = 3;
ms = 8;
colormult = .35;

expi = 5;
expj = 41;


idxcurr1 = find(dayidx==1);
idxcurr2 = find(dayidx==2);

for triali = 1:min(numel(idxcurr1),numel(idxcurr2)),

  hfig = triali;

  expi = idxcurr1(triali);
  if rawdata(expi).auto_Grab_success == 0,
    continue;
  end
  expj = idxcurr2(triali);
  if rawdata(expj).auto_Grab_success == 1,
    continue;
  end
  if isnan(rawdata(expj).auto_GS00_Grab_0),
    continue;
  end
  
  figure(hfig);
  clf;
  set(hfig,'Position',[100,100,1600,900]);
  
  hax = createsubplots(1,2,.05);
  
  hold(hax(1),'on');
  hold(hax(2),'on');
  axis(hax,'equal');
  set(hax,'YDir','reverse');
  
  tsi = nan(1,nbehaviors);
  for i = 1:nbehaviors,
    tsi(i) = rawdata(expi).(behaviorfns{i});
  end
  
  tsj = nan(1,nbehaviors);
  for i = 1:nbehaviors,
    tsj(i) = rawdata(expj).(behaviorfns{i});
  end
  
  dt = max(max(tsi(1:nbehsplot+1)-tsi(1)),max(tsj(1:nbehsplot+1)-tsj(1)));
  offy = 50;
  
  h = nan(2,nbehaviors);
  
  for i = 1:nbehsplot,
    t1 = tsi(i);
    t2 = tsi(i+1);
    if t1 - tsi(1) > dt,
      break;
    end
    if isnan(t1),
      break;
    end
    if isnan(t2),
      color = [1,1,1]*(1-colormult);
      t2 = tsi(1) + dt - 1;
    else
      color = colors(i,:)*(1-colormult);
    end
    h(1,i) = plot(hax(1),trxdata(expi).x1(t1:t2),trxdata(expi).y1(t1:t2),'o-','Color',color,'LineWidth',lw,'MarkerSize',ms);
    plot(hax(2),trxdata(expi).x2(t1:t2),trxdata(expi).y2(t1:t2),'o-','Color',color,'LineWidth',lw,'MarkerSize',ms);
  end
  
  for i = 1:nbehsplot,
    t1 = tsj(i);
    t2 = tsj(i+1);
    if t1 - tsj(1) > dt,
      break;
    end
    if isnan(t1),
      break;
    end
    if t2 < t1,
      break;
    end
    
    if isnan(t2),
      color = [1,1,1];
      t2 = tsj(1) + dt - 1;
    else
      color = colors(i,:)*(1-colormult)+colormult;
    end
    h(2,i) = plot(hax(1),trxdata(expj).x1(t1:t2),trxdata(expj).y1(t1:t2)+offy,'o-','Color',color,'LineWidth',lw,'MarkerSize',ms);
    plot(hax(2),trxdata(expj).x2(t1:t2),trxdata(expj).y2(t1:t2)+offy,'o-','Color',color,'LineWidth',lw,'MarkerSize',ms);
  end
  
  for toff = 0:10:dt-1,
    ti = tsi(1)+toff;
    tj = tsj(1)+toff;
    
    i = find(tsi<=ti,1,'last');
    if i > nbehsplot,
      continue;
    end
    if isempty(i),
      colori = [0,0,0];
    else
      colori = colors(i,:)*(1-colormult);
    end
    
    j = find(tsj<=tj,1,'last');
    if j > nbehsplot,
      continue;
    end
    if isempty(j),
      colorj = [0,0,0]+colormult;
    else
      colorj = colors(j,:)*(1-colormult)+colormult;
    end
    
    mux = mean([trxdata(expi).x1(ti),trxdata(expj).x1(tj)]);
    muy = mean([trxdata(expi).y1(ti),trxdata(expj).y1(tj)+offy]);
    
    plot(hax(1),[trxdata(expi).x1(ti),mux],[trxdata(expi).y1(ti),muy],'-','Color',colori);
    plot(hax(1),[trxdata(expj).x1(tj),mux],[trxdata(expj).y1(tj)+offy,muy],'-','Color',colorj);
    
    
    mux = mean([trxdata(expi).x2(ti),trxdata(expj).x2(tj)]);
    muy = mean([trxdata(expi).y2(ti),trxdata(expj).y2(tj)+offy]);
    
    plot(hax(2),[trxdata(expi).x2(ti),mux],[trxdata(expi).y2(ti),muy],'-','Color',colori);
    plot(hax(2),[trxdata(expj).x2(tj),mux],[trxdata(expj).y2(tj)+offy,muy],'-','Color',colorj);
    
    
  end
  
  set(hax,'Color','k');
  axisalmosttight([],hax(1));
  axisalmosttight([],hax(2));
  
  s = cell(2,nbehaviors-1);
  for i = 1:nbehaviors-1,
    s{1,i} = sprintf('%s-%s, control',behaviors{i},behaviors{i+1});
    s{2,i} = sprintf('%s-%s, CNO',behaviors{i},behaviors{i+1});
  end
  hl = legend(h(ishandle(h)),s(ishandle(h)),'Location','NorthWest','TextColor','w');
  
  set(gcf,'InvertHardCopy','off');
  
  title(hax(1),sprintf('Trial %d',triali));
  
  savefig(fullfile(savefigdir,sprintf('AlignedTrajectoriesTrial%d.png',triali)),hfig,'png');
  
end

%%

colors = [1,0,1
  0,1,1
  1,1,0
  0,0,1
  1,0,1];

nbehsplot = 2;

lw = 1;
ms = 8;

behi = 2;

hfig = 1;
figure(hfig);
clf;

hax = createsubplots(ndays,2,.05);
hax = reshape(hax,[ndays,2]);
set(hax,'YDir','reverse','Color','k');

for dayi = 1:ndays,
  idxcurr = find(dayidx==dayi);

  iscno = unique([trxdata(idxcurr).iscno]);
  assert(numel(iscno)==1);
  
  imagesc(im(:,1:imsz(2)/2,:),'Parent',hax(dayi,1));
  imagesc(im(:,imsz(2)/2+1:end,:),'Parent',hax(dayi,2));
  
  hold(hax(dayi,1),'on');
  hold(hax(dayi,2),'on');
  axis(hax(dayi,1),'image','off');
  axis(hax(dayi,2),'image','off');  
  
  for behi = 1:nbehsplot,
  
    for expii = 1:numel(idxcurr),
      expi = idxcurr(expii);

      ts = nan(1,nbehaviors);
      for i = 1:nbehaviors,
        ts(i) = rawdata(expi).(behaviorfns{i});
      end
      if any(isnan(ts(1:nbehsplot+1))),
        continue;
      end
      
%       if ~iscno && (rawdata(expi).auto_Grab_success == 0),
%         continue;
%       end
%       if iscno && (rawdata(expi).auto_Grab_success == 1),
%         continue;
%       end
  
      t1 = ts(behi);
      t2 = ts(behi+1);
      color = colors(behi,:);
      h(behi) = plot(hax(dayi,1),trxdata(expi).x1(t1:t2),trxdata(expi).y1(t1:t2),'-','Color',color,'LineWidth',lw,'MarkerSize',ms);
      plot(hax(dayi,2),trxdata(expi).x2(t1:t2)-imsz(2)/2,trxdata(expi).y2(t1:t2),'-','Color',color,'LineWidth',lw,'MarkerSize',ms);
      
      if behi == nbehsplot,
        plot(hax(dayi,1),trxdata(expi).x1(t1),trxdata(expi).y1(t1),'o','Color',colors(behi+1,:),'LineWidth',3,'MarkerSize',ms);
        plot(hax(dayi,2),trxdata(expi).x2(t1)-imsz(2)/2,trxdata(expi).y2(t1),'o','Color',colors(behi+1,:),'LineWidth',3,'MarkerSize',ms);
      end
    end
    
    if iscno,
      title(hax(dayi,1),sprintf('%s, CNO',days{dayi}));
      title(hax(dayi,2),sprintf('%s, CNO',days{dayi}));
    else
      title(hax(dayi,1),sprintf('%s',days{dayi}));
      title(hax(dayi,2),sprintf('%s',days{dayi}));
    end
    
  end
  
  
end

savefig(fullfile(savefigdir,'HandOpenTrajectories.png'),gcf,'png');

%% 

fps = 500;
toff = round(50/1000*fps);

imshandopen = cell(2,ndays);
for dayi = 1:ndays,
  idxcurr = find(dayidx==dayi);

  for expii = 1:numel(idxcurr),
    expi = idxcurr(expii);

    ts = nan(1,nbehaviors);
    for i = 1:nbehaviors,
      ts(i) = rawdata(expi).(behaviorfns{i});
    end
    if any(isnan(ts(1:nbehsplot+1))),
      continue;
    end
      
%       if ~iscno && (rawdata(expi).auto_Grab_success == 0),
%         continue;
%       end
%       if iscno && (rawdata(expi).auto_Grab_success == 1),
%         continue;
%       end

    behi = nbehsplot;
    t1 = ts(behi);
    
    readframe = get_readframe_fcn(fullfile(trxdata(expi).expdir,'movie_comb.avi'));
    imopen = readframe(t1);
    impre = readframe(t1-toff);
    if isempty(imshandopen{1,dayi}),
      imshandopen{1,dayi} = impre(:,:,1);
      imshandopen{2,dayi} = imopen(:,:,1);
    else
      imshandopen{1,dayi}(:,:,end+1) = impre(:,:,1);
      imshandopen{2,dayi}(:,:,end+1) = imopen(:,:,1);
    end
  end
  
end

for dayi = 1:ndays,
  figure(dayi);
  clf;
  
  nax = size(imshandopen{1,dayi},3);
  nc = ceil(sqrt(nax));
  nr = ceil(nax/nc);
  hax = createsubplots(nr,nc,.01);
  delete(hax(nax+1:end));
  for i = 1:nax,
    imagesc(imshandopen{1,dayi}(:,1:imsz(2)/2,i),'Parent',hax(i),[0,255]);
    axis(hax(i),'image','off');
  end
  
  colormap gray;
  savefig(fullfile(savefigdir,sprintf('HandOpenSensoryInfoDay%s.png',days{dayi})),dayi,'png');
  
end


%% plot aligned trajectories


colors = [1,0,1
  0,1,1
  1,1,0
  0,0,1
  0,1,0];

nbehsplot = 2;

lw = 5;
ms = 20;
colormult = .35;

expi = 1;
expj = 50;


idxcurr1 = find(dayidx==1|dayidx==3);
idxcurr2 = find(dayidx==2);

for ii = 1:numel(idxcurr1),
  expi = idxcurr1(ii);
  if rawdata(expi).auto_Grab_success == 0,
    continue;
  end
  
  for jj = 1:numel(idxcurr2),

    expj = idxcurr2(jj);
    if rawdata(expj).auto_Grab_success == 1,
      continue;
    end
    if isnan(rawdata(expj).auto_GS00_Grab_0),
      continue;
    end

    
    tsi = nan(1,nbehaviors);
    for i = 1:nbehaviors,
      tsi(i) = rawdata(expi).(behaviorfns{i});
    end
    
    tsj = nan(1,nbehaviors);
    for i = 1:nbehaviors,
      tsj(i) = rawdata(expj).(behaviorfns{i});
    end
    
%     if tsj(nbehsplot)-tsj(1)-(tsi(nbehsplot)-tsi(1)) < 40/1000*fps || ...
%         tsj(nbehsplot)-tsj(1)-(tsi(nbehsplot)-tsi(1)) > 200/1000*fps,
%       continue;
%     end
    
    hfig = 1;
    figure(hfig);
    clf;
    
    set(hfig,'Units','pixels','Position',[100,100,1600,900],'Color','k','InvertHardCopy','off');
    
    hax = createsubplots(2,1,.05);
    
    for i = 1:2,
      hold(hax(i),'on');
    end
    
    dt = max(max(tsi(1:nbehsplot+1)-tsi(1)),max(tsj(1:nbehsplot+1)-tsj(1)));
    offy = 50;
    
    h = nan(2,nbehaviors);
    
    for i = 1:nbehaviors-1,
      t1 = tsi(i);
      t2 = tsi(i+1);
      if t1 - tsi(1) > dt,
        break;
      end
      if isnan(t1),
        break;
      end
      if t1 > t2,
        continue;
      end
      if isnan(t2),
        color = [1,1,1];
        t2 = tsi(1) + dt - 1;
      else
        color = colors(i,:);
      end
      h(1,i) = plot(hax(1),(t1-tsi(1):t2-tsi(1))/fps*1000,trxdata(expi).x1(t1:t2),'-','Color',color,'LineWidth',lw,'MarkerSize',ms);
    end
    
    minx = min(min(trxdata(expi).x1(tsi(1):tsi(1)+dt)),min(trxdata(expj).x1(tsj(1):tsj(1)+dt)));
    maxx = max(max(trxdata(expi).x1(tsi(1):tsi(1)+dt)),max(trxdata(expj).x1(tsj(1):tsj(1)+dt)));
    
    for i = 1:nbehaviors-1,
      t1 = tsj(i);
      t2 = tsj(i+1);
      
      if t1 > t2,
        continue;
      end
      
      if t1 - tsj(1) >= dt,
        break;
      end
      if isnan(t1),
        break;
      end
      if t2 < t1,
        break;
      end
      
      if isnan(t2),
        color = [1,1,1];
        t2 = tsj(1) + dt - 1;
      else
        color = colors(i,:);
      end
      h(2,i) = plot(hax(2),(t1-tsj(1):t2-tsj(1))/fps*1000,trxdata(expj).x1(t1:t2),'-','Color',color,'LineWidth',lw,'MarkerSize',ms);
    end
    
    dx = maxx-minx;
    ylim = [minx,maxx]+.05*dx*[-1,1];
    set(hax,'Color','k','XLim',[0,dt/fps*1000],'YLim',ylim);
    
    tgrab = (tsi(nbehsplot)-tsi(1))/fps*1000;
    
    plot(hax(1),[tgrab,tgrab],ylim,'w--','LineWidth',lw);
    plot(hax(2),[tgrab,tgrab],ylim,'w--','LineWidth',lw);
    
    plot(hax(1),tgrab,trxdata(expi).x1(tsi(nbehsplot)),'wo','MarkerFaceColor',colors(nbehsplot,:),'LineWidth',lw,'MarkerSize',ms);
    
    tpre = tgrab - 50;
    %plot(hax(1),[tpre,tpre],ylim,'w--');
    plot(hax(1),tpre,trxdata(expi).x1(tsi(nbehsplot)-round(50*fps/1000)),'wo','LineWidth',lw,'MarkerSize',ms);
    
    tgrab = (tsj(nbehsplot)-tsj(1))/fps*1000;
    %plot(hax(2),[tgrab,tgrab],ylim,'w-');
    plot(hax(2),tgrab,trxdata(expj).x1(tsj(nbehsplot)),'wo','MarkerFaceColor',colors(nbehsplot,:),'LineWidth',lw,'MarkerSize',ms);
    tpre = tgrab - 50;
    plot(hax(2),tpre,trxdata(expj).x1(tsj(nbehsplot)-round(50*fps/1000)),'wo','LineWidth',lw,'MarkerSize',ms);
    xlabel(hax(2),'Time since lift (ms)');
    ylabel(hax(2),'Paw x-coordinate (px)');
    set(hax,'XColor','w','YColor','w');
    
    drawnow

    savefig(fullfile(savefigdir,sprintf('TwoExampleTrajectoriesAligned_%02d_%02d.png',expi,expj)),hfig,'png');
    
  end
  
end
    
    
%% 

ntrials = numel(rawdata);

fnsspecial = {'auto_Grab_success','auto_Grab_successtype'};
fns = fieldnames(rawdata);
datafns = setdiff(fns(~cellfun(@isempty,regexp(fns,'^auto'))),fnsspecial);
nstats = numel(datafns);

% some behaviors don't have counts
tmp = regexp(datafns,'^auto_([^_]+)_0$','tokens','once');
behaviors = unique([tmp{:}]);
nbehaviors = numel(behaviors);
for i = 1:nbehaviors,
  j = 0;
  fnprev = '';
  fncount = sprintf('auto_%s_num',behaviors{i});
  if ismember(fncount,datafns),
    continue;
  end
  count = zeros(ntrials,1);
  while true,
    fn = sprintf('auto_%s_%d',behaviors{i},j);
    if ~ismember(fn,datafns),
      break;
    end
    tmp = ~isnan([rawdata.(fn)]);
    if j > 0,
      tmp = tmp & [rawdata.(fn)]~=[rawdata.(fnprev)];
    end
    count(tmp) = count(tmp) + 1;
    j = j+1;
    fnprev = fn;
  end
  for j = 1:ntrials,
    rawdata(j).(fncount) = count(j); %#ok<SAGROW>
  end
end

datafns = setdiff(fns(~cellfun(@isempty,regexp(fns,'^auto'))),fnsspecial);
nstats = numel(datafns);

% set to nan when del = 0
for i = find(~cellfun(@isempty,regexp(datafns,'^auto_del','once'))),
  fn = datafns{i};
  for j = 1:ntrials,
    if rawdata(j).(fn) == 0,
      rawdata(j).(fn) = nan;
    end
  end
end

X = nan(ntrials,nstats);
for i = 1:nstats,
  X(:,i) = [rawdata.(datafns{i})];
end

% find redundant stats
d = inf(nstats,nstats);
n = nan(nstats,1);
idxreplace = zeros(1,nstats);
for i = 1:nstats,
  idx = ~isnan(X(:,i));
  n(i) = nnz(idx);
  for j = 1:nstats,
    if i == j,
      continue;
    end
    d(i,j) = nnz(X(idx,j) ~= X(idx,i)) / nnz(idx);
    if idxreplace(j) > 0,
      continue;
    end
    if d(i,j) == 0,
      idxreplace(i) = j;
      idxreplace(idxreplace==i) = j;
      fprintf('Can compute %d: %s from %d: %s\n',i,datafns{i},idxreplace(i),datafns{idxreplace(i)});
    end
  end
end

% is there any information in the nan pattern?
idx = ~cellfun(@isempty,regexp(datafns,'num$','once'));
countfns = datafns(idx);

Z = nan(ntrials,numel(countfns)+2);
for i = 1:numel(countfns),
  Z(:,i) = [rawdata.(countfns{i})];
end
Z(:,end-1) = double([rawdata.success]);
[successtypes,~,Z(:,end)] = unique({rawdata.successtype});

mindiff = inf(1,nstats);
idxnanpattern = zeros(6,nstats);
for i = find(idxreplace>0),
  isdata = ~isnan(X(:,i));
  
  for j1 = 1:size(Z,2),
    
    for j2 = 1:size(Z,2),
      
      for k1 = min(Z(:,j1)):max(Z(:,j1)),

        if j1 == size(Z,2),
          isdata1 = Z(:,j1) == k1;
        else
          isdata1 = (Z(:,j1) >= k1);
        end
        
        for k2 = min(Z(:,j2)):max(Z(:,j2)),

          if j2 == size(Z,2),
            isdata2 = Z(:,j2) == k2;
          else
            isdata2 = (Z(:,j2) >= k2);
          end

          for l1 = 1:2,
            
            if l1 == 1,
              isdata1x = isdata1;
            else
              isdata1x = ~isdata1;
            end
            for l2 = 1:2,
              if l2 == 1,
                isdata2x = isdata2;
              else
                isdata2x = ~isdata2;
              end

              isdata3 = isdata1x & isdata2x;
              
              nmismatch = nnz(isdata ~= isdata3);
              if nmismatch < mindiff(i),
                idxnanpattern(:,i) = [j1,j2,k1,k2,l1,l2];
                mindiff(i) = nmismatch;
              end
            end
          end
        end
      end
    end
  end
end

fprintf('Max info in nan patterns = %d / %d\n',max(mindiff(idxreplace>0)),ntrials);

datafns = datafns(idxreplace==0);
nstats = numel(datafns);

X = nan(ntrials,nstats);
for i = 1:nstats,
  X(:,i) = [rawdata.(datafns{i})];
end

%% compare features and various statistics

% success or not
y = [rawdata.success];

hfig = CompareStatsAcrossGroups(X,y,datafns,{'Failure','Success'},'title','Success vs failure, all trials');

savefig(sprintf('SortedStatisticsPredictingSuccess_AllTrials%s.pdf',datestr(timestamp,'yyyymmdd')),hfig,'pdf');
savefig(sprintf('SortedStatisticsPredictingSuccess_AllTrials%s.png',datestr(timestamp,'yyyymmdd')),hfig,'png');

% cno vs not cno
y = [rawdata.iscno];
hfig = CompareStatsAcrossGroups(X,y,datafns,{'No CNO','CNO'},'title','No CNO vs CNO, all trials');

savefig(sprintf('SortedStatisticsPredictingCNO_AllTrials%s.pdf',datestr(timestamp,'yyyymmdd')),hfig,'pdf');
savefig(sprintf('SortedStatisticsPredictingCNO_AllTrials%s.png',datestr(timestamp,'yyyymmdd')),hfig,'png');

% success CNO vs success not CNO
idx = [rawdata.success] == 1;
y = [rawdata.iscno];
hfig = CompareStatsAcrossGroups(X(idx,:),y(idx),datafns,{'No CNO','CNO'},'title','No CNO vs CNO, success trials');

savefig(sprintf('SortedStatisticsPredictingCNO_SuccessTrials%s.pdf',datestr(timestamp,'yyyymmdd')),hfig,'pdf');
savefig(sprintf('SortedStatisticsPredictingCNO_SuccessTrials%s.png',datestr(timestamp,'yyyymmdd')),hfig,'png');

% failure CNO vs failure not CNO
idx = [rawdata.success] == 0;
y = [rawdata.iscno];
hfig = CompareStatsAcrossGroups(X(idx,:),y(idx),datafns,{'No CNO','CNO'},'title','No CNO vs CNO, failure trials');

savefig(sprintf('SortedStatisticsPredictingCNO_FailureTrials%s.pdf',datestr(timestamp,'yyyymmdd')),hfig,'pdf');
savefig(sprintf('SortedStatisticsPredictingCNO_FailureTrials%s.png',datestr(timestamp,'yyyymmdd')),hfig,'png');

% correlate trial to statistics
X1 = [X,double([rawdata.success])'];
datafns1 = [datafns,{'success'}];
y = [rawdata.trial];

hfig = MeasureStatDependencies(X1,y,datafns1,'Effect of trial number, all trials');

savefig(sprintf('SortedStatisticsVsTrialNumber_AllTrials%s.pdf',datestr(timestamp,'yyyymmdd')),hfig,'pdf');
savefig(sprintf('SortedStatisticsVsTrialNumber_AllTrials%s.png',datestr(timestamp,'yyyymmdd')),hfig,'png');

idx = [rawdata.success]==1;
hfig = MeasureStatDependencies(X(idx,:),y(idx),datafns,'Effect of trial number, success trials');

savefig(sprintf('SortedStatisticsVsTrialNumber_SuccessTrials%s.pdf',datestr(timestamp,'yyyymmdd')),hfig,'pdf');
savefig(sprintf('SortedStatisticsVsTrialNumber_SuccessTrials%s.png',datestr(timestamp,'yyyymmdd')),hfig,'png');

idx = [rawdata.success]==0;
hfig = MeasureStatDependencies(X(idx,:),y(idx),datafns,'Effect of trial number, failure trials');

savefig(sprintf('SortedStatisticsVsTrialNumber_FailureTrials%s.pdf',datestr(timestamp,'yyyymmdd')),hfig,'pdf');
savefig(sprintf('SortedStatisticsVsTrialNumber_FailureTrials%s.png',datestr(timestamp,'yyyymmdd')),hfig,'png');

idx = [rawdata.iscno]==1;
hfig = MeasureStatDependencies(X1(idx,:),y(idx),datafns1,'Effect of trial number, CNO trials');

savefig(sprintf('SortedStatisticsVsTrialNumber_CNOTrials%s.pdf',datestr(timestamp,'yyyymmdd')),hfig,'pdf');
savefig(sprintf('SortedStatisticsVsTrialNumber_CNOTrials%s.png',datestr(timestamp,'yyyymmdd')),hfig,'png');

idx = [rawdata.iscno]==0;
hfig = MeasureStatDependencies(X1(idx,:),y(idx),datafns1,'Effect of trial number, no CNO trials');

savefig(sprintf('SortedStatisticsVsTrialNumber_noCNOTrials%s.pdf',datestr(timestamp,'yyyymmdd')),hfig,'pdf');
savefig(sprintf('SortedStatisticsVsTrialNumber_noCNOTrials%s.png',datestr(timestamp,'yyyymmdd')),hfig,'png');

%% how much does a single date contribute?

% correlate trial to statistics
X1 = [X,double([rawdata.success])'];
datafns1 = [datafns,{'success'}];
y = [rawdata.trial];

[uniquedates,~,dateidx] = unique({rawdata.date});

idx = [rawdata.iscno]==1;
[tmp,~,groupidx] = unique(dateidx(idx));
groupnames = uniquedates(tmp);

hfig = MeasureStatDependenciesPerGroup(X1(idx,:),y(idx),groupidx,datafns1,groupnames,'Effect of trial number per date, CNO trials');

savefig(sprintf('SortedStatisticsVsTrialNumber_PerGroup_CNOTrials%s.pdf',datestr(timestamp,'yyyymmdd')),hfig,'pdf');
savefig(sprintf('SortedStatisticsVsTrialNumber_PerGroup_CNOTrials%s.png',datestr(timestamp,'yyyymmdd')),hfig,'png');


idx = [rawdata.iscno]==0;
[tmp,~,groupidx] = unique(dateidx(idx));
groupnames = uniquedates(tmp);

hfig = MeasureStatDependenciesPerGroup(X1(idx,:),y(idx),groupidx,datafns1,groupnames,'Effect of trial number per date, no CNO trials');

savefig(sprintf('SortedStatisticsVsTrialNumber_PerGroup_noCNOTrials%s.pdf',datestr(timestamp,'yyyymmdd')),hfig,'pdf');
savefig(sprintf('SortedStatisticsVsTrialNumber_PerGroup_noCNOTrials%s.png',datestr(timestamp,'yyyymmdd')),hfig,'png');

%% make a table showing all statistics, all groups

ys = nan(0,numel(rawdata));
idxs = false(0,numel(rawdata));
ynames = {};
idxnames = {};

ys(end+1,:) = [rawdata.success];
ynames{end+1} = 'Success';
ys(end+1,:) = [rawdata.iscno];
ynames{end+1} = 'CNO';

idxs(end+1,:) = true(1,numel(rawdata));
idxnames{end+1} = 'All trials';
idxs(end+1,:) = [rawdata.success];
idxnames{end+1} = 'Success trials';
idxs(end+1,:) = ~[rawdata.success];
idxnames{end+1} = 'Failure trials';
idxs(end+1,:) = [rawdata.iscno];
idxnames{end+1} = 'CNO trials';
idxs(end+1,:) = ~[rawdata.iscno];
idxnames{end+1} = 'No CNO trials';

pairsignore = {'Success','Success trials'
  'Success','Failure trials'
  'CNO','CNO trials'
  'CNO','No CNO trials'};

hfig = CompareStatsAcrossMultipleGroups(X,ys,ynames,idxs,idxnames,datafns,pairsignore);
