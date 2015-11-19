%% compare stimuli when a given behavior occurs

%% set up path

if ispc,
  addpath ../../FlyBowlAnalysis;
  rootdatadir = 'C:\Data\JAABA\groundtruth_pBDPGAL4U_data';
else
  addpath /groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis;
  rootdatadir = '/groups/branson/bransonlab/projects/JAABA/data/groundtruth_pBDPGAL4U_data';
end
addpath ../perframe;
addpath ../perframe/compute_perframe_features;
addpath ../misc;
addpath ../filehandling;
outfigdir = '../figures/WildtypeOut';
if ~exist(outfigdir,'dir'),
  mkdir(outfigdir);
end

%% data locations

expnames = dir(fullfile(rootdatadir,'*2012*'));
expnames(~[expnames.isdir]) = [];
expnames = {expnames.name};
expdirs = cellfun(@(x) fullfile(rootdatadir,x),expnames,'UniformOutput',false);
nexps = numel(expdirs);


for i = 1:nexps,
  
  scoresfilestrs1 = dir(fullfile(expdirs{i},'*core*.mat'));
  scoresfilestrs1 = {scoresfilestrs1.name};
  scoresfilestrs1 = scoresfilestrs1(cellfun(@isempty,regexp(scoresfilestrs1,'bak','once')));
  
  if i == 1,
    scoresfilestrs = scoresfilestrs1;
  else
    scoresfilestrs = intersect(scoresfilestrs,scoresfilestrs1);
  end
  
end

scoresfns = cellfun(@(x) x(1:end-4),scoresfilestrs,'UniformOutput',false);
labelfns = regexprep(scoresfns,'scores','labels','preservecase','once');
nbehaviors = numel(labelfns);

tmp = trx.a_mm;
meana = mean([tmp{:}]);
tmp = trx.b_mm;
meanb = mean([tmp{:}]);

%% parameters

labelfn = 'labelsTouch';
framesplot = 'during';
arena_center_mm_x = 0;
arena_center_mm_y = 0; 
arena_radius_mm = 65;

%% load in data

trx = Trx('trxfilestr','registered_trx.mat');
for i = 1:nexps,
  trx.AddExpDir(expdirs{i},'openmovie',false);
end

%% plot within entire arena

hfig = 1;
figure(hfig);
clf;
hold on;

for fly = 1:trx.nflies,
  [starts,ends] = get_interval_ends(trx(fly).(labelfn));
  x = [];
  y = [];
  c = [];
  for i = 1:numel(starts),
    x = [x,trx(fly).x_mm(starts(i):ends(i)),nan];
    y = [y,trx(fly).y_mm(starts(i):ends(i)),nan];
    %c = [c,log(1:ends(i)-starts(i)+1),nan];
    c = [c,linspace(0,1,ends(i)-starts(i)+1),0];
  end
  patch(x,y,c,'EdgeColor','interp','FaceColor','none','linewidth',1);
end

axis equal;

%% plot within one part of the arena

hfig = 2;
figure(hfig);
clf;
hold on;

tmp = linspace(-pi/4,pi/4,100);
plot(arena_radius_mm*cos(tmp),arena_radius_mm*sin(tmp),'k-');

xlim = [arena_radius_mm,-arena_radius_mm];
ylim = [arena_radius_mm,-arena_radius_mm];
for fly = 1:trx.nflies,
  [starts,ends] = get_interval_ends(trx(fly).(labelfn));
  x = [];
  y = [];
  c = [];
  for i = 1:numel(starts),
    xcurr = trx(fly).x_mm(starts(i):ends(i));
    ycurr = trx(fly).y_mm(starts(i):ends(i));
    mux = xcurr(1)-arena_center_mm_x; muy = ycurr(1)-arena_center_mm_y;
    angle = atan2(muy,mux);
    if numel(xcurr) == 1,
      dir = 0;
    else
      angle2 = atan2(ycurr(end)-arena_center_mm_y,xcurr(end)-arena_center_mm_x);
      dangle = modrange(angle2-angle,-pi,pi);
      dir = dangle < 0;
    end
    
    R = [cos(angle),sin(angle);-sin(angle),cos(angle)];
    if dir,
      R = [1,0;0,-1]*R;
    end
    pcurr = R*[xcurr-arena_center_mm_x;ycurr-arena_center_mm_y];
    x = [x,pcurr(1,:),nan];
    y = [y,pcurr(2,:),nan];
    c = [c,linspace(0,1,ends(i)-starts(i)+1),0];
    %c = [c,log(1:ends(i)-starts(i)+1),nan];
    %c = [c,1:ends(i)-starts(i)+1,nan];
    xlim = [min(xlim(1),min(x)),max(xlim(2),max(x))];
    ylim = [min(ylim(1),min(y)),max(ylim(2),max(y))];
  end
  patch(x,y,c,'EdgeColor','interp','FaceColor','none');
  %scatter(x(~isnan(x)),y(~isnan(x)),[],c,'.');
end

xlim(2) = arena_radius_mm;
dx = diff(xlim);
dy = diff(ylim);
xlim = xlim+.05*dx*[-1,1];
ylim = ylim+.05*dy*[-1,1];

axis equal;
switch labelfn,
  case 'labelsBackup',
    axis([47.205168027870400  65.224169839417883  -1.024244770031297   5.700111994984872]);
end

%% locations of other flies

hfig = 3;
figure(hfig);
clf;
hold on;

nbins = 15;
binradius = 6;
binedges = linspace(-binradius,binradius,nbins+1);
bincenters = (binedges(1:end-1)+binedges(2:end))/2;

counts = zeros(nbins,nbins);
meandu = zeros(nbins,nbins);
meandv = zeros(nbins,nbins);

xlim = [arena_radius_mm,-arena_radius_mm];
ylim = [arena_radius_mm,-arena_radius_mm];
for fly = 1:trx.nflies,
  [starts,ends] = get_interval_ends(trx(fly).(labelfn));
  x = [];
  y = [];
  c = [];
  for i = 1:numel(starts),
    fly2 = trx(fly).closestfly_nose2ell_angle_min30to30(starts(i));
    t0 = starts(i)+trx(fly).firstframe-1;
    t1 = ends(i)+1+trx(fly).firstframe-1;
    if t0 < trx(fly2).firstframe || t1 > trx(fly2).endframe,
      continue;
    end
    i0 = t0 + trx(fly2).off;
    i1 = t1 + trx(fly2).off;

    xcurr = trx(fly).x_mm(starts(i):ends(i)+1);
    ycurr = trx(fly).y_mm(starts(i):ends(i)+1);
    thetacurr = trx(fly).theta_mm(starts(i):ends(i)+1);
    xcurr2 = trx(fly2).x_mm(i0:i1) - xcurr;
    ycurr2 = trx(fly2).y_mm(i0:i1) - ycurr;
    ct = cos(thetacurr);
    st = sin(thetacurr);
    u = xcurr2.*ct + ycurr2.*st;
    v = -xcurr2.*st + ycurr2.*ct;
    du = diff(u);
    dv = diff(v);
    [~,idxx] = histc(u(1:end-1),binedges);
    [~,idxy] = histc(v(1:end-1),binedges);
    idxx(idxx>nbins) = nbins;
    idxy(idxy>nbins) = nbins;
    idx = nan(size(idxx));
    tmp = idxx > 0 & idxy > 0;
    if ~any(tmp), continue; end
    idx(tmp) = sub2ind([nbins,nbins],idxy(tmp),idxx(tmp));
    for j = unique(idx(tmp)),
      idxcurr = idx == j;
      meandu(j) = meandu(j) + sum(du(idxcurr));
      meandv(j) = meandv(j) + sum(dv(idxcurr));
      counts(j) = counts(j) + nnz(idxcurr);
    end

%     x = [x,u,nan];
%     y = [y,v,nan];
%     c = [c,linspace(0,1,ends(i)-starts(i)+1),nan];
    %c = [c,log(1:ends(i)-starts(i)+1),nan];
    %c = [c,1:ends(i)-starts(i)+1,nan];
  end
%   xlim = [min(xlim(1),min(x)),max(xlim(2),max(x))];
%   ylim = [min(ylim(1),min(y)),max(ylim(2),max(y))];
%   patch(x,y,c,'EdgeColor','interp','FaceColor','none');
  %scatter(x(~isnan(x)),y(~isnan(x)),[],c,'.');
end

meandu = meandu ./ counts;
meandv = meandv ./ counts;
nodata = counts <= 2;
meanduplot = meandu; meandvplot = meandv;
meanduplot(nodata) = 0;
meandvplot(nodata) = 0;

clf;
[xplot1,yplot1] = meshgrid(bincenters,bincenters);
xplot = [xplot1(:),xplot1(:)+meanduplot(:)*4,nan(numel(xplot1),1)];
yplot = [yplot1(:),yplot1(:)+meandvplot(:)*4,nan(numel(yplot1),1)];
cplot = [counts(:),counts(:),nan(nbins^2,1)];
patch(xplot1(:),yplot1(:),counts(:),'EdgeColor','None','FaceColor','None','marker','.','MarkerEdgeColor','flat','MarkerFaceColor','flat');
patch(xplot',yplot',cplot','EdgeColor','interp','FaceColor','none');

colormap(flipud(gray(256)));

%% 

hfig = 4;
figure(hfig);
clf;
nr = 3;
nc = 4;
hax = createsubplots(nr,nc,.05);


nbins = 15;
binradius = 6;
binedges = linspace(-binradius,binradius,nbins+1);
bincenters = (binedges(1:end-1)+binedges(2:end))/2;
allcounts = cell(1,numel(labelfns));

for behi = 1:numel(labelfns),
  labelfn = labelfns{behi};
  counts = zeros(nbins,nbins);
  allu = [];
  allv = [];
  for fly = 1:trx.nflies,
    [starts,ends] = get_interval_ends(trx(fly).(labelfn));
    for i = 1:numel(starts),
      tmp = trx(fly).closestfly_nose2ell_angle_min30to30;
      if numel(tmp) < ends(i), continue; end
      fly2 = tmp(starts(i));
      t0 = starts(i)+trx(fly).firstframe-1;
      t1 = ends(i)+trx(fly).firstframe-1;
      if t0 < trx(fly2).firstframe || t1 > trx(fly2).endframe,
        continue;
      end
      i0 = t0 + trx(fly2).off;
      i1 = t1 + trx(fly2).off;
      
      xcurr = trx(fly).x_mm(starts(i):ends(i));
      ycurr = trx(fly).y_mm(starts(i):ends(i));
      thetacurr = trx(fly).theta_mm(starts(i):ends(i));
      xcurr2 = trx(fly2).x_mm(i0:i1) - xcurr;
      ycurr2 = trx(fly2).y_mm(i0:i1) - ycurr;
      ct = cos(thetacurr);
      st = sin(thetacurr);
      u = xcurr2.*ct + ycurr2.*st;
      v = -xcurr2.*st + ycurr2.*ct;
      countscurr = hist3([v(:),u(:)],'Edges',{binedges,binedges});
      countscurr = countscurr(1:end-1,1:end-1);
      z = sum(countscurr(:));
      if z == 0,
        continue;
      end
      counts = counts + countscurr/z;
    end    
  end

  allcounts{behi} = counts;
  
%   
%   hfig = 4;
% figure(hfig);
% clf;
% nr = 3;
% nc = 4;
% hax = createsubplots(nr,nc,[.02,.01;.04,.04]);
% 
% 
%   for behi = 1:numel(labelfns),
  
  axes(hax(behi));
  imagesc([bincenters(1),bincenters(end)],[bincenters(1),bincenters(end)],allcounts{behi}/sum(allcounts{behi}(:)));
  axis image;
  hold on
  h(behi) = drawflyo(0,0,0,meana,meanb);
  set(h(behi),'Color',[.99,.99,.99]);
  hti = title(labelfns{behi});
  set(hti,'Interpreter','none');
  colorbar;
  drawnow;
  if mod(behi,nr) > 0,
    set(hax(behi),'XTickLabel',{});
  end
  if behi > nr,
    set(hax(behi),'YTickLabel',{});
  end
  
end