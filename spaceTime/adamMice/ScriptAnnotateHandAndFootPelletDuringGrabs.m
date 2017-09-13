setuppaths;
dbstop if error;

windowradius = 10;

[annotations,expdirs,annotationi,tmpfile] = AnnotateHandAddFoodPelletDuringGrabs('windowradius',windowradius);

%% break up based on parent directory

parentdirs = cell(1,numel(expdirs));
expnames = cell(1,numel(expdirs));
if ~ispc,
  for i = 1:numel(expdirs),
    [parentdirs{i},expnames{i}] = fileparts(strrep(expdirs{i},'\','/'));
  end
else
  for i = 1:numel(expdirs),
    [parentdirs{i},expnames{i}] = fileparts(expdirs{i});
  end
end

[uniqueparentdirs,~,groupidx] = unique(parentdirs);
groupnames = cell(size(uniqueparentdirs));
for i = 1:numel(uniqueparentdirs),
  [~,groupnames{i}] = fileparts(uniqueparentdirs{i});
end

%% plot

hfig = 1;
figure(hfig);
clf;

colors = lines(numel(groupnames));
h = [];
for groupi = 1:numel(groupnames),
  subplot(1,2,1);
  axis equal;
  hold on;
  pawside = cat(1,annotations(groupidx==groupi).PawSide);
  pelletside = cat(1,annotations(groupidx==groupi).PelletSide);
  d = pawside - pelletside;
  mu = mean(d,1);
  sig = cov(d);
  [a,b,theta] = cov2ell(sig);
  htmp = ellipsedrawpatch(a/2,b/2,mu(1),mu(2),theta);
  set(htmp,'FaceColor',colors(groupi,:),'FaceAlpha',.1,'EdgeColor',colors(groupi,:));
  plot(mu(1),mu(2),'o','Color','k','MarkerFaceColor',colors(groupi,:));
  plot(d(:,1),d(:,2),'.','Color',colors(groupi,:));
  axisalmosttight;
  set(gca,'YDir','reverse');
  title('Side View');

  subplot(1,2,2);
  axis equal;
  hold on;
  pawfront = cat(1,annotations(groupidx==groupi).PawFront);
  pelletfront = cat(1,annotations(groupidx==groupi).PelletFront);
  d = pawfront - pelletfront;
  mu = mean(d,1);
  sig = cov(d);
  [a,b,theta] = cov2ell(sig);
  htmp = ellipsedrawpatch(a/2,b/2,mu(1),mu(2),theta);
  set(htmp,'FaceColor',colors(groupi,:),'FaceAlpha',.1,'EdgeColor',colors(groupi,:));
  plot(mu(1),mu(2),'o','Color','k','MarkerFaceColor',colors(groupi,:));
  h(groupi) = plot(d(:,1),d(:,2),'.','Color',colors(groupi,:));
  axisalmosttight;
  set(gca,'YDir','reverse');
  title('Front View');
end

legend(h,groupnames);

%% plot

hfig = 2;
figure(hfig);
clf;

colors = lines(numel(groupnames));
h = [];
isbadann_pawside = cellfun(@(x) numel(x) ~= 2,{annotations.PawSide});
isbadann_pawfront = cellfun(@(x) numel(x) ~= 2,{annotations.PawFront});

for groupi = 1:numel(groupnames),
  subplot(1,2,1);
  axis equal;
  hold on;
  pawside = cat(1,annotations(groupidx(:)==groupi&~isbadann_pawside(:)).PawSide);
  pelletside = cat(1,annotations(groupidx(:)==groupi&~isbadann_pawside(:)).ip_food);
  d = pawside - pelletside;
  mu = mean(d,1);
  sig = cov(d);
  [a,b,theta] = cov2ell(sig);
  htmp = ellipsedrawpatch(a/2,b/2,mu(1),mu(2),theta);
  set(htmp,'FaceColor',colors(groupi,:),'FaceAlpha',.1,'EdgeColor',colors(groupi,:));
  plot(mu(1),mu(2),'o','Color','k','MarkerFaceColor',colors(groupi,:));
  plot(d(:,1),d(:,2),'.','Color',colors(groupi,:));
  axisalmosttight;
  set(gca,'YDir','reverse');
  title('Side View');

  subplot(1,2,2);
  axis equal;
  hold on;
  pawfront = cat(1,annotations(groupidx(:)==groupi&~isbadann_pawfront(:)).PawFront);
  pelletfront = cat(1,annotations(groupidx(:)==groupi&~isbadann_pawfront(:)).ip_foodfront);
  d = pawfront - pelletfront;
  mu = mean(d,1);
  sig = cov(d);
  [a,b,theta] = cov2ell(sig);
  htmp = ellipsedrawpatch(a/2,b/2,mu(1),mu(2),theta);
  set(htmp,'FaceColor',colors(groupi,:),'FaceAlpha',.1,'EdgeColor',colors(groupi,:));
  plot(mu(1),mu(2),'o','Color','k','MarkerFaceColor',colors(groupi,:));
  h(groupi) = plot(d(:,1),d(:,2),'.','Color',colors(groupi,:));
  axisalmosttight;
  set(gca,'YDir','reverse');
  title('Front View');
end

legend(h,groupnames);

%% for each group, show order of experiments

hfig = 3;

figure(hfig);
clf;

colors = lines(numel(groupnames));
h = [];
isbadann_pawside = cellfun(@(x) numel(x) ~= 2,{annotations.PawSide});
isbadann_pawfront = cellfun(@(x) numel(x) ~= 2,{annotations.PawFront});

for groupi = 1:numel(groupnames),
  

  subplot(1,2,1);
  axis equal;
  hold on;

  groupidxcurr = find(groupidx(:)==groupi&~isbadann_pawside(:));
  pawside = cat(1,annotations(groupidxcurr).PawSide);
  pelletside = cat(1,annotations(groupidxcurr).ip_food);
  
  trialnum = regexp(expnames(groupidxcurr),'v(\d+)$','once','tokens');
  trialnum = str2double([trialnum{:}]);
  [~,ordercurr] = sort(trialnum);
  maxtrial = max(trialnum);
  alpha = (0:maxtrial-1)/(1.1*maxtrial);

  d = pawside - pelletside;  
  mu = mean(d,1);
  sig = cov(d);
  [a,b,theta] = cov2ell(sig);
  htmp = ellipsedrawpatch(a/2,b/2,mu(1),mu(2),theta);
  set(htmp,'FaceColor',colors(groupi,:),'FaceAlpha',.1,'EdgeColor',colors(groupi,:));
  plot(mu(1),mu(2),'o','Color','k','MarkerFaceColor',colors(groupi,:));
  %plot(d(ordercurr,1),d(ordercurr,2),'-','Color',colors(groupi,:)*.2+.8);
  for j = 1:size(d,1),
    plot(d(j,1),d(j,2),'s','MarkerFaceColor',colors(groupi,:)*(1-alpha(j))+alpha(j),'Color',colors(groupi,:),'MarkerSize',15);
    text(d(j,1),d(j,2),num2str(trialnum(j)),'HorizontalAlignment','center');
  end
  %plot(d(:,1),d(:,2),'o','Color',colors(groupi,:),'MarkerSize',16);
  %text(d(:,1),d(:,2),trialnum,'Color',colors(groupi,:),'HorizontalAlignment','center','VerticalAlignment','middle');
  axisalmosttight;
  set(gca,'YDir','reverse');
  title('Side View');

  subplot(1,2,2);
  axis equal;
  hold on;
  
  groupidxcurr = find(groupidx(:)==groupi&~isbadann_pawfront(:));
  pawfront = cat(1,annotations(groupidxcurr).PawFront);
  pelletfront = cat(1,annotations(groupidxcurr).ip_foodfront);
  
  trialnum = regexp(expnames(groupidxcurr),'v(\d+)$','once','tokens');
  trialnum = str2double([trialnum{:}]);
  [~,ordercurr] = sort(trialnum);
  maxtrial = max(trialnum);
  alpha = (0:maxtrial-1)/(1.1*maxtrial);
  
  d = pawfront - pelletfront;
  mu = mean(d,1);
  sig = cov(d);
  [a,b,theta] = cov2ell(sig);
  htmp = ellipsedrawpatch(a/2,b/2,mu(1),mu(2),theta);
  set(htmp,'FaceColor',colors(groupi,:),'FaceAlpha',.1,'EdgeColor',colors(groupi,:));
  h(groupi) = plot(mu(1),mu(2),'o','Color','k','MarkerFaceColor',colors(groupi,:));
  for j = 1:size(d,1),
    plot(d(j,1),d(j,2),'s','MarkerFaceColor',colors(groupi,:)*(1-alpha(j))+alpha(j),'Color',colors(groupi,:),'MarkerSize',14);
    text(d(j,1),d(j,2),num2str(trialnum(j)),'HorizontalAlignment','center');
  end
  axisalmosttight;
  set(gca,'YDir','reverse');
  title('Front View');
end


legend(h,groupnames);
