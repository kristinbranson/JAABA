%% compare behavior detections for different wild-type strains

%% set up path

if ispc,
  addpath ../../FlyBowlAnalysis;
else
  addpath /groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis;
end
addpath ../perframe;
addpath ../perframe/compute_perframe_features;
addpath ../misc;
addpath ../filehandling;
outfigdir = '../figures/WildtypeOut';
if ~exist(outfigdir,'dir'),
  mkdir(outfigdir);
end

[~,computername] = system('hostname');
computername = strtrim(computername);

switch computername,
  
  case 'kabram-ws.janelia.priv',
    rootdatadir= '/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/experiments/wildtype/results';
  
  otherwise,

    rootdatadir = '../experiments/wildtype/results';
  
end


if ispc,
  addpath C:\Code\FlyBowlAnalysis;
else
  addpath /groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis;
end

%% data locations

expnames = dir(fullfile(rootdatadir,'*201206*'));
expnames(~[expnames.isdir]) = [];
expnames = {expnames.name};
expdirs = cellfun(@(x) fullfile(rootdatadir,x),expnames,'UniformOutput',false);
nexps = numel(expdirs);

scoresfilestrs = dir(fullfile(expdirs{1},'*cores*.mat'));
scoresfilestrs = {scoresfilestrs.name};
scoresfns = cellfun(@(x) x(1:end-4),scoresfilestrs,'UniformOutput',false);
labelfns = regexprep(scoresfns,'scores','labels','preservecase','once');
nbehaviors = numel(labelfns);


%% parameters

linenamecodes = {
  'EXT_CSMH'             'CantonS MH'
  'EXT_CantonS_1220002'  'CantonS UH'
  'EXT_DL'               'Dickinson'
  'FCF_cantons_1500002'  'CantonS GR'
  'UAH_R-1220001'        'OregonR'
  'UAH_R-1220003'        'Berlin'
  };

behaviorcodes = {
  'labels_Walk'          'Walk'
  'labels_Stops'         'Stop'
  'labels_Jump'          'Jump'
  'labels_Righting'      'Righting'
  'labelsWingGrooming'   'Wing grooming'
  'labelsCrabwalk'       'Crabwalk'
  'labelsBackup'         'Backup'
  'labels_pivot_tail'    'Tail pivot turn'
  'labelsTouch'          'Touch'
  'labels_Chasev7'       'Chase'
  };


print_linenames_order = {'Berlin','CantonS MH','CantonS UH','CantonS GR','Dickinson','OregonR'};

%% compute statistics of behavior detections

meanfractime = nan(nexps,nbehaviors);
meanfractime_male = nan(nexps,nbehaviors);
meanfractime_female = nan(nexps,nbehaviors);

for expi = 1:nexps,
  
  % current experiment
  try
    trx = Trx('trxfilestr','registered_trx.mat','moviefilestr','movie.ufmf','perframedir','perframe');
    trx.AddExpDir(expdirs{expi},'openmovie',false);
  catch ME
    warning(getReport(ME));
    continue;
  end
  
  % number of frames analyzed
  nframes_total_perfly = trx.nframes;
  nframes_total = sum(nframes_total_perfly);
  nframes_male_perfly = cellfun(@(x) nnz(strcmp(x,'M')),trx.sex);
  nframes_female_perfly = cellfun(@(x) nnz(strcmp(x,'F')),trx.sex);
  nframes_male = sum(nframes_male_perfly);
  nframes_female = sum(nframes_female_perfly);

  % in per-fly std computation, ignore trajectories that are short
  ignoreidx = nframes_total_perfly <= 500;

  fprintf('%s:\n',expnames{expi});

  
  for behi = 1:nbehaviors,

    % mean, both sexes
    nframes_pos = nnz(cell2mat(trx.(labelfns{behi})));
    nframes_pos_perfly = cellfun(@nnz,trx.(labelfns{behi}));
    fractime_perfly = nframes_pos_perfly ./ nframes_total;
    meanfractime(expi,behi) = nframes_pos / nframes_total;

    % mean, per sex
    nframes_male_pos_perfly = nan(1,trx.nflies);
    nframes_female_pos_perfly = nan(1,trx.nflies);
    for fly = 1:trx.nflies,
      idxmale = strcmp(trx(fly).sex,'M');
      idxfemale = strcmp(trx(fly).sex,'F');
      nframes_male_pos_perfly(fly) = nnz(trx(fly).(labelfns{behi})(idxmale));
      nframes_female_pos_perfly(fly) = nnz(trx(fly).(labelfns{behi})(idxfemale));
    end
    nframes_male_pos = sum(nframes_male_pos_perfly);
    nframes_female_pos = sum(nframes_female_pos_perfly);
    meanfractime_male(expi,behi) = nframes_male_pos / nframes_male;
    meanfractime_female(expi,behi) = nframes_female_pos / nframes_female;

    % std TODO
    
    fprintf('%s: b = %f, m = %f, f = %f\n',labelfns{behi},...
      meanfractime(expi,behi),...
      meanfractime_male(expi,behi),meanfractime_female(expi,behi));
    
  end
end

%% z-score the data so that we can plot on the same axes

zmu = nanmean(meanfractime,1);
zsigma = nanstd(meanfractime,1,1);
zmeanfractime = bsxfun(@rdivide,bsxfun(@minus,meanfractime,zmu),zsigma);
zmu_male = nanmean(meanfractime_male,1);
zsigma_male = nanstd(meanfractime_male,1,1);
zmeanfractime_male = bsxfun(@rdivide,bsxfun(@minus,meanfractime_male,zmu_male),zsigma_male);
zmu_female = nanmean(meanfractime_female,1);
zsigma_female = nanstd(meanfractime_female,1,1);
zmeanfractime_female = bsxfun(@rdivide,bsxfun(@minus,meanfractime_female,zmu_female),zsigma_female);

%% compute per-strain data

for i = 1:nexps,
  metadatacurr = parseExpDir(expdirs{i});
  if i == 1,
    metadata = repmat(metadatacurr,[1,nexps]);
  else
    metadata(i) = metadatacurr;
  end
end

line_names = {metadata.line};
[unique_linenames,~,lineidx] = unique(line_names);

nlines = numel(unique_linenames);

meanfractime_line = nan(nlines,nbehaviors);
meanfractime_male_line = nan(nlines,nbehaviors);
meanfractime_female_line = nan(nlines,nbehaviors);
stdfractime_line = nan(nlines,nbehaviors);
stdfractime_male_line = nan(nlines,nbehaviors);
stdfractime_female_line = nan(nlines,nbehaviors);

for linei = 1:nlines,
  meanfractime_line(linei,:) = nanmean(meanfractime(lineidx==linei,:),1);
  stdfractime_line(linei,:) = nanstd(meanfractime(lineidx==linei,:),1,1);
  meanfractime_male_line(linei,:) = nanmean(meanfractime_male(lineidx==linei,:),1);
  stdfractime_male_line(linei,:) = nanstd(meanfractime_male(lineidx==linei,:),1,1);
  meanfractime_female_line(linei,:) = nanmean(meanfractime_female(lineidx==linei,:),1);
  stdfractime_female_line(linei,:) = nanstd(meanfractime_female(lineidx==linei,:),1,1);
end

zmeanfractime_line = bsxfun(@rdivide,bsxfun(@minus,meanfractime_line,zmu),zsigma);
zmeanfractime_male_line = bsxfun(@rdivide,bsxfun(@minus,meanfractime_male_line,zmu_male),zsigma_male);
zmeanfractime_female_line = bsxfun(@rdivide,bsxfun(@minus,meanfractime_female_line,zmu_female),zsigma_female);
zstdfractime_line = bsxfun(@rdivide,stdfractime_line,zsigma);
zstdfractime_male_line = bsxfun(@rdivide,stdfractime_male_line,zsigma_male);
zstdfractime_female_line = bsxfun(@rdivide,stdfractime_female_line,zsigma_female);

%% plot

hfig = 1;
figure(hfig);
clf;
set(hfig,'Position',[80 94 1716 854]);
colors = jet(nlines)*.7;

behaviornames = regexprep(labelfns,'(labels_)|(Labels_)|(labels)|(Labels)','');

combinedzmeanfractime = reshape([zmeanfractime;zmeanfractime_male;zmeanfractime_female],[nexps,nbehaviors*3]);
combinedzmeanfractime_line = reshape([zmeanfractime_line;zmeanfractime_male_line;zmeanfractime_female_line],[nlines,nbehaviors*3]);
combinedzstdfractime_line = reshape([zstdfractime_line;zstdfractime_male_line;zstdfractime_female_line],[nlines,nbehaviors*3]);
combinedbehaviornames = [behaviornames
  cellfun(@(x) [x,'_male'],behaviornames,'UniformOutput',false)
  cellfun(@(x) [x,'_female'],behaviornames,'UniformOutput',false)];
combinedbehaviornames = reshape(combinedbehaviornames,[1,nbehaviors*3]);

hdata = nan(1,nlines);
hold on;
for i = 1:nlines,
  plot((1:3*nbehaviors)-.025+.05*(i-1),combinedzmeanfractime(lineidx==i,:),'.','Color',colors(i,:)*.5+.5);
end
for i = 1:nlines,
  plot(repmat(1:3*nbehaviors,[2,1])-.025+.05*(i-1),...
    bsxfun(@plus,combinedzmeanfractime_line(i,:),[combinedzstdfractime_line(i,:);-combinedzstdfractime_line(i,:)]),...
    '-','Color',colors(i,:));
end
for i = 1:nlines,
  hdata(i) = plot((1:3*nbehaviors)-.025+.05*(i-1),combinedzmeanfractime_line(i,:),'o',...
    'Color',colors(i,:),'MarkerFaceColor',colors(i,:));
end
set(gca,'XTick',1:nbehaviors*3,'XTickLabel',combinedbehaviornames);
rotateticklabel(gca,90)

legend(hdata,unique_linenames,'interpreter','none');
set(gca,'XLim',[0,3*nbehaviors+1]);

ylabel('Z-scored fraction of time');

SaveFigLotsOfWays(hfig,fullfile(outfigdir,'WildtypeBehaviorComparison'));

%% plot by weekend

idxplot = ~cellfun(@isempty,regexp({metadata.date},'^2012060','once'));
meanfractime2 = meanfractime(idxplot,:);
meanfractime_male2 = meanfractime_male(idxplot,:);
meanfractime_female2 = meanfractime_female(idxplot,:);
metadata2 = metadata(idxplot);
nexps2 = nnz(idxplot);

%% z-score the data so that we can plot on the same axes

zmu = nanmean(meanfractime2,1);
zsigma = nanstd(meanfractime2,1,1);
zmeanfractime = bsxfun(@rdivide,bsxfun(@minus,meanfractime2,zmu),zsigma);
zmu_male = nanmean(meanfractime_male2,1);
zsigma_male = nanstd(meanfractime_male2,1,1);
zmeanfractime_male = bsxfun(@rdivide,bsxfun(@minus,meanfractime_male2,zmu_male),zsigma_male);
zmu_female = nanmean(meanfractime_female2,1);
zsigma_female = nanstd(meanfractime_female2,1,1);
zmeanfractime_female = bsxfun(@rdivide,bsxfun(@minus,meanfractime_female2,zmu_female),zsigma_female);

%% compute per-strain data

line_names = {metadata2.line};
[unique_linenames,~,lineidx] = unique(line_names);

nlines = numel(unique_linenames);

meanfractime_line = nan(nlines,nbehaviors);
meanfractime_male_line = nan(nlines,nbehaviors);
meanfractime_female_line = nan(nlines,nbehaviors);
stdfractime_line = nan(nlines,nbehaviors);
stdfractime_male_line = nan(nlines,nbehaviors);
stdfractime_female_line = nan(nlines,nbehaviors);

for linei = 1:nlines,
  meanfractime_line(linei,:) = nanmean(meanfractime2(lineidx==linei,:),1);
  stdfractime_line(linei,:) = nanstd(meanfractime2(lineidx==linei,:),1,1);
  meanfractime_male_line(linei,:) = nanmean(meanfractime_male2(lineidx==linei,:),1);
  stdfractime_male_line(linei,:) = nanstd(meanfractime_male2(lineidx==linei,:),1,1);
  meanfractime_female_line(linei,:) = nanmean(meanfractime_female2(lineidx==linei,:),1);
  stdfractime_female_line(linei,:) = nanstd(meanfractime_female2(lineidx==linei,:),1,1);
end

zmeanfractime_line = bsxfun(@rdivide,bsxfun(@minus,meanfractime_line,zmu),zsigma);
zmeanfractime_male_line = bsxfun(@rdivide,bsxfun(@minus,meanfractime_male_line,zmu_male),zsigma_male);
zmeanfractime_female_line = bsxfun(@rdivide,bsxfun(@minus,meanfractime_female_line,zmu_female),zsigma_female);
zstdfractime_line = bsxfun(@rdivide,stdfractime_line,zsigma);
zstdfractime_male_line = bsxfun(@rdivide,stdfractime_male_line,zsigma_male);
zstdfractime_female_line = bsxfun(@rdivide,stdfractime_female_line,zsigma_female);

%% plot

hfig = 2;
figure(hfig);
clf;
set(hfig,'Units','pixels','Position',[80 94 1716 854]);
%colors = jet(nlines)*.7;
colors = lines(nlines);

behaviornames = regexprep(labelfns,'(labels_)|(Labels_)|(labels)|(Labels)','');

combinedzmeanfractime = reshape([zmeanfractime;zmeanfractime_male;zmeanfractime_female],[nexps2,nbehaviors*3]);
combinedzmeanfractime_line = reshape([zmeanfractime_line;zmeanfractime_male_line;zmeanfractime_female_line],[nlines,nbehaviors*3]);
combinedzstdfractime_line = reshape([zstdfractime_line;zstdfractime_male_line;zstdfractime_female_line],[nlines,nbehaviors*3]);
combinedbehaviornames = [behaviornames
  cellfun(@(x) [x,'_male'],behaviornames,'UniformOutput',false)
  cellfun(@(x) [x,'_female'],behaviornames,'UniformOutput',false)];
combinedbehaviornames = reshape(combinedbehaviornames,[1,nbehaviors*3]);

hdata = nan(1,nlines);
hold on;
for i = 1:nlines,
  plot((1:3*nbehaviors)-.025+.05*(i-1),combinedzmeanfractime(lineidx==i,:),'.','Color',colors(i,:)*.5+.5);
end
for i = 1:nlines,
  plot(repmat(1:3*nbehaviors,[2,1])-.025+.05*(i-1),...
    bsxfun(@plus,combinedzmeanfractime_line(i,:),[combinedzstdfractime_line(i,:);-combinedzstdfractime_line(i,:)]),...
    '-','Color',colors(i,:));
end
for i = 1:nlines,
  hdata(i) = plot((1:3*nbehaviors)-.025+.05*(i-1),combinedzmeanfractime_line(i,:),'o',...
    'Color',colors(i,:),'MarkerFaceColor',colors(i,:));
end
set(gca,'XTick',1:nbehaviors*3,'XTickLabel',combinedbehaviornames);
rotateticklabel(gca,90)

legend(hdata,unique_linenames,'interpreter','none');
set(gca,'XLim',[0,3*nbehaviors+1]);

ylabel('Z-scored fraction of time');

SaveFigLotsOfWays(hfig,fullfile(outfigdir,'WildtypeBehaviorComparison_Weekend2'));

%% plot each behavior separately

[~,idx] = ismember(unique_linenames,linenamecodes(:,1));
print_linenames = linenamecodes(idx,2)';

[~,idx] = ismember(labelfns,behaviorcodes(:,1));
behaviors = behaviorcodes(idx,2)';

hfig = 3;
figure(hfig);
clf;
set(hfig,'Units','pixels','Position',[100,100,1400 1000]);
hax = createsubplots(5,2,[[.05,.05];[.05,.01]]);

hbar = nan(3,nbehaviors);

datacolors = [0,0,0;0,0,.5;.5,0,0];
datanames = {'Both','Male','Female'};
grp1 = repmat(print_linenames(lineidx)',[1,3]);
grp2 = repmat(datanames,[nexps2,1]);
grporder = cell(3,nlines);
for i = 1:nlines,
  for j = 1:3,
    grporder{j,i} = sprintf('%s,%s',print_linenames_order{i},datanames{j});
  end
end
grporder = grporder(:);

for behi = 1:nbehaviors,

  axes(hax(behi));
  
%   v = [meanfractime_line(:,behi),meanfractime_male_line(:,behi),meanfractime_female_line(:,behi)];
%   s = [stdfractime_line(:,behi),stdfractime_male_line(:,behi),stdfractime_female_line(:,behi)];
  p = [meanfractime2(:,behi),meanfractime_male2(:,behi),meanfractime_female2(:,behi)];
%   hbar(:,behi) = bar(1:nlines,v,'grouped');
%   for i = 1:3,
%     set(hbar(i,behi),'FaceColor',datacolors(i,:));
%   end
%   hold on;
%   for i = 1:3,
%     x = mean(get(get(hbar(i,behi),'Children'),'XData'),1);
%     herr(i,behi) = errorbar(x,v(:,i),s(:,i),'linestyle','none','color',[.5,.5,.5]);
%     hexp(i,behi) = plot(x(lineidx),p(:,i),'.','color',[0,.5,.5]);
%   end
  
  
  hb = boxplot(p(:),{grp1(:),grp2(:)},...
  'ColorGroup',grp2(:),...
  'GroupOrder',grporder,...
  'PlotStyle','traditional',...
  'FactorDirection','list',...
  'labelverbosity','minor',...
  'colors',datacolors,...
  'symbol','+');

  hxlabel = findobj(get(hb(1),'Parent'),'Type','text');
  delete(hxlabel);
    
  if behi == 1,
    for i = 1:3,
      hdcurr = findobj(get(hb(1),'Parent'),'Type','line','Color',datacolors(i,:),'Marker','none','LineStyle','-');
      hd(i) = hdcurr(1);
    end
    legend(hd,datanames);
  end
  ylabel(sprintf('Frac. time %s',lower(behaviors{behi})));
  
  if mod(behi,5) == 0,
    set(gca,'XTick',2:3:3*nlines,'XTickLabel',print_linenames_order);
  else
    set(gca,'XTick',2:3:3*nlines,'XTickLabel',{});
  end
  
  %set(gca,'YLim',[0,max(v(:)+s(:))*1.025]);
  %set(gca,'XLim',[.5,nlines+.5]);
  
end


SaveFigLotsOfWays(hfig,fullfile(outfigdir,'WildtypeBehaviorComparison_Weekend2_PerBehaviorBoxPlots'));

%% plot each behavior separately, only plot both statistics

[~,idx] = ismember(unique_linenames,linenamecodes(:,1));
print_linenames = linenamecodes(idx,2)';

% [~,idx] = ismember(labelfns,behaviorcodes(:,1));
% behaviors = behaviorcodes(idx,2)';

ndx2meanfrac = [];
behaviors = {};
for ndx = 1:numel(behaviorcodes(:,1))
  ndx2labelfns = find(cellfun(@(x) strcmp(x,behaviorcodes{ndx,1}),labelfns));
  behaviors{ndx} = behaviorcodes{ndx,2};
  ndx2meanfrac(ndx) = ndx2labelfns;
end
nbehaviors = numel(behaviors);

hfig = 4;
figure(hfig);
clf;
set(hfig,'Units','pixels','Position',[100,100,1400 1000]);
hax = createsubplots(5,2,[[.05,.05];[.05,.01]]);

hbar = nan(3,nbehaviors);

%datacolors = [0,0,0;0,0,.5;.5,0,0];
%datanames = {'Both','Male','Female'};
grp1 = print_linenames(lineidx)';
%grp1 = repmat(print_linenames(lineidx)',[1,3]);
%grp2 = repmat(datanames,[nexps2,1]);
grporder = print_linenames_order;
% grporder = cell(3,nlines);
% for i = 1:nlines,
%   for j = 1:3,
%     grporder{j,i} = sprintf('%s,%s',print_linenames_order{i},datanames{j});
%   end
% end
% grporder = grporder(:);

for behi = 1:nbehaviors,

  axes(hax(behi));
  
%   v = [meanfractime_line(:,behi),meanfractime_male_line(:,behi),meanfractime_female_line(:,behi)];
%   s = [stdfractime_line(:,behi),stdfractime_male_line(:,behi),stdfractime_female_line(:,behi)];
  %p = [meanfractime2(:,behi),meanfractime_male2(:,behi),meanfractime_female2(:,behi)];
  p = [meanfractime2(:,ndx2meanfrac(behi))];
%   hbar(:,behi) = bar(1:nlines,v,'grouped');
%   for i = 1:3,
%     set(hbar(i,behi),'FaceColor',datacolors(i,:));
%   end
%   hold on;
%   for i = 1:3,
%     x = mean(get(get(hbar(i,behi),'Children'),'XData'),1);
%     herr(i,behi) = errorbar(x,v(:,i),s(:,i),'linestyle','none','color',[.5,.5,.5]);
%     hexp(i,behi) = plot(x(lineidx),p(:,i),'.','color',[0,.5,.5]);
%   end
  
  
  hb = boxplot(p(:),grp1(:),...
  'GroupOrder',grporder,...
  'colors',[0,0,0],...
  'symbol','+');

  hxlabel = findobj(get(hb(1),'Parent'),'Type','text');
  delete(hxlabel);
    
%   if behi == 1,
%     for i = 1:3,
%       hdcurr = findobj(get(hb(1),'Parent'),'Type','line','Color',datacolors(i,:),'Marker','none','LineStyle','-');
%       hd(i) = hdcurr(1);
%     end
%     legend(hd,datanames);
%   end
  ylabel(sprintf('Frac. time %s',lower(behaviors{behi})));
  
  if mod(behi,5) == 0,
    set(gca,'XTick',1:nlines,'XTickLabel',print_linenames_order);
  else
    set(gca,'XTick',1:nlines,'XTickLabel',{});
  end
  set(gca,'box','off');
  set(get(gca,'XLabel'),'FontSize',16);
  set(get(gca,'YLabel'),'FontSize',16);
  set(gca,'FontSize',16);
  
  %set(gca,'YLim',[0,max(v(:)+s(:))*1.025]);
  %set(gca,'XLim',[.5,nlines+.5]);
  
end


SaveFigLotsOfWays(hfig,fullfile(outfigdir,'WildtypeBehaviorComparison_Weekend2_PerBehaviorBoxPlots_BothOnly'));

%% Save the Data..

save(fullfile(outfigdir,['WildtypeBehaviorComparison_ResultsData']),'behaviornames','labelfns',...
  'meanfractime2','metadata2','scoresfns');

%% use logistic regression to predict the strain given the behavior classification results for each pair

% reorder line names
print_linenames_perexp = print_linenames(lineidx);
[~,orderedlineidx] = ismember(print_linenames_perexp,print_linenames_order);

% X is nexps x nbehaviors
X = meanfractime2;
y = orderedlineidx';
isbaddata = any(isnan(X),2);
X(isbaddata,:) = [];
y(isbaddata) = [];
errorrate = zeros(nlines,nlines);
errorrate(:) = .5;
for line1 = 1:nlines,
  idx1 = line1 == y;
  for line2 = line1+1:nlines,
    idx2 = line2 == y;
    idxcurr = idx1 | idx2;
    ncurr = nnz(idxcurr);
    Xcurr = X(idxcurr,:);
    % 1 means line1, 0 means line2
    ycurr = idx1(idxcurr);
    % loop over all examples for cross validation
    
    yfitcurr = nan(ncurr,1);
    for i = 1:ncurr,
      
      % logistic regression
      idxtrain = true(ncurr,1); idxtrain(i) = false;
      coeffscurr = glmfit(Xcurr(idxtrain,:),ycurr(idxtrain),'binomial','link','logit');
      yfitcurr(i) = glmval(coeffscurr,Xcurr(i,:),'logit');
      
    end
    
    errorcurr = nnz( (yfitcurr>.5) ~= ycurr );
    errorrate(line1,line2) = errorcurr / ncurr;
    errorrate(line2,line1) = errorrate(line1,line2);
  end
end

hfig = 5;
figure(hfig);
clf;

imagesc(errorrate,[0,.5]);
colormap gray;
set(gca,'XTick',1:nlines,'XTickLabel',print_linenames_order,...
  'YTick',1:nlines,'YTickLabel',print_linenames_order,'Box','off');
axis image;
colorbar;
title('Pairwise error rate');

errorrate0 = errorrate;

%% use regularized logistic regression to predict the strain given the behavior classification results for each pair

% reorder line names
print_linenames_perexp = print_linenames(lineidx);
[~,orderedlineidx] = ismember(print_linenames_perexp,print_linenames_order);
nlambda = 10;

% X is nexps x nbehaviors
X = meanfractime2;
y = orderedlineidx';
isbaddata = any(isnan(X),2);
X(isbaddata,:) = [];
y(isbaddata) = [];
errorrate = zeros(nlines,nlines);
errorrate(:) = .5;
for line1 = 1:nlines,
  idx1 = line1 == y;
  for line2 = line1+1:nlines,
    idx2 = line2 == y;
    idxcurr = idx1 | idx2;
    ncurr = nnz(idxcurr);
    Xcurr = X(idxcurr,:);
    % 1 means line1, 0 means line2
    ycurr = idx1(idxcurr);
    % loop over all examples for cross validation
    
    yfitcurr = nan(ncurr,1);
    for i = 1:ncurr,
      
      % logistic regression
      idxtrain = true(ncurr,1); idxtrain(i) = false;
      [coeffscurr,fitinfo] = lassoglm(Xcurr(idxtrain,:),ycurr(idxtrain),'binomial','link','logit',...
        'NumLambda',nlambda,'CV',ncurr-1);
      j = fitinfo.Index1SE;
      coeffscurr0 = coeffscurr(:,j);
      cnst = fitinfo.Intercept(j);
      coeffscurr1 = [cnst;coeffscurr0];
      yfitcurr(i) = glmval(coeffscurr1,Xcurr(i,:),'logit');      
    end
    
    errorcurr = nnz( (yfitcurr>.5) ~= ycurr );
    errorrate(line1,line2) = errorcurr / ncurr;
    errorrate(line2,line1) = errorrate(line1,line2);
  end
end

hfig = 6;
figure(hfig);
clf;

imagesc(errorrate,[0,.5]);
colormap gray;
set(gca,'XTick',1:nlines,'XTickLabel',print_linenames_order,...
  'YTick',1:nlines,'YTickLabel',print_linenames_order,'Box','off');
xlabel('Wildtype strains');ylabel('Wildtype strains');
axis image;
colorbar;
title('Pairwise error rate');

error_rate = errorrate;
axis image;
box off;
numel(error_rate)
size(error_rate)
[Row,Col]= size(error_rate);
for i = 1:Row
  for j = 1:Col
    if i==j; continue; end
    if errorrate(i,j) ==0; continue; end
    ht(i,j) = text(j,i,sprintf('%.2f',error_rate(i,j)),...
      'FontUnits','pixels','FontWeight','bold','Color',[.7,0,0],...
      'HorizontalAlignment','center','VerticalAlignment','middle',...
      'FontSize',24);
  end
end

SaveFigLotsOfWays(hfig,fullfile(outfigdir,'WildTypeStrainPrediction_PairwiseErrorRate'));

save('WildTypeStrainPredictionData.mat','errorrate','print_linenames_order');

%     'Berlin'    'CantonS MH'    'CantonS UH'    'CantonS GR'    'Dickinson'    'OregonR'
% errorrate =
%     0.5000         0    0.0714    0.0625         0         0
%          0    0.5000    0.0714         0         0         0
%     0.0714    0.0714    0.5000         0         0         0
%     0.0625         0         0    0.5000    0.0625         0
%          0         0         0    0.0625    0.5000    0.0588
%          0         0         0         0    0.0588    0.5000