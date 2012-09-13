%% compare behavior detections for different ages

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
outfigdir = '../figures/AgeOut';
if ~exist(outfigdir,'dir'),
  mkdir(outfigdir);
end

[~,computername] = system('hostname');
computername = strtrim(computername);

switch computername,
  
  case 'kabram-ws.janelia.priv',
    rootdatadir0 = '/groups/branson/home/kabram/flyMovies/Age';
  
  otherwise,

  rootdatadir0 = '../experiments/age';
  
end

if ispc,
  addpath C:\Code\FlyBowlAnalysis;
else
  addpath /groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis;
end

%% parameters

conditionnamecodes = {
  'age02days'      '2'
  'age06days'      '6'
  'age13days'      '13'
  'age20days'      '20'
  };

behaviorcodes = {
  'labelsBackup'         'Backup'
  'labelsCrabwalk'       'Crabwalk'
  'labelsTouch'          'Touch'
  'labelsWingGrooming'   'Wing grooming'
  'labels_Chasev7'       'Chase'
  'labels_Jump'          'Jump'
  'labels_Righting'      'Righting'
  'labels_Stops'         'Stop'
  'labels_Walk'          'Walk'
  'labels_pivot_tail'    'Tail pivot turn'
  };

% weekend = '1';


%% data locations

rootdatadir1s = dir(fullfile(rootdatadir0,'age*days'));
rootdatadir1s = cellfun(@(x)fullfile(rootdatadir0,x),{rootdatadir1s([rootdatadir1s.isdir]).name},'UniformOutput',false);

expnames = {};
expdirs = {};

for rooti = 1:numel(rootdatadir1s),
  rootdatadir = rootdatadir1s{rooti};
  % only use weekend 2
%   switch weekend,
%     case 'all',
%       expnamescurr = dir(fullfile(rootdatadir,'*2012*'));
%     case '1',
%       expnamescurr = dir(fullfile(rootdatadir,'*201205*'));
%     case '2',
%       expnamescurr = dir(fullfile(rootdatadir,'*201206*'));
%   end
  % Only use 0830 and 0904 experiments.
  expnamescurr = [];
  expnamescurr = [expnamescurr dir(fullfile(rootdatadir,'*20120830*'))];
  expnamescurr = [expnamescurr dir(fullfile(rootdatadir,'*20120904*'))];
  expnamescurr(~[expnamescurr.isdir]) = [];
  expnamescurr = {expnamescurr.name};
  expnames = [expnames,expnamescurr];
  expdirscurr = cellfun(@(x) fullfile(rootdatadir,x),expnamescurr,'UniformOutput',false);
  expdirs = [expdirs,expdirscurr];
end
nexps = numel(expdirs);

scoresfilestrs = dir(fullfile(expdirs{1},'*cores*.mat'));
scoresfilestrs = {scoresfilestrs.name};
scoresfns = cellfun(@(x) x(1:end-4),scoresfilestrs,'UniformOutput',false);
labelfns = regexprep(scoresfns,'scores','labels','preservecase','once');
nbehaviors = numel(labelfns);


%% compute statistics of behavior detections

meanfractime = nan(nexps,nbehaviors);
meanfractime_male = nan(nexps,nbehaviors);
meanfractime_female = nan(nexps,nbehaviors);

clear trx;

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

%% compute per-condition data

for i = 1:nexps,
  [~,conditionname] = fileparts(fileparts(expdirs{i}));
  metadatacurr = parseExpDir(expdirs{i});
  metadatacurr.condition = conditionname;
  if i == 1,
    metadata = repmat(metadatacurr,[1,nexps]);
  else
    metadata(i) = metadatacurr;
  end
end

condition_names = {metadata.condition};
[unique_conditionnames,~,conditionidx] = unique(condition_names);

nconditions = numel(unique_conditionnames);

meanfractime_condition = nan(nconditions,nbehaviors);
meanfractime_male_condition = nan(nconditions,nbehaviors);
meanfractime_female_condition = nan(nconditions,nbehaviors);
stdfractime_condition = nan(nconditions,nbehaviors);
stdfractime_male_condition = nan(nconditions,nbehaviors);
stdfractime_female_condition = nan(nconditions,nbehaviors);

for conditioni = 1:nconditions,
  meanfractime_condition(conditioni,:) = nanmean(meanfractime(conditionidx==conditioni,:),1);
  stdfractime_condition(conditioni,:) = nanstd(meanfractime(conditionidx==conditioni,:),1,1);
  meanfractime_male_condition(conditioni,:) = nanmean(meanfractime_male(conditionidx==conditioni,:),1);
  stdfractime_male_condition(conditioni,:) = nanstd(meanfractime_male(conditionidx==conditioni,:),1,1);
  meanfractime_female_condition(conditioni,:) = nanmean(meanfractime_female(conditionidx==conditioni,:),1);
  stdfractime_female_condition(conditioni,:) = nanstd(meanfractime_female(conditionidx==conditioni,:),1,1);
end

zmeanfractime_condition = bsxfun(@rdivide,bsxfun(@minus,meanfractime_condition,zmu),zsigma);
zmeanfractime_male_condition = bsxfun(@rdivide,bsxfun(@minus,meanfractime_male_condition,zmu_male),zsigma_male);
zmeanfractime_female_condition = bsxfun(@rdivide,bsxfun(@minus,meanfractime_female_condition,zmu_female),zsigma_female);
zstdfractime_condition = bsxfun(@rdivide,stdfractime_condition,zsigma);
zstdfractime_male_condition = bsxfun(@rdivide,stdfractime_male_condition,zsigma_male);
zstdfractime_female_condition = bsxfun(@rdivide,stdfractime_female_condition,zsigma_female);

%% plot

hfig = 1;
figure(hfig);
clf;
set(hfig,'Units','pixels','Position',[80 94 1716 854]);
tmp = jet(256)*.7;
colors = tmp(round(linspace(1,256,nconditions)),:);
%colors = jet(nconditions)*.7;

behaviornames = regexprep(labelfns,'(labels_)|(Labels_)|(labels)|(Labels)','');

combinedzmeanfractime = reshape([zmeanfractime;zmeanfractime_male;zmeanfractime_female],[nexps,nbehaviors*3]);
combinedzmeanfractime_condition = reshape([zmeanfractime_condition;zmeanfractime_male_condition;zmeanfractime_female_condition],[nconditions,nbehaviors*3]);
combinedzstdfractime_condition = reshape([zstdfractime_condition;zstdfractime_male_condition;zstdfractime_female_condition],[nconditions,nbehaviors*3]);
combinedbehaviornames = [behaviornames
  cellfun(@(x) [x,'_male'],behaviornames,'UniformOutput',false)
  cellfun(@(x) [x,'_female'],behaviornames,'UniformOutput',false)];
combinedbehaviornames = reshape(combinedbehaviornames,[1,nbehaviors*3]);

hdata = nan(1,nconditions);
hold on;
for i = 1:nconditions,
  plot((1:3*nbehaviors)-.025+.05*(i-1),combinedzmeanfractime(conditionidx==i,:),'.','Color',colors(i,:)*.5+.5);
end
for i = 1:nconditions,
  plot(repmat(1:3*nbehaviors,[2,1])-.025+.05*(i-1),...
    bsxfun(@plus,combinedzmeanfractime_condition(i,:),[combinedzstdfractime_condition(i,:);-combinedzstdfractime_condition(i,:)]),...
    '-','Color',colors(i,:));
end
for i = 1:nconditions,
  hdata(i) = plot((1:3*nbehaviors)-.025+.05*(i-1),combinedzmeanfractime_condition(i,:),'o',...
    'Color',colors(i,:),'MarkerFaceColor',colors(i,:));
end
set(gca,'XTick',1:nbehaviors*3,'XTickLabel',combinedbehaviornames);
rotateticklabel(gca,90);

legend(hdata,unique_conditionnames,'interpreter','none');
set(gca,'XLim',[0,3*nbehaviors+1]);

ylabel('Z-scored fraction of time');

SaveFigLotsOfWays(hfig,fullfile(outfigdir,['AgeBehaviorComparison_wk',weekend]));

%% plot each behavior separately

[~,idx] = ismember(unique_conditionnames,conditionnamecodes(:,1));
print_conditionnames = conditionnamecodes(idx,2)';
[~,idx] = ismember(behaviorcodes(:,1),labelfns);
behaviors = labelfns(idx)';

hfig = 3;
figure(hfig);
clf;
set(hfig,'Units','pixels','Position',[100,100,1400 1000]);
hax = createsubplots(5,2,[[.05,.05];[.05,.01]]);

hbar = nan(3,nbehaviors);

datacolors = [0,0,0;0,0,.5;.5,0,0];
datanames = {'Both','Male','Female'};
grp1 = repmat(print_conditionnames(conditionidx)',[1,3]);
grp2 = repmat(datanames,[nexps,1]);

grporder = cell(3,nconditions);
for i = 1:nconditions,
  for j = 1:3,
    grporder{j,i} = sprintf('%s,%s',print_conditionnames{i},datanames{j});
  end
end
grporder = grporder(:);

for behi = 1:nbehaviors,

  axes(hax(behi));
  
%   v = [meanfractime_condition(:,behi),meanfractime_male_condition(:,behi),meanfractime_female_condition(:,behi)];
%   s = [stdfractime_condition(:,behi),stdfractime_male_condition(:,behi),stdfractime_female_condition(:,behi)];
  p = [meanfractime(:,behi),meanfractime_male(:,behi),meanfractime_female(:,behi)];
%   hbar(:,behi) = bar(1:nconditions,v,'grouped');
%   for i = 1:3,
%     set(hbar(i,behi),'FaceColor',datacolors(i,:));
%   end
%   hold on;
%   for i = 1:3,
%     x = mean(get(get(hbar(i,behi),'Children'),'XData'),1);
%     herr(i,behi) = errorbar(x,v(:,i),s(:,i),'conditionstyle','none','color',[.5,.5,.5]);
%     hexp(i,behi) = plot(x(conditionidx),p(:,i),'.','color',[0,.5,.5]);
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
    set(gca,'XTick',2:3:3*nconditions,'XTickLabel',print_conditionnames);
    xlabel('Age');
  else
    set(gca,'XTick',2:3:3*nconditions,'XTickLabel',{});
  end
  
  %set(gca,'YLim',[0,max(v(:)+s(:))*1.025]);
  %set(gca,'XLim',[.5,nconditions+.5]);
  
end


SaveFigLotsOfWays(hfig,fullfile(outfigdir,['AgeBehaviorComparison_PerBehaviorBoxPlots_wk',weekend]));


%% plot each behavior separately, only plot both statistics

[~,idx] = ismember(unique_conditionnames,conditionnamecodes(:,1));
print_conditionnames = conditionnamecodes(idx,2)';

[~,idx] = ismember(labelfns,behaviorcodes(:,1));
behaviors = behaviorcodes(idx,2)';

hfig = 4;
figure(hfig);
clf;
set(hfig,'Units','pixels','Position',[100,100,1400 1000]);
hax = createsubplots(5,2,[[.05,.05];[.05,.01]]);

hbar = nan(3,nbehaviors);

%datacolors = [0,0,0;0,0,.5;.5,0,0];
%datanames = {'Both','Male','Female'};
grp1 = print_conditionnames(conditionidx)';
%grp1 = repmat(print_conditionnames(conditionidx)',[1,3]);
%grp2 = repmat(datanames,[nexps2,1]);
grporder = print_conditionnames;
% grporder = cell(3,nconditions);
% for i = 1:nconditions,
%   for j = 1:3,
%     grporder{j,i} = sprintf('%s,%s',print_conditionnames_order{i},datanames{j});
%   end
% end
% grporder = grporder(:);

for behi = 1:nbehaviors,

  axes(hax(behi));
  
%   v = [meanfractime_condition(:,behi),meanfractime_male_condition(:,behi),meanfractime_female_condition(:,behi)];
%   s = [stdfractime_condition(:,behi),stdfractime_male_condition(:,behi),stdfractime_female_condition(:,behi)];
  %p = [meanfractime2(:,behi),meanfractime_male2(:,behi),meanfractime_female2(:,behi)];
  p = [meanfractime(:,behi)];
%   hbar(:,behi) = bar(1:nconditions,v,'grouped');
%   for i = 1:3,
%     set(hbar(i,behi),'FaceColor',datacolors(i,:));
%   end
%   hold on;
%   for i = 1:3,
%     x = mean(get(get(hbar(i,behi),'Children'),'XData'),1);
%     herr(i,behi) = errorbar(x,v(:,i),s(:,i),'linestyle','none','color',[.5,.5,.5]);
%     hexp(i,behi) = plot(x(conditionidx),p(:,i),'.','color',[0,.5,.5]);
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
    set(gca,'XTick',1:nconditions,'XTickLabel',print_conditionnames);
  else
    set(gca,'XTick',1:nconditions,'XTickLabel',{});
  end
  
  %set(gca,'YLim',[0,max(v(:)+s(:))*1.025]);
  %set(gca,'XLim',[.5,nconditions+.5]);
  
end


SaveFigLotsOfWays(hfig,fullfile(outfigdir,['AgeBehaviorComparison_PerBehaviorBoxPlots_BothOnly_wk',weekend]));


%% use logistic regression to predict the strain given the behavior classification results for each pair

% reorder condition names
print_conditionnames_perexp = print_conditionnames(conditionidx);
[~,orderedconditionidx] = ismember(print_conditionnames_perexp,print_conditionnames);

[~,idx] = ismember(labelfns,behaviorcodes(:,1));
behaviors = behaviorcodes(idx,2)';

% X is nexps x nbehaviors
X = meanfractime(:,idx);
y = orderedconditionidx';
isbaddata = any(isnan(X),2);
X(isbaddata,:) = [];
y(isbaddata) = [];
errorrate = zeros(nconditions,nconditions);
errorrate(:) = .5;
for condition1 = 1:nconditions,
  idx1 = condition1 == y;
  for condition2 = condition1+1:nconditions,
    idx2 = condition2 == y;
    idxcurr = idx1 | idx2;
    ncurr = nnz(idxcurr);
    Xcurr = X(idxcurr,:);
    % 1 means condition1, 0 means condition2
    ycurr = idx1(idxcurr);
    % loop over all examples for cross validation
    
    yfitcurr = nan(ncurr,1);
    
    for i = 1:ncurr,
            
      idxtrain = true(ncurr,1); idxtrain(i) = false;
      
      % equal weight for positives and negatives
      weightscurr = nan(nnz(idxtrain),1);
      n0 = nnz(ycurr(idxtrain)==0);
      n1 = nnz(ycurr(idxtrain)==1);
      weightscurr(ycurr(idxtrain)==0) = 1 / n0;
      weightscurr(ycurr(idxtrain)==1) = 1 / n1;

      
      % logistic regression
      coeffscurr = glmfit(Xcurr(idxtrain,:),ycurr(idxtrain),'binomial','link','logit','weights',weightscurr);
      yfitcurr(i) = glmval(coeffscurr,Xcurr(i,:),'logit');
      
    end
    
    errorcurr = nnz( (yfitcurr>.5) ~= ycurr );
    errorrate(condition1,condition2) = errorcurr / ncurr;
    errorrate(condition2,condition1) = errorrate(condition1,condition2);
  end
end

hfig = 5;
figure(hfig);
clf;

imagesc(errorrate,[0,.5]);
colormap gray;
set(gca,'XTick',1:nconditions,'XTickLabel',print_conditionnames,...
  'YTick',1:nconditions,'YTickLabel',print_conditionnames,'Box','off');
axis image;
colorbar;
title('Pairwise error rate');

errorrate0 = errorrate;

%% use regularized logistic regression to predict the strain given the behavior classification results for each pair

% reorder condition names
print_conditionnames_perexp = print_conditionnames(conditionidx);
[~,orderedconditionidx] = ismember(print_conditionnames_perexp,print_conditionnames);
nlambda = 10;

[idx,~] = ismember(labelfns,behaviorcodes(:,1));

% X is nexps x nbehaviors
X = meanfractime(:,idx);

y = orderedconditionidx';
isbaddata = any(isnan(X),2);
X(isbaddata,:) = [];
y(isbaddata) = [];
errorrate = zeros(nconditions,nconditions);
errorrate(:) = .5;
for condition1 = 1:nconditions,
  idx1 = condition1 == y;
  for condition2 = condition1+1:nconditions,
    
    fprintf('condition1 = %d, condition2 = %d\n',condition1,condition2);
    
    idx2 = condition2 == y;
    idxcurr = idx1 | idx2;
    ncurr = nnz(idxcurr);
    Xcurr = X(idxcurr,:);
    % 1 means condition1, 0 means condition2
    ycurr = idx1(idxcurr);
    % loop over all examples for cross validation
    
    yfitcurr = nan(ncurr,1);
    for i = 1:ncurr,
      
      idxtrain = true(ncurr,1); idxtrain(i) = false;
      
      % equal weight for positives and negatives
      weightscurr = nan(nnz(idxtrain),1);
      n0 = nnz(ycurr(idxtrain)==0);
      n1 = nnz(ycurr(idxtrain)==1);
      weightscurr(ycurr(idxtrain)==0) = 1 / n0;
      weightscurr(ycurr(idxtrain)==1) = 1 / n1;

      
      % logistic regression
      [coeffscurr,fitinfo] = lassoglm(Xcurr(idxtrain,:),ycurr(idxtrain),'binomial','link','logit',...
        'NumLambda',nlambda,'CV',ncurr-1,'Weights',weightscurr);
      j = fitinfo.Index1SE;
      coeffscurr0 = coeffscurr(:,j);
      cnst = fitinfo.Intercept(j);
      coeffscurr1 = [cnst;coeffscurr0];
      yfitcurr(i) = glmval(coeffscurr1,Xcurr(i,:),'logit');      
    end
    
    errorcurr = nnz( (yfitcurr>.5) ~= ycurr );
    errorrate(condition1,condition2) = errorcurr / ncurr;
    errorrate(condition2,condition1) = errorrate(condition1,condition2);
  end
end

hfig = 6;
figure(hfig);
clf;

imagesc(errorrate,[0,.5]);
colormap gray;
set(gca,'XTick',1:nconditions,'XTickLabel',print_conditionnames,...
  'YTick',1:nconditions,'YTickLabel',print_conditionnames,'Box','off');
axis image;
colorbar;
title('Pairwise error rate');

SaveFigLotsOfWays(hfig,fullfile(outfigdir,'AgePrediction_PairwiseErrorRate'));

save(['AgePredictionData.mat'],'errorrate','print_conditionnames');

% all:
%     2         4         5         9         12
%     0.5000         0    1.0000    0.1667         0
%          0    0.5000    0.2500    0.1667         0
%     1.0000    0.2500    0.5000    0.0833         0
%     0.1667    0.1667    0.0833    0.5000    0.1667
%          0         0         0    0.1667    0.5000

% wk 2:
%     4         9         12

% Aug 30 and Sep 4 2012 data with all the 13 behaviors.
%        2        6         13      20
%     0.5000    0.0833    0.0833         0
%     0.0833    0.5000    0.1250    0.1667
%     0.0833    0.1250    0.5000    0.1667
%          0    0.1667    0.1667    0.5000

% Aug 30 and Sep 4 2012 data with 10 old behaviors.
%        2        6         13      20
%     0.5000    0.1250    0.0833         0
%     0.1250    0.5000    0.1250    0.1667
%     0.0833    0.1250    0.5000    0.0417
%          0    0.1667    0.0417    0.5000

%% Logistic regression parameters -- Mayank Sep 12.

% reorder condition names
print_conditionnames_perexp = print_conditionnames(conditionidx);
[~,orderedconditionidx] = ismember(print_conditionnames_perexp,print_conditionnames);
nlambda = 10;

[idx,~] = ismember(labelfns,behaviorcodes(:,1));

% X is nexps x nbehaviors
X = meanfractime(:,idx);
y = orderedconditionidx';
isbaddata = any(isnan(X),2);
X(isbaddata,:) = [];
y(isbaddata) = [];
coeffsall = {};
fitinfoall = {};
topbehavior = {};
for condition1 = 1:nconditions,
  idx1 = condition1 == y;
  for condition2 = condition1+1:nconditions,
    
    fprintf('condition1 = %d, condition2 = %d\n',condition1,condition2);
    
    idx2 = condition2 == y;
    idxcurr = idx1 | idx2;
    ncurr = nnz(idxcurr);
    Xcurr = X(idxcurr,:);
    % 1 means condition1, 0 means condition2
    ycurr = idx1(idxcurr);
    % loop over all examples for cross validation
    
    yfitcurr = nan(ncurr,1);
    idxtrain = true(ncurr,1); 
      
    % equal weight for positives and negatives
    weightscurr = nan(nnz(idxtrain),1);
    n0 = nnz(ycurr(idxtrain)==0);
    n1 = nnz(ycurr(idxtrain)==1);
    weightscurr(ycurr(idxtrain)==0) = 1 / n0;
    weightscurr(ycurr(idxtrain)==1) = 1 / n1;
    
    
    % logistic regression
    [coeffscurr,fitinfo] = lassoglm(Xcurr(idxtrain,:),ycurr(idxtrain),'binomial','link','logit',...
      'NumLambda',nlambda,'CV',ncurr-1,'Weights',weightscurr);
    j = fitinfo.Index1SE;
    coeffsall{condition1,condition2} = coeffscurr(:,j);
    [~,tb] = max(abs(coeffscurr(:,j)));
    topbehavior{condition1,condition2} = behaviors{tb};
    fitinfoall{condition1,condition2} = fitinfo;
  end
end



%% plot first two principal components

X = [meanfractime];

% z-score
mu = mean(X,1);
sigma = std(X,1,1);
X = bsxfun(@rdivide,bsxfun(@minus,X,mu),sigma);

[coeff,scores] = princomp(X);

hfig = 7;
figure(hfig);
clf;

hold on;
tmp = jet(256);
colors = tmp(round(linspace(1,256,nconditions)),:)*.7;

h = nan(1,nconditions);
for i = 1:nconditions,
  idxcurr = orderedconditionidx == i;
  h(i) = plot(scores(idxcurr,1),scores(idxcurr,2),'o','color',colors(i,:),'markerfacecolor',colors(i,:));
end

legend(h,cellfun(@(x) [x,' days'],print_conditionnames,'UniformOutput',false));
axisalmosttight;
xlabel('Principal component 1');
ylabel('Principal component 2');
set(gca,'Box','off');

SaveFigLotsOfWays(hfig,fullfile(outfigdir,['AgePrincipalComponents_wk',weekend]));

%% plot linear discriminant projection

X = [meanfractime];


% z-score
mu = mean(X,1);
sigma = std(X,1,1);
X = bsxfun(@rdivide,bsxfun(@minus,X,mu),sigma);

Sb = zeros(size(X,2));
Sw = zeros(size(X,2));
mutotal = mean(X,1);
for i = 1:nconditions,
  idxcurr = orderedconditionidx == i;
  mucurr = mean(X(idxcurr,:),1);
  xcurr = bsxfun(@minus,X(idxcurr,:),mucurr);
  Sw = Sw + xcurr'*xcurr;
  Sb = Sb + nnz(idxcurr)*(mucurr-mutotal)'*(mucurr-mutotal);
end
[v,d] = eigs(Sb,Sw,2);
proj = X*v;

hfig = 8;
figure(hfig);
clf;

hold on;
tmp = jet(256);
colors = tmp(round(linspace(1,256,nconditions)),:)*.7;

h = nan(1,nconditions);
for i = 1:nconditions,
  idxcurr = orderedconditionidx == i;
  h(i) = plot(proj(idxcurr,1),proj(idxcurr,2),'o','color',colors(i,:),'markerfacecolor',colors(i,:));
end

legend(h,cellfun(@(x) [x,' days'],print_conditionnames,'UniformOutput',false));
axisalmosttight;
xlabel('Discriminant 1');
ylabel('Discriminant 2');
set(gca,'Box','off');

SaveFigLotsOfWays(hfig,fullfile(outfigdir,['AgeLDA_wk',weekend]));