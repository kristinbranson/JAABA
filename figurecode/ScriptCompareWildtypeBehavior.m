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
rootdatadir = '../experiments/wildtype/results';

if ispc,
  addpath C:\Code\FlyBowlAnalysis;
else
  addpath /groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis;
end

%% data locations

expnames = dir(fullfile(rootdatadir,'*2012*'));
expnames(~[expnames.isdir]) = [];
expnames = {expnames.name};
expdirs = cellfun(@(x) fullfile(rootdatadir,x),expnames,'UniformOutput',false);
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

for expi = expi:nexps,
  
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