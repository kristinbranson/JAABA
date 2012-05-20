function [hfigs,errordata,figfilenames,outmatfilename] = ...
  Compute_GroundTruthError_vs_TrainingSetSize(behavior,configfilename,train_expdirs,groundtruth_expdirs,varargin)

ifpcolor = [0,0,.7];
fpcolor = [.3,.3,1];
ifncolor = [.7,0,0];
fncolor = [1,.3,.3];
ierrorcolor = [0,0,0];
errorcolor = [.3,.3,.3];
chancecolor = [.7,.7,.7];
outfigdir = '.';

[hfigs,figpos,ifpcolor,fpcolor,ifncolor,fncolor,...
  ierrorcolor,errorcolor,chancecolor,outfigdir] = myparse(varargin,'hfigs',[],...
  'figpos',[95 550 1246 443],...
  'ifpcolor',ifpcolor,...
  'fpcolor',fpcolor,...
  'ifncolor',ifncolor,...
  'fncolor',fncolor,...
  'ierrorcolor',ierrorcolor,...
  'errorcolor',errorcolor,...
  'chancecolor',chancecolor,...
  'outfigdir',outfigdir);

if ~exist(outfigdir,'dir'),
  mkdir(outfigdir);
end

%%

Jtrain = JLabelData(configfilename,...
  'setstatusfn',@fprintf_wrapper,'clearstatusfn',@() fprintf('Done.\n'));
for i = 1:numel(train_expdirs),
  Jtrain.AddExpDir(train_expdirs{i});
end

%% set groundtruthing data

Jgt = JLabelData(configfilename,...
  'setstatusfn',@fprintf_wrapper,'clearstatusfn',@() fprintf('Done.\n'));
Jgt.SetGTMode(true);
for i = 1:numel(groundtruth_expdirs),
  Jgt.AddExpDir(groundtruth_expdirs{i});
end

%% train on increasing amounts of training data

islabeled = (Jtrain.windowdata.labelidx_new ~= 0) & (Jtrain.windowdata.labelidx_imp);
label_timestamp = Jtrain.GetLabelTimestamps(Jtrain.windowdata.exp(islabeled),Jtrain.windowdata.flies(islabeled,:),Jtrain.windowdata.t(islabeled));
[unique_timestamps] = unique(label_timestamp);
ntrain_bouts = numel(unique_timestamps);

nframestrain = nan(1,ntrain_bouts);

crosserror = nan(4,3,ntrain_bouts);
for bouti = 1:ntrain_bouts,
  
  fprintf('bout %d / %d...\n',bouti,ntrain_bouts);
  
  if bouti == ntrain_bouts,
    timerange = [];
  else
    timerange = [unique_timestamps(1),unique_timestamps(bouti+1)];
  end
  Jtrain.Train(false,timerange);
  
  Jgt.classifier = Jtrain.classifier;
  Jgt.windowdata.isvalidprediction(:) = false;

  Jgt.windowdata.scoreNorm = 1;
  tmp = Jgt.GetGTPerformance();
  
  crosserror(:,:,bouti) = tmp.numbers;
  if bouti == ntrain_bouts,
    nframestrain(bouti) = nnz(label_timestamp);
  else
    nframestrain(bouti) = nnz(label_timestamp < unique_timestamps(bouti+1));
  end
  
  if mod(bouti,10) == 0,
    save('tmp.mat','crosserror','nframestrain','bouti');
  end
  
end

false_positive_important_rate = reshape(crosserror(3,1,:) ./ (crosserror(3,1,:) + crosserror(3,3,:)),[1,ntrain_bouts]);
false_positive_rate = reshape(crosserror(4,1,:) ./ (crosserror(4,1,:) + crosserror(4,3,:)),[1,ntrain_bouts]);
false_negative_important_rate = reshape(crosserror(1,3,:) ./ (crosserror(1,3,:) + crosserror(1,1,:)),[1,ntrain_bouts]);
false_negative_rate = reshape(crosserror(2,3,:) ./ (crosserror(2,3,:) + crosserror(2,1,:)),[1,ntrain_bouts]);

important_error_rate = (false_positive_important_rate+false_negative_important_rate)/2;
unimportant_error_rate = (false_positive_rate+false_negative_rate)/2;

timestamp = max(label_timestamp);
ds = datestr(timestamp,'yyyymmddTHHMMSS');

outmatfilename = fullfile(outfigdir,sprintf('GroundTruthError_TrainingSetSize_%s_%s.mat',behavior,ds));
save(outmatfilename,'-regexp','^false_|crosserror|nframestrain|label_timestamp|^.*important_error_rate');

errordata = struct;
errordata.crosserror = crosserror;
errordata.false_positive_important_rate = false_positive_important_rate;
errordata.label_timestamp = label_timestamp;
errordata.false_negative_important_rate = false_negative_important_rate;
errordata.false_positive_rate = false_positive_rate;
errordata.nframestrain = nframestrain;
errordata.false_negative_rate = false_negative_rate;
errordata.important_error_rate = important_error_rate;
errordata.unimportant_error_rate = unimportant_error_rate;
errordata.ds = ds;

%% plot

if numel(hfigs) >= 1 && ~isnan(hfigs(1)),
  hfig = hfigs(1);
  figure(hfig);
else
  hfig = figure;
  hfigs(1) = hfig;
end
clf;
set(hfig,'Units','pixels','Position',figpos);

xlim = [0,nframestrain(end)+nframestrain(1)];
ylim = [0,1];

legends = {};

plot(nframestrain([1,ntrain_bouts]),[.5,.5],'--','color',chancecolor);
legends{end+1} = 'Chance';

hold on;

plot(nframestrain,false_positive_rate,'.-','color',fpcolor);
legends{end+1} = 'Non-important negatives classified as positives';

plot(nframestrain,false_negative_rate,'.-','color',fncolor);
legends{end+1} = 'Non-important positives classified as negatives';

plot(nframestrain,false_negative_important_rate,'.-','color',ifncolor);
legends{end+1} = 'Important positives classified as negatives';

plot(nframestrain,false_positive_important_rate,'.-','color',ifpcolor);
legends{end+1} = 'Important negatives classified as positives';

plot(nframestrain,unimportant_error_rate,'.-','color',errorcolor);
legends{end+1} = 'Non-important error rate';

plot(nframestrain,important_error_rate,'.-','color',ierrorcolor);
legends{end+1} = 'Important error rate';

axis([xlim,ylim]);

legend(legends);
xlabel('Training set size (frames)');
ylabel('Error rate');

figfilenames{1} = fullfile(outfigdir,sprintf('GroundTruthErrorRateDetails_%s_%s',behavior,ds));
SaveFigLotsOfWays(hfig,figfilenames{1});

%% only plot the error, not false positives and false negatives

if numel(hfigs) >= 2 && ~isnan(hfigs(2)),
  hfig = hfigs(2);
  figure(hfig);
else
  hfig = figure;
  hfigs(2) = hfig;
end
clf;
set(hfig,'Units','pixels','Position',[95 550 1246 443]);

legends = {};

plot(nframestrain([1,ntrain_bouts]),[.5,.5],'--','color',chancecolor);
legends{end+1} = 'Chance';

hold on;

plot(nframestrain,unimportant_error_rate,'.-','color',errorcolor);
legends{end+1} = 'Non-important error rate';

plot(nframestrain,important_error_rate,'.-','color',ierrorcolor);
legends{end+1} = 'Important error rate';

axis([xlim,ylim]);

legend(legends);
xlabel('Training set size (frames)');
ylabel('Error rate');

figfilenames{2} = fullfile(outfigdir,sprintf('GroundTruthErrorRate_%s_%s',behavior,ds));

SaveFigLotsOfWays(hfig,figfilenames{2});