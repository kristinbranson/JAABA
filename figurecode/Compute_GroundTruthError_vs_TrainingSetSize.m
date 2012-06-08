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
fracframes = logspace(-2.25,0,20);

[hfigs,figpos,ifpcolor,fpcolor,ifncolor,fncolor,...
  ierrorcolor,errorcolor,chancecolor,outfigdir,fracframes] = myparse(varargin,'hfigs',[],...
  'figpos',[95 550 1246 443],...
  'ifpcolor',ifpcolor,...
  'fpcolor',fpcolor,...
  'ifncolor',ifncolor,...
  'fncolor',fncolor,...
  'ierrorcolor',ierrorcolor,...
  'errorcolor',errorcolor,...
  'chancecolor',chancecolor,...
  'outfigdir',outfigdir,...
  'fracframes',fracframes);
fracframes = sort(fracframes);

if ~exist(outfigdir,'dir'),
  mkdir(outfigdir);
end

%%

Jtrain = JLabelData(configfilename,...
  'setstatusfn',@fprintf_wrapper,'clearstatusfn',@() fprintf('Done.\n'));
for i = 1:numel(train_expdirs),
  Jtrain.AddExpDir(train_expdirs{i});
end

Jtrain.TrimWindowData(true);

%% set groundtruthing data

Jgt = JLabelData(configfilename,...
  'setstatusfn',@fprintf_wrapper,'clearstatusfn',@() fprintf('Done.\n'));
Jgt.SetGTMode(true);
for i = 1:numel(groundtruth_expdirs),
  Jgt.AddExpDir(groundtruth_expdirs{i});
end
%Jgt.TrimWindowData(true);

%% train on increasing amounts of training data

% find bouts of training data
islabeled = (Jtrain.windowdata.labelidx_new ~= 0) & (Jtrain.windowdata.labelidx_imp);
label_timestamp = Jtrain.GetLabelTimestamps(Jtrain.windowdata.exp(islabeled),Jtrain.windowdata.flies(islabeled,:),Jtrain.windowdata.t(islabeled));
[unique_timestamps,~,timestampidx] = unique(label_timestamp);
if min(unique_timestamps) == 0,
  error('timestamp = 0!!\n');
end

ntrain_bouts = numel(unique_timestamps);

% check if there is a unique time stamp for each bout
idxlabeled = find(islabeled);
for i = 1:ntrain_bouts,
  idxcurr = idxlabeled(timestampidx==i);
  if numel(unique(Jtrain.windowdata.exp(idxcurr))) > 1 || ...
      numel(unique(Jtrain.windowdata.flies(idxcurr))) > 1,
    error('Same timestamp for different experiments of flies');
  end
%   ts = sort(Jtrain.windowdata.t(idxcurr));
%   if any(diff(ts) > 1),
%     warning('Multiple bouts for same timestamp');
%   end
end

% choose limits on frames based on fracframes
nframes_train_total = nnz(islabeled);
frame_lims = fracframes*nframes_train_total;
nframes_per_bout = nan(1,ntrain_bouts);
for i = 1:ntrain_bouts,
  nframes_per_bout(i) = nnz(label_timestamp == unique_timestamps(i));
end
nframestrain = cumsum(nframes_per_bout);

% choose the closest bout end to frame limits
bouti_lims = [];
for i = 1:numel(fracframes),  
  [~,bouti] = min(abs(frame_lims(i)-nframestrain));
  if i > 1 && bouti <= lastbouti && frame_lims(i) >= nframestrain(bouti) && bouti < ntrain_bouts,
    bouti = bouti + 1;
  end
  if i > 1 && bouti <= lastbouti,
    continue;
  end
  bouti_lims(end+1) = bouti; %#ok<AGROW>
  lastbouti = bouti;
end

niters = numel(bouti_lims);


nframestrain = nan(1,niters);

crosserror = nan(4,3,niters);
for i = 1:niters,

  bouti = bouti_lims(i);
  if bouti == ntrain_bouts,
    nframestrain(i) = nnz(label_timestamp);
  else
    nframestrain(i) = nnz(label_timestamp < unique_timestamps(bouti+1));
  end
  fprintf('Iteration %d / %d, training on %d / %d frames...\n',i,niters,nframestrain(i),nframes_train_total);
  
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
  
  crosserror(:,:,i) = tmp.numbers;
  
  if mod(i,10) == 0,
    save('tmp.mat','crosserror','nframestrain','bouti');
  end
  
end

false_positive_important_rate = reshape(crosserror(3,1,:) ./ (crosserror(3,1,:) + crosserror(3,3,:)),[1,niters]);
false_positive_rate = reshape(crosserror(4,1,:) ./ (crosserror(4,1,:) + crosserror(4,3,:)),[1,niters]);
false_negative_important_rate = reshape(crosserror(1,3,:) ./ (crosserror(1,3,:) + crosserror(1,1,:)),[1,niters]);
false_negative_rate = reshape(crosserror(2,3,:) ./ (crosserror(2,3,:) + crosserror(2,1,:)),[1,niters]);

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

plot(nframestrain([1,niters]),[.5,.5],'--','color',chancecolor);
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

figfilenames{1} = fullfile(outfigdir,sprintf('GroundTruthErrorRateDetails1_%s_%s',behavior,ds));
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

plot(nframestrain([1,niters]),[.5,.5],'--','color',chancecolor);
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

figfilenames{2} = fullfile(outfigdir,sprintf('GroundTruthErrorRate1_%s_%s',behavior,ds));

SaveFigLotsOfWays(hfig,figfilenames{2});