% script for trying out random fern code

%% set up path

[~,computername] = system('hostname');
computername = strtrim(computername);

switch computername,
  
  case 'bransonk-lw1',

    JCtrax_path = 'E:\Code\JCtrax';
    FlyBowlAnalysis_path = 'E:\Code\FlyBowlAnalysis';
    rootdatadir = 'E:\Code\Jdetect\larva\fly_data\TrainingData';
    settingsdir = 'E:\Code\FlyBowlAnalysis\settings';
    mRMR_path = 'E:\Code\mRMR_0.9_compiled';

  case 'bransonk-lw2',

    JCtrax_path = 'C:\Code\JCtrax';
    FlyBowlAnalysis_path = 'C:\Code\FlyBowlAnalysis';
    rootdatadir = 'C:\Code\Jdetect\larva\fly_data\TrainingData';
    settingsdir = 'C:\Code\FlyBowlAnalysis\settings';
    mRMR_path = 'C:\Code\mrmr';

  case 'bransonk-desktop',

    JCtrax_path = '/groups/branson/home/bransonk/tracking/code/JCtrax';
    FlyBowlAnalysis_path = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis';
    rootdatadir = '/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/larva/fly_data/TrainingData/';
    settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
    mRMR_path = '/groups/branson/home/bransonk/behavioranalysis/code/mRMR_0.9_compiled';
    
  otherwise
    
    warning('Unknown computer %s. Paths may not be setup correctly',computername);
    JCtrax_path = '/groups/branson/home/bransonk/tracking/code/JCtrax';
    FlyBowlAnalysis_path = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis';
    rootdatadir = '/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/larva/fly_data/TrainingData/';
    settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
    mRMR_path = '/groups/branson/home/bransonk/behavioranalysis/code/mRMR_0.9_compiled';

end

addpath(genpath(fullfile(JCtrax_path,'pdollar_toolbox')));
addpath(fullfile(JCtrax_path,'misc'));
addpath(fullfile(JCtrax_path,'filehandling'));
addpath(FlyBowlAnalysis_path);
addpath(genpath(mRMR_path));


%% data, parameter locations

labelmatnames = {
  'perframe_labeledsharpturns_movie01_pBDPGAL4U_TrpA_Rig1Plate10BowlA_20110202T105734_fly1_11.mat'
  'perframe_labeledsharpturns_movie01_pBDPGAL4U_TrpA_Rig1Plate10BowlA_20110202T105734_fly2_08.mat'
  'perframe_labeledsharpturns_movie02_pBDPGAL4U_TrpA_Rig2Plate14BowlA_20110202T110111_fly1_01.mat'
  'perframe_labeledsharpturns_movie02_pBDPGAL4U_TrpA_Rig2Plate14BowlA_20110202T110111_fly2_15.mat'
  'perframe_labeledsharpturns_movie03_pBDPGAL4U_TrpA_Rig2Plate14BowlB_20110202T110116_fly1_12.mat'
  'perframe_labeledsharpturns_movie03_pBDPGAL4U_TrpA_Rig2Plate14BowlB_20110202T110116_fly2_03.mat'
  'perframe_labeledsharpturns_movie07_GMR_14G03_AE_01_TrpA_Rig2Plate14BowlA_20110202T102516_fly2_17.mat'
  'perframe_labeledsharpturns_movie08_GMR_14H07_AE_01_TrpA_Rig1Plate10BowlB_20110202T084231_fly1_09.mat'
  'perframe_labeledsharpturns_movie09_GMR_16E02_AE_01_TrpA_Rig1Plate10BowlC_20110202T140143_fly1_09.mat'
  };
nexps = numel(labelmatnames);

paramsfilename = 'SharpturnPerFrameParams.xml';

%% parameters

% parameters files to use when evaluating fern parameters & associated mat
% files for caching results
windowfeaturesparamsfile_lots = 'WindowFeatures_Lots.xml';
windowfeaturesmatfile_lots = 'WindowFeatures20110608.mat';
windowfeaturesparamsfile_few = 'SharpturnPerFrameParams.xml';
windowfeaturesmatfile_few = 'SharpturnWindowFeatures20110616.mat';
paramsfilenames = {windowfeaturesparamsfile_lots,windowfeaturesparamsfile_few};
windowfeaturesmatfiles = {windowsfeaturesmatfile_lots,windowfeaturesmatfile_few};

% nferns to test
n_nferns_tests_lots = 100;
nferns_test_lots = round(logspace(1,4,n_nferns_tests_lots));
n_nferns_tests_few = 5;
nferns_test_few = round(logspace(1,4,n_nferns_tests_lots));

% fern depths to test
fern_depth_tests_lots = 3:15;
n_fern_depth_tests_lots = numel(fern_depth_tests_lots);
fern_depth_tests_few = [5,10,15];
n_fern_depth_tests_few = numel(fern_depth_tests_few);

% one setting of fern parameters
nferns = 100;
fern_depth = 10;
threshold_range = [-4,4];
maxdiff_logprob0 = 0;
mindiff_logprob1 = 0;
fern_params = struct('S',fern_depth,'M',nferns,'thrr',threshold_range);

% whether to recompute window features
docomputewindowfeatures = false;

%% load data

X = cell(1,nexps);
Y = cell(1,nexps);
feature_names = {};

for expi = 1:nexps,

  labelfilename = fullfile(rootdatadir,labelmatnames{expi});
  if ~exist(labelfilename,'file'),
    warning('File %s does not exist, skipping',labelfilename);
    continue;
  end
  labeldata = load(labelfilename);
  for i = 1:numel(fns),
    fn = fns{i};
    j = find(strcmp(labeldata.perframefns,fn));
    if isempty(j),
      error('Unknown perframe field %s',fn);
    end
    [x_curr,feature_names_curr] = ...
      ComputeWindowFeatures(labeldata.perframedata{j},cellparams.(fn){:});
    if expi == 1,
      feature_names = [feature_names,cellfun(@(s) [{fn},s],feature_names_curr,'UniformOutput',false)];
    end
    
    % crop out labeled part of video
    x_curr = x_curr(:,labeldata.trk_labelstart:labeldata.trk_labelend-1);

    % store
    X{expi} = [X{expi};x_curr];
  end
  
  % store labels
  Y{expi} = labeldata.label;
  
end
nfeatures = numel(feature_names);

%% display names for stats

statnames = cell(1,nfeatures);
for i = 1:nfeatures,
  fn = '';
  for j = 1:numel(feature_names{i}),
    if j == 1,
    elseif mod(j,2) == 0,
      fn = [fn,', '];
    else
      fn = [fn,'='];
    end
    if ischar(feature_names{i}{j}),
      fn = [fn,feature_names{i}{j}];
    else
      fn = [fn,num2str(feature_names{i}{j})];
    end
  end
  statnames{i} = fn;
end

%% combine exps

X_combine = X{1};
Y_combine = Y{1};
expi_combine = ones(1,size(X{1},2));
for expi = 2:nexps,
  if isempty(X{expi}),
    continue;
  end
  X_combine = [X_combine,nan([nfeatures,1]),X{expi}];
  Y_combine = [Y_combine,nan,Y{expi}];
  expi_combine = [expi_combine,nan,expi+zeros(1,size(X{expi},2))];
end
N = numel(Y_combine);

%% remove unknowns

idx = ~isnan(Y_combine);
X_known = X_combine(:,idx);
Y_known = Y_combine(idx);
expi_known = expi_combine(idx);

%% train & test on everything
X_train = X_known;
Y_train = Y_known;

X_test = X_known;
Y_test = Y_known;

%% z-score data
mu = nanmean(X_train,2);
sig = nanstd(X_train,1,2);
X_train_zscore = bsxfun(@rdivide,bsxfun(@minus,X_train,mu),sig);
X_test_zscore = bsxfun(@rdivide,bsxfun(@minus,X_test,mu),sig);

%% train

nferns = 100;
fern_depth = 10;
threshold_range = [-4,4];
maxdiff_logprob0 = 0;
mindiff_logprob1 = 0;

params = struct('S',fern_depth,'M',nferns,'thrr',threshold_range);
fern_model = fernsClfTrain( X_train_zscore',Y_train+1,params);

[Y_test_pred2,logprob] = fernsClfApply(X_test_zscore',fern_model);
Y_test_pred2 = Y_test_pred2-1;
Y_test_pred = 2+zeros(1,size(logprob,1));
dlogprob = logprob(:,2) - logprob(:,1);
Y_test_pred(dlogprob>=mindiff_logprob1) = 1;
Y_test_pred(dlogprob<=maxdiff_logprob0) = 0;

%% precision-recall curves on training set

nsamples = 100;
dlogprob_lim = [min(dlogprob),max(dlogprob)];
dlogprob_samples = SelectHistEdges(nsamples-1,dlogprob_lim,'logabs');

precision0 = nan(1,nsamples);
precision1 = nan(1,nsamples);
recall0 = nan(1,nsamples);
recall1 = nan(1,nsamples);

n0 = nnz(Y_test==0);
n1 = nnz(Y_test==1);
for i = 1:nsamples,
  thresh = dlogprob_samples(i);
  idx1 = dlogprob>=thresh;
  idx0 = dlogprob<=thresh;
  n0curr = nnz(idx0' & Y_test==0);
  n1curr = nnz(idx1' & Y_test==1);
  recall0(i) = n0curr/n0;
  recall1(i) = n1curr/n1;
  precision0(i) = n0curr/nnz(idx0);
  precision1(i) = n1curr/nnz(idx1);
end

handles = struct;
handles.fig = 3;
figure(handles.fig);
clf;
handles.pr0 = plot(recall0,precision0,'-','color',[0,.5,0],'linewidth',3);
hold on;
handles.pr1 = plot(recall1,precision1,'-','color',[0,0,0],'linewidth',3);
scatter(recall0,precision0,[],dlogprob_samples,'.');
scatter(recall1,precision1,[],dlogprob_samples,'.');
legend([handles.pr0,handles.pr1],{'none','turn'});
xlabel('Recall (|ypred==c & ytrue==c| / |ytrue==c|)');
ylabel('Precision (|ypred==c & ytrue==c| / |ypred==c|)');
colorbar;
colormap(jet(nsamples)*.7);
title('Precision-Recall, training set');

%% train & test on different sets

nferns = 100;
fern_depth = 10;
threshold_range = [-4,4];
params = struct('S',fern_depth,'M',nferns,'thrr',threshold_range);

% all outputs
Y_test = Y_known;
dlogprob = nan(numel(Y_known),1);

for expi_test = 1:nexps,

  idx_train = expi_known~=expi_test;
  idx_test = ~idx_train;
  X_train = X_known(:,idx_train);
  Y_train = Y_known(idx_train);
  
  X_test = X_known(:,idx_test);

  %% z-score data
  mu = nanmean(X_train,2);
  sig = nanstd(X_train,1,2);
  X_train_zscore = bsxfun(@rdivide,bsxfun(@minus,X_train,mu),sig);
  X_test_zscore = bsxfun(@rdivide,bsxfun(@minus,X_test,mu),sig);

  %% train
  fern_model = fernsClfTrain( X_train_zscore',Y_train+1,params);

  [~,logprob] = fernsClfApply(X_test_zscore',fern_model);
  dlogprob(idx_test) = logprob(:,2) - logprob(:,1);

end

X_test = X_known;
  
%% precision-recall curves, cross-validation

nsamples = 100;
dlogprob_lim = [min(dlogprob),max(dlogprob)];
dlogprob_samples = SelectHistEdges(nsamples-1,dlogprob_lim,'logabs');

precision0 = nan(1,nsamples);
precision1 = nan(1,nsamples);
recall0 = nan(1,nsamples);
recall1 = nan(1,nsamples);

n0 = nnz(Y_test==0);
n1 = nnz(Y_test==1);
for i = 1:nsamples,
  thresh = dlogprob_samples(i);
  idx1 = dlogprob>=thresh;
  idx0 = dlogprob<=thresh;
  n0curr = nnz(idx0' & Y_test==0);
  n1curr = nnz(idx1' & Y_test==1);
  recall0(i) = n0curr/n0;
  recall1(i) = n1curr/n1;
  precision0(i) = n0curr/nnz(idx0);
  precision1(i) = n1curr/nnz(idx1);
end

handles = struct;
handles.fig = 4;
figure(handles.fig);
clf;
handles.pr0 = plot(recall0,precision0,'-','color',[0,.5,0],'linewidth',3);
hold on;
handles.pr1 = plot(recall1,precision1,'-','color',[0,0,0],'linewidth',3);
scatter(recall0,precision0,[],dlogprob_samples,'.');
scatter(recall1,precision1,[],dlogprob_samples,'.');
legend([handles.pr0,handles.pr1],{'none','turn'});
xlabel('Recall (|ypred==c & ytrue==c| / |ytrue==c|)');
ylabel('Precision (|ypred==c & ytrue==c| / |ypred==c|)');
colorbar;
colormap(jet(nsamples)*.7);
title('Precision-Recall, Cross-Validation');

%% interactive figure in which we can plot some of the window features

handles = struct;
handles.fig = 1;
axespos = [.05,.1,.7,.85];
border = .01;
border1 = .02;
sliderwidth = .0125;
textheight = .05;
sliderpos0 = [axespos(1)+axespos(3)+border1,...
  axespos(2),...
  sliderwidth,...
  axespos(4)];
textpos0 = [axespos(1)+axespos(3),...
  axespos(2)-textheight,...
  sliderwidth+2*border1,...
  textheight];
sliderpos1 = [sliderpos0(1)+sliderpos0(3)+border1,...
  axespos(2),...
  sliderwidth,...
  axespos(4)];
textpos1 = [sliderpos0(1)+sliderpos0(3),...
  axespos(2)+axespos(4),...
  sliderwidth+2*border1,...
  textheight];
listpos = [sliderpos1(1)+sliderpos1(3)+border1,...
  axespos(2),...
  1-(sliderpos1(1)+sliderpos1(3)+border+border1),...
  axespos(4)];
isfig = ishandle(handles.fig);
figure(handles.fig);
clf;
if ~isfig,
  pos = get(handles.fig,'Position');
  pos([3,4]) = [1600,600];
  set(handles.fig,'Position',pos);
end
slidermin = floor(min(dlogprob));
slidermax = ceil(max(dlogprob));
sliderstep = [.1,2]/(slidermax-slidermin);
handles.ax = axes('Parent',handles.fig,'Units','Normalized','Position',axespos);
handles.slider0 = uicontrol(handles.fig,'Style','slider',...
  'Units','normalized','Position',sliderpos0,...
  'Min',slidermin,'Max',slidermax,...
  'Value',maxdiff_logprob0,...
  'SliderStep',sliderstep);
handles.text0 = uicontrol(handles.fig,'Style','text',...
  'Units','normalized','Position',textpos0,...
  'String',sprintf('none<=%.1f',maxdiff_logprob0));
handles.text1 = uicontrol(handles.fig,'Style','text',...
  'Units','normalized','Position',textpos1,...
  'String',sprintf('%.1f<=turn',mindiff_logprob1));
handles.slider1 = uicontrol(handles.fig,'Style','slider',...
  'Units','normalized','Position',sliderpos1,...
  'Min',slidermin,'Max',slidermax,...
  'Value',mindiff_logprob1,...
  'SliderStep',sliderstep);
handles.listbox = uicontrol(handles.fig,'Style','listbox',...
  'Units','normalized','Position',listpos,...
  'String',statnames,...
  'Max',nfeatures,'Min',0,'Value',[]);

xlim = [0,size(X_known,2)];
ylim = [min(X_known(:)),max(X_known(:))];
handles.ylim = ylim;
set(handles.fig,'ToolBar','figure');
set(handles.ax,'XLim',xlim,'YLim',ylim);
hold(handles.ax,'on');
colors = lines(5)*.5+.5;
handles.patch = nan(1,5);

for v = 1:5,
  handles.patch(v) = patch(nan(1,4),nan(1,4),...
    colors(v,:),'EdgeColor','none');
end
% legend(handles.patch,'l=none, d=turn','l=none, d=unknown',...
%   'l=turn, d=none','l=turn, d=turn','l=turn,d=unknown');

handles.title = title('');
handles.data = [];
setappdata(handles.listbox,'handles',handles);
setappdata(handles.listbox,'x',X_known);
setappdata(handles.fig,'handles',handles);
setappdata(handles.fig,'ytrue',Y_test);
setappdata(handles.fig,'dlogprob',dlogprob);
DebugComputeWindowFeatures_SliderCallback(handles.slider0);
set(handles.listbox,'Callback',@DebugComputeWindowFeatures_ListBoxCallback);
set(handles.slider0,'Callback',@DebugComputeWindowFeatures_SliderCallback);
set(handles.slider1,'Callback',@DebugComputeWindowFeatures_SliderCallback);

%% interactive figure in which we can histogram some of the window features

nbins = 100;
edges = nan(nfeatures,nbins+1);
centers = nan(nfeatures,nbins);
frac = nan(nfeatures,nbins,4);
for i = 1:nfeatures,
  edges(i,:) = linspace(min(X_known(i,:)),max(X_known(i,:)),nbins+1);
  centers(i,:) = (edges(i,1:end-1)+edges(i,2:end))/2;
  edges(i,end) = inf;
end

handles = struct;
handles.fig = 2;
axespos = [.05,.1,.7,.85];
border = .01;
border1 = .02;
sliderwidth = .0125;
textheight = .05;
sliderpos0 = [axespos(1)+axespos(3)+border1,...
  axespos(2),...
  sliderwidth,...
  axespos(4)];
textpos0 = [axespos(1)+axespos(3),...
  axespos(2)-textheight,...
  sliderwidth+2*border1,...
  textheight];
sliderpos1 = [sliderpos0(1)+sliderpos0(3)+border1,...
  axespos(2),...
  sliderwidth,...
  axespos(4)];
textpos1 = [sliderpos0(1)+sliderpos0(3),...
  axespos(2)+axespos(4),...
  sliderwidth+2*border1,...
  textheight];
listpos = [sliderpos1(1)+sliderpos1(3)+border1,...
  axespos(2),...
  1-(sliderpos1(1)+sliderpos1(3)+border+border1),...
  axespos(4)];
isfig = ishandle(handles.fig);
figure(handles.fig);
clf;
if ~isfig,
  pos = get(handles.fig,'Position');
  pos([3,4]) = [1600,600];
  set(handles.fig,'Position',pos);
end
slidermin = floor(min(dlogprob));
slidermax = ceil(max(dlogprob));
sliderstep = [.1,2]/(slidermax-slidermin);
handles.ax = axes('Parent',handles.fig,'Units','Normalized','Position',axespos);
handles.slider0 = uicontrol(handles.fig,'Style','slider',...
  'Units','normalized','Position',sliderpos0,...
  'Min',slidermin,'Max',slidermax,...
  'Value',maxdiff_logprob0,...
  'SliderStep',sliderstep);
handles.text0 = uicontrol(handles.fig,'Style','text',...
  'Units','normalized','Position',textpos0,...
  'String',sprintf('none<=%.1f',maxdiff_logprob0));
handles.text1 = uicontrol(handles.fig,'Style','text',...
  'Units','normalized','Position',textpos1,...
  'String',sprintf('%.1f<=turn',mindiff_logprob1));
handles.slider1 = uicontrol(handles.fig,'Style','slider',...
  'Units','normalized','Position',sliderpos1,...
  'Min',slidermin,'Max',slidermax,...
  'Value',mindiff_logprob1,...
  'SliderStep',sliderstep);
handles.listbox = uicontrol(handles.fig,'Style','listbox',...
  'Units','normalized','Position',listpos,...
  'String',statnames,...
  'Max',1,'Min',0,'Value',1);
handles.data = [];

data = struct;
data.ytrue = Y_test;
data.dlogprob = dlogprob;
data.centers = centers;
data.edges = edges;
data.x = X_known;
data.colors = lines(6);
setappdata(handles.fig,'handles',handles);
setappdata(handles.fig,'data',data);
setappdata(handles.listbox,'handles',handles);
setappdata(handles.listbox,'centers',centers);
setappdata(handles.listbox,'legends',{'l=none','d=none','l=turn','d=turn'});
DebugComputeWindowFeatures_SliderCallback_Hist(handles.slider0);
DebugComputeWindowFeatures_ListBoxCallback_Hist(handles.listbox);
set(handles.listbox,'Callback',@DebugComputeWindowFeatures_ListBoxCallback_Hist);
set(handles.slider0,'Callback',@DebugComputeWindowFeatures_SliderCallback_Hist);
set(handles.slider1,'Callback',@DebugComputeWindowFeatures_SliderCallback_Hist);

%% look at effect of nferns

nsamples = 50;
min_nferns = 100;
max_nferns = 100000;
nferns_sample = unique(round(logspace(log10(min_nferns),log10(max_nferns),nsamples)));
fern_depth = 10;
threshold_range = [-4,4];

mindiff_logprob1s = [-10,0,10];
maxdiff_logprob0s = [-10,0,10];

% loop over all possible nferns
Y_test = Y_known;
n0 = nnz(Y_test==0);
n1 = nnz(Y_test==1);
precision0 = nan(numel(maxdiff_logprob0),nsamples);
precision1 = nan(numel(mindiff_logprob1),nsamples);
recall0 = nan(numel(maxdiff_logprob0),nsamples);
recall1 = nan(numel(mindiff_logprob1),nsamples);
for samplei = 1:nsamples,
  
  nferns = nferns_sample(samplei);
  fprintf('Nferns = %d...\n',nferns);
  params = struct('S',fern_depth,'M',nferns,'thrr',threshold_range);

  % all outputs
  dlogprob = nan(numel(Y_known),1);
  
  for expi_test = 1:nexps,

    fprintf('Experiment %d...\n',expi_test);
    
    idx_train = expi_known~=expi_test;
    idx_test = ~idx_train;

    Y_train = Y_known(idx_train);
    X_train = X_known(:,idx_train);
    X_test = X_known(:,idx_test);
    
    % z-score
    mu = nanmean(X_train,2);
    sig = nanstd(X_train,1,2);
    X_train_zscore = bsxfun(@rdivide,bsxfun(@minus,X_train,mu),sig);
    X_test_zscore = bsxfun(@rdivide,bsxfun(@minus,X_test,mu),sig);

    % train
    fern_model = fernsClfTrain( X_train_zscore',Y_train+1,params);

    [~,logprob] = fernsClfApply(X_test_zscore',fern_model);
    dlogprob(idx_test) = logprob(:,2) - logprob(:,1);

  end

  % compute precision, recalls
  for i = 1:numel(maxdiff_logprob0s),
    thresh = maxdiff_logprob0s(i);
    idx0 = dlogprob<=thresh;
    n0curr = nnz(idx0' & Y_test==0);
    recall0(i,samplei) = n0curr/n0;
    precision0(i,samplei) = n0curr/nnz(idx0);
  end

  for i = 1:numel(mindiff_logprob1s),
    thresh = maxdiff_logprob0s(i);
    idx1 = dlogprob>=thresh;
    n1curr = nnz(idx1' & Y_test==1);
    recall1(i,samplei) = n1curr/n1;
    precision1(i,samplei) = n1curr/nnz(idx1);
  end
  
end