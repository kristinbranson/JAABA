% script for computing window features

% set up path
if ispc,
  JCtrax_path = 'E:\Code\JCtrax';
  FlyBowlAnalysis_path = 'E:\Code\FlyBowlAnalysis';
  rootdatadir = 'E:\Code\Jdetect\larva\fly_data\TrainingData';
  settingsdir = 'E:\Code\FlyBowlAnalysis\settings';
else
  JCtrax_path = '/groups/branson/home/bransonk/tracking/code/JCtrax';
  FlyBowlAnalysis_path = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis';
  settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
end
addpath(fullfile(JCtrax_path,'misc'));
addpath(fullfile(JCtrax_path,'filehandling'));
addpath(FlyBowlAnalysis_path);

%% load test data

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

paramsfilename = 'SharpturnPerFrameParams_AllParams.xml';

%% set parameters
[params,cellparams] = ReadPerFrameParams(paramsfilename);
fns = fieldnames(cellparams);

%% do it

X = cell(1,nexps);
Y = cell(1,nexps);
feature_names = {};

for expi = 1:nexps,

  labelfilename = fullfile(rootdatadir,labelmatnames{expi});
  fprintf('Computing window features for file %s...\n',labelmatnames{expi});
  if ~exist(labelfilename,'file'),
    warning('File %s does not exist, skipping',labelfilename);
    continue;
  end
  labeldata = load(labelfilename);
  for i = 1:numel(fns),
    fn = fns{i};
    fprintf('Per-frame feature %s...\n',fn);
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

%% combine stats to plot everything

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

%% interactive figure in which we can plot some of the window features

handles = struct;
handles.fig = 1;
axespos = [.05,.1,.7,.85];
border = .01;
listpos = [axespos(1)+axespos(3)+border,...
  axespos(2),...
  1-(axespos(1)+axespos(3)+2*border),...
  axespos(4)];
figure(handles.fig);
clf;
handles.ax = axes('Parent',handles.fig,'Units','Normalized','Position',axespos);
handles.listbox = uicontrol(handles.fig,'Style','listbox',...
  'Units','normalized','Position',listpos,...
  'String',statnames,...
  'Max',nfeatures,'Min',0,'Value',[]);
xlim = [0,N+1];
ylim = [min(X_combine(:)),max(X_combine(:))];
set(handles.fig,'ToolBar','figure');
set(handles.ax,'XLim',xlim,'YLim',ylim);
hold(handles.ax,'on');
[i0,i1] = get_interval_ends(Y_combine==1);
i1 = i1 - 1;
for i = 1:numel(i0),
  patch([i0(i)-.5,i1(i)+.5,i1(i)+.5,i0(i)-.5,i0(i)-.5],ylim([1,1,2,2,1]),...
    [1,.7,.7],'EdgeColor','none');
end
[i0,i1] = get_interval_ends(isnan(Y_combine));
i1 = i1 - 1;
for i = 1:numel(i0),
  patch([i0(i)-.5,i1(i)+.5,i1(i)+.5,i0(i)-.5,i0(i)-.5],ylim([1,1,2,2,1]),...
    [.5,.5,.5],'EdgeColor','none');
end
handles.title = title('');
handles.data = [];
setappdata(handles.listbox,'handles',handles);
setappdata(handles.listbox,'x',X_combine);
set(handles.listbox,'Callback',@DebugComputeWindowFeatures_ListBoxCallback);

%% interactive figure in which we can histogram some of the window features

nbins = 100;
edges = nan(nfeatures,nbins+1);
centers = nan(nfeatures,nbins);
frac = nan(nfeatures,nbins,2);
for i = 1:nfeatures,
  edges(i,:) = linspace(min(X_combine(i,:)),max(X_combine(i,:)),nbins+1);
  centers(i,:) = (edges(i,1:end-1)+edges(i,2:end))/2;
  edges(i,end) = inf;
  for j = 1:2,
    counts = histc(X_combine(i,Y_combine==j-1),edges(i,:));
    frac(i,:,j) = counts(1:end-1) / sum(counts(1:end-1));
  end
end

handles = struct;
handles.fig = 2;
axespos = [.05,.1,.7,.85];
border = .01;
listpos = [axespos(1)+axespos(3)+border,...
  axespos(2),...
  1-(axespos(1)+axespos(3)+2*border),...
  axespos(4)];
figure(handles.fig);
clf;
handles.ax = axes('Parent',handles.fig,'Units','Normalized','Position',axespos);
handles.listbox = uicontrol(handles.fig,'Style','listbox',...
  'Units','normalized','Position',listpos,...
  'String',statnames,...
  'Max',1,'Min',0,'Value',1);
set(handles.fig,'ToolBar','figure');
hold(handles.ax,'on');
handles.data = [];
setappdata(handles.listbox,'handles',handles);
setappdata(handles.listbox,'frac',frac);
setappdata(handles.listbox,'centers',centers);
setappdata(handles.listbox,'legends',{'none','sharpturn'});
DebugComputeWindowFeatures_ListBoxCallback_Hist(handles.listbox);
set(handles.listbox,'Callback',@DebugComputeWindowFeatures_ListBoxCallback_Hist);

