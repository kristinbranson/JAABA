%% feature selection script

%% set up path

if ispc,
  JCtrax_path = 'E:\Code\JCtrax';
  FlyBowlAnalysis_path = 'E:\Code\FlyBowlAnalysis';
  rootdatadir = 'E:\Code\Jdetect\larva\fly_data\TrainingData';
  settingsdir = 'E:\Code\FlyBowlAnalysis\settings';
  mRMR_path = 'E:\Code\mRMR_0.9_compiled';
else
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

%% parameters

windowfeaturesparamsfile_lots = 'WindowFeatures_Lots.xml';


%  d - a N*M matrix, indicating N samples, each having M dimensions. Must be integers.
%  f - a N*1 matrix (vector), indicating the class/category of the N samples. Must be categorical.
%  K - the number of features need to be selected

%% data locations

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

%% load per-frame data

perframedata = cell(1,nexps);
perframefns = cell(1,nexps);
for expi = 1:nexps,

  labelfilename = fullfile(rootdatadir,labelmatnames{expi});
  if ~exist(labelfilename,'file'),
    warning('File %s does not exist, skipping',labelfilename);
    continue;
  end
  labeldata = load(labelfilename);
  perframedata{expi} = labeldata.perframedata;
  perframefns{expi} = labeldata.perframefns;
end

%% histogram

fns = unique([perframefns{:}]);
nbins = 10;
edges = struct;
LARGENUM = 1000000000;
for i = 1:numel(fns),
  fn = fns{i};
  data = [];
  for expi = 1:nexps,
    j = find(strcmp(perframefns{expi},fn),1);
    data = [data,perframedata{expi}{j}];
  end
  
  minv = min(data);
  maxv = max(data);
  if minv < 0,
    [edges.(fn),centers.(fn)] = SelectHistEdges(nbins,[minv,maxv],'logabs');
    edges.(fn)(1) = -LARGENUM;
    edges.(fn)(end) = LARGENUM;
  else
    [edges.(fn),centers.(fn)] = SelectHistEdges(nbins,[minv,maxv],'log');
    edges.(fn)(1) = -LARGENUM;
    edges.(fn)(end) = LARGENUM;    
  end
    
end

%%

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