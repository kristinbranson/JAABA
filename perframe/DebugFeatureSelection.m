%% feature selection script

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

%% parameters

windowfeaturesparamsfile_lots = 'WindowFeatures_Lots.xml';

nferns = 100;
fern_depth = 10;
threshold_range = [-4,4];
maxdiff_logprob0 = 0;
mindiff_logprob1 = 0;

fern_params = struct('S',fern_depth,'M',nferns,'thrr',threshold_range);
docomputewindowfeatures = false;
windowfeaturesmatfile = 'WindowFeatures20110608.mat';

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

%% choose histogram edges
% 
% fns = unique([perframefns{:}]);
% nbins = 10;
% edges = struct;
% LARGENUM = 1000000000;
% prctiles_edges = [0,linspace(1,99,nbins-1),100];
% for i = 1:numel(fns),
%   fn = fns{i};
%   data = [];
%   for expi = 1:nexps,
%     j = find(strcmp(perframefns{expi},fn),1);
%     data = [data,perframedata{expi}{j}];
%   end
% 
%   edges.(fn) = prctile(data,prctiles_edges);
%   edges.(fn)(1) = -LARGENUM;
%   edges.(fn)(end) = LARGENUM;
%     
% end
% 
% for i = 1:numel(fns),
%   fprintf('%s: ',fns{i});
%   fprintf('%f,',edges.(fns{i}));
%   fprintf('\n');
% end
% 
% % dcenter: -1000000000.000000,2.354018,3.505172,4.450626,5.725226,7.387971,9.806358,13.386464,19.999789,46.094717,1000000000.000000,
% % dell2nose: -1000000000.000000,0.908980,2.323822,3.511097,4.842711,6.509972,8.906403,12.438393,19.160503,45.235669,1000000000.000000,
% % dist2wall: -1000000000.000000,2.594928,5.118107,5.982994,6.774351,7.669668,8.790653,10.456486,14.438847,45.176156,1000000000.000000,
% % dnose2ell: -1000000000.000000,0.990906,2.418441,3.584737,4.850066,6.514153,8.883609,12.455677,19.033511,45.398485,1000000000.000000,
% % dtheta: -1000000000.000000,-8.215147,-1.983003,-0.804740,-0.281295,0.000394,0.282245,0.804283,1.974842,8.075577,1000000000.000000,
% % velmag: -1000000000.000000,0.000000,0.084922,0.485632,3.067131,7.485459,11.659870,15.393652,19.694228,33.017928,1000000000.000000,

%% read parameters

[params,cellparams] = ReadPerFrameParams(windowfeaturesparamsfile_lots);
fns = fieldnames(cellparams);

%% compute window features

if docomputewindowfeatures,

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
    fprintf('experimend %d, feature %s...\n',expi,fn);
    j = find(strcmp(labeldata.perframefns,fn));
    if isempty(j),
      error('Unknown perframe field %s for experiment %s',fn,labelmatnames{expi});
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

%% save window features

save(windowfeaturesmatfile,'X','Y','feature_names','windowfeaturesparamsfile_lots');

else
  
%% load pre-computed window features

load(windowfeaturesmatfile);

end

nfeatures = numel(feature_names);

%% choose features to add that most improve accuracy without cross-validation in a greedy manner

Xall = cat(2,X{:});
mu = nanmean(Xall,2);
sig = nanstd(Xall,1,2);
Xall_zscore = bsxfun(@rdivide,bsxfun(@minus,Xall,mu),sig);
Yall = cat(2,Y{:});
unknownidx = isnan(Yall);
Xall_zscore(:,unknownidx) = [];
Xall(:,unknownidx) = [];
Yall(unknownidx) = [];

matlabpool local 4;

features_selected = [];
mincosts = nan(1,nfeatures);
costs = nan(nfeatures,nfeatures);

fern_model_prev = [];
data = Xall_zscore';
hs = Yall'+1;

for iteri = iteri:nfeatures,
    
  data_prev = data(:,features_selected(1:iteri-1));
  if iteri > 1,
    fern_model_prev = fernsClfTrain(data_prev,hs,fern_params);
  else
    fern_model_prev = [];
  end

  costs_curr = nan(1,nfeatures);
  parfor featurei = 1:nfeatures,
    if ismember(featurei,features_selected),
      continue;
    end

    if ~isnan(costs_curr(featurei)), continue; end
    data_curr = [data_prev,data(:,featurei)];
    
    if iteri == 1,
      fern_model = fernsClfTrain(data_curr,hs,fern_params);
    else
      fern_model = fernsClfAddFeature(data_curr,hs,fern_model_prev,iteri,1:iteri-1);
    end
    
    [~,logprob] = fernsClfApply(data_curr,fern_model);
    dlogprob = logprob(:,2) - logprob(:,1);
    Ypred = 2+zeros(size(dlogprob));
    Ypred(dlogprob>=mindiff_logprob1) = 1;
    Ypred(dlogprob<=maxdiff_logprob0) = 0;
    costs_curr(featurei) = nnz(Ypred ~= 2 & Yall' ~= 2 & Ypred ~= Yall') + .5*nnz(Ypred == 2 & Ypred ~= Yall');
    fprintf('Cost = %d for feature %d: ',costs_curr(featurei),featurei); disp(feature_names{featurei});
  end
  [mincosts(iteri),features_selected(iteri)] = min(costs_curr); %#ok<SAGROW>
  costs(iteri,:) = costs_curr;
  fprintf('\n\n** Selected feature: '); disp(feature_names{features_selected(iteri)}); 
  if iteri == 1,
    fprintf('Cost = %d\n',mincosts(iteri));
  else
    fprintf('Decreased cost by %d - %d = %d to %d\n',...
      mincosts(iteri-1),mincosts(iteri),mincosts(iteri-1)-mincosts(iteri),mincosts(iteri));
  end
  save FeatureSelection_RandomFerns_20110612 iteri features_selected  mincosts fern_params mindiff_logprob1 maxdiff_logprob0;

end

%% test adding training data

N = size(data,1);
new_data_idx = N-5:N;
old_data_idx = setdiff(1:N,new_data_idx);

rng(0);
[ferns0,hsPr0] = fernsClfTrain(data,hs,fern_params);

rng(0);
[ferns1,hsPr1] = fernsClfTrain(data(old_data_idx,:),hs(old_data_idx),fern_params);
[ferns2,hsPr2] = fernsClfAddTrainingData(data,hs,new_data_idx,ferns1);

assert(all(ferns2.pFern(:)==ferns0.pFern(:)))
assert(all(ferns2.inds(:)==ferns0.inds(:)))
assert(all(hsPr2==hsPr0))

%% test removing training data

N = size(data,1);
remove_data_idx = N-5:N;
keep_data_idx = setdiff(1:N,remove_data_idx);

rng(0);
[ferns0,hsPr0] = fernsClfTrain(data(keep_data_idx,:),hs(keep_data_idx,:),fern_params);

rng(0);
[ferns1,hsPr1] = fernsClfTrain(data,hs,fern_params);
[ferns2,hsPr2] = fernsClfRemoveTrainingData(hs,remove_data_idx,ferns1);

assert(all(ferns2.pFern(:)==ferns0.pFern(:)))
assert(all(ferns2.inds(:)==ferns0.inds(:)))
assert(all(hsPr2==hsPr0))

%% test relabeling training data

N = size(data,1);
change_idx = 1:N-6;
hs_new = hs;
hs_new(change_idx) = 3 - hs(change_idx);

rng(0);
[ferns0,hsPr0] = fernsClfTrain(data,hs_new,fern_params);

rng(0);
[ferns1,hsPr1] = fernsClfTrain(data,hs,fern_params);
[ferns2,hsPr2] = fernsClfRelabelTrainingData(hs,hs_new,ferns1);

assert(all(ferns2.counts(:)==ferns0.counts(:)))
assert(all(ferns2.pFern(:)==ferns0.pFern(:)))
assert(all(ferns2.inds(:)==ferns0.inds(:)))
assert(all(hsPr2==hsPr0))

%% test changing number of ferns

Mnew = fern_params.M - 25;
tmp_fern_params = fern_params;
tmp_fern_params.M = Mnew;

rng(0);
[ferns0,hsPr0] = fernsClfTrain(data,hs,tmp_fern_params);

rng(0);
[ferns1,hsPr1] = fernsClfTrain(data,hs,fern_params);
[ferns2,hsPr2] = fernsClfChangeNFerns(data,hs,ferns1,Mnew,struct('thrr',fern_params.thrr));

% these won't turn out the same cuzza random number generation order
%assert(all(ferns2.counts(:)==ferns0.counts(:)))
%assert(all(ferns2.pFern(:)==ferns0.pFern(:)))
%assert(all(ferns2.inds(:)==ferns0.inds(:)))
%assert(all(hsPr2==hsPr0))

%% test cross-validation

Xall = cat(2,X{:});
mu = nanmean(Xall,2);
sig = nanstd(Xall,1,2);
Xall_zscore = bsxfun(@rdivide,bsxfun(@minus,Xall,mu),sig);
Yall = cat(2,Y{:});
cv_sets = nan(size(Yall));
off = 0;
for expi = 1:nexps,
  cv_sets(off+1:off+numel(Y{expi})) = expi;
  off = off + numel(Y{expi});
end
unknownidx = isnan(Yall);
Xall_zscore(:,unknownidx) = [];
Xall(:,unknownidx) = [];
Yall(unknownidx) = [];
cv_sets(unknownidx) = [];

fern_params_cv = fern_params;
fern_params_cv.cv_sets = cv_sets;
[fracwrong_cv,fracwrong,hsPr_cv,probs_cv,ferns,hsPr,probs] = fernsClfCrossValidation(Xall',Yall+1,fern_params_cv);

%% test feature replacement

fernsClfRemoveFeature(Xall',Yall+1,ferns,1)

%% try mrmr feature selection

% discretize the data
nbins_mrmr = 10;
prctiles = linspace(0,100,nbins_mrmr+1);
bin_edges = prctile(Xall',prctiles);
Xbin = nan(size(Xall));
for featurei = 1:nfeatures,
  [~,Xbin(featurei,:)] = histc(Xall(featurei,:),bin_edges(:,featurei)');
end
Xbin(Xbin < 1) = 1;
Xbin(Xbin > nbins_mrmr) = nbins_mrmr;
features_selected_mrmr = mrmr_mid_d(Xbin', Yall', nfeatures);
%
% The MID scheme of minimum redundancy maximal relevance (mRMR) feature selection
% 
% The parameters:
%  d - a N*M matrix, indicating N samples, each having M dimensions. Must be integers.
%  f - a N*1 matrix (vector), indicating the class/category of the N samples. Must be categorical.
%  K - the number of features need to be selected
%