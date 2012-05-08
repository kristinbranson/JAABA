% script for testing the effects of random fern parameters

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
windowfeaturesparamsfiles = {windowfeaturesparamsfile_lots,windowfeaturesparamsfile_few};
windowfeaturesmatfiles = {windowfeaturesmatfile_lots,windowfeaturesmatfile_few};
nparamsfiles = numel(windowfeaturesparamsfiles);

% nferns to test
n_nferns_tests_lots = 50;
nferns_tests_lots = round(logspace(1,4,n_nferns_tests_lots));
n_nferns_tests_few = 5;
nferns_tests_few = round(logspace(log10(50),log10(5000),n_nferns_tests_few));

% fern depths to test
fern_depth_tests_lots = 1:15;
n_fern_depth_tests_lots = numel(fern_depth_tests_lots);
fern_depth_tests_few = [3,6,10];
n_fern_depth_tests_few = numel(fern_depth_tests_few);

% threshold range
threshold_range = [-4,4];

% thresholds for positives, negatives
maxdiff_logprob0 = 0;
mindiff_logprob1 = 0;

% whether to recompute window features
docomputewindowfeatures = false;

% fraction of cost to ignore for time speed up
epsilon = .05;

%% load per-frame data

perframedata = cell(1,nexps);
perframefns = cell(1,nexps);
Y = cell(1,nexps);
labelstarts = nan(1,nexps);
labelends = nan(1,nexps);
for expi = 1:nexps,

  labelfilename = fullfile(rootdatadir,labelmatnames{expi});
  if ~exist(labelfilename,'file'),
    warning('File %s does not exist, skipping',labelfilename);
    continue;
  end
  labeldata = load(labelfilename);
  perframedata{expi} = labeldata.perframedata;
  perframefns{expi} = labeldata.perframefns;
  labelstarts(expi) = labeldata.trk_labelstart;
  labelends(expi) = labeldata.trk_labelend;
  Y{expi} = labeldata.label;
end

%% read parameters

windowfeaturesparams = cell(1,nparamsfiles);
windowfeaturescellparams = cell(1,nparamsfiles);
windowfeaturefns = cell(1,nparamsfiles);
for i = 1:nparamsfiles,
  [windowfeaturesparams{i},windowfeaturescellparams{i}] = ReadPerFrameParams(windowfeaturesparamsfiles{i});
  windowfeaturefns{i} = fieldnames(windowfeaturescellparams{i});
end

%% compute window features

window_feature_names = cell(1,nparamsfiles);
windowfeatures = cell(1,nparamsfiles);
nfeatures = nan(1,nparamsfiles);
all_data = cell(1,nparamsfiles);
all_hs = cell(1,nparamsfiles);
all_expidx = cell(1,nparamsfiles);

for filei = 1:nparamsfiles,
  
  if docomputewindowfeatures || ~exist(windowfeaturesmatfiles{filei},'file'),

    X = cell(1,nexps);
    feature_names = {};
    cellparams = windowfeaturescellparams{filei};
    windowparamsfile = windowfeaturesparamsfiles{filei};
    fns = windowfeaturefns{filei};
    
    for expi = 1:nexps,
      
      for i = 1:numel(fns),
        fn = fns{i};
        fprintf('experimend %d, feature %s...\n',expi,fn);
        j = find(strcmp(perframefns{expi},fn));
        if isempty(j),
          error('Unknown perframe field %s for experiment %s',fn,labelmatnames{expi});
        end
        [x_curr,feature_names_curr] = ...
          ComputeWindowFeatures(perframedata{expi}{j},cellparams.(fn){:});
        if expi == 1,
          feature_names = [feature_names,cellfun(@(s) [{fn},s],feature_names_curr,'UniformOutput',false)]; %#ok<AGROW>
        end
        
        % crop out labeled part of video
        x_curr = x_curr(:,labelstarts(expi):labelends(expi)-1);
        
        % store
        X{expi} = [X{expi};x_curr];
      end
      
    end

    % save window features
    save(windowfeaturesmatfiles{filei},'X','Y','feature_names','windowparamsfile');
    
  else
  
    % load pre-computed window features
    load(windowfeaturesmatfiles{filei},'X','feature_names');

  end
  
  windowfeatures{filei} = X;
  window_feature_names{filei} = feature_names;
  nfeatures(filei) = numel(feature_names);
  
  % z-score
  Xall = cat(2,X{:});
  mu = nanmean(Xall,2);
  sig = nanstd(Xall,1,2);
  Xall_zscore = bsxfun(@rdivide,bsxfun(@minus,Xall,mu),sig);
  
  % which experiment
  expidx = nan(1,size(Xall,2));
  off = 0;
  for expi = 1:nexps,
    expidx(off+1:off+size(X{expi},2)) = expi;
    off = off + size(X{expi},2);
  end
  
  % remove unknowns
  Yall = cat(2,Y{:});
  unknownidx = isnan(Yall);
  Xall_zscore(:,unknownidx) = [];
  Yall(unknownidx) = [];
  expidx(unknownidx) = [];
  
  % store in format used in ferns
  all_data{filei} = Xall_zscore';
  all_hs{filei} = Yall'+1;
  all_expidx{filei} = expidx;
    
end

%% test the effect of number of ferns

confusion_matrix_training = nan([2,3,n_nferns_tests_lots,n_fern_depth_tests_few,nparamsfiles]);
confusion_matrix_cv = nan([2,3,n_nferns_tests_lots,n_fern_depth_tests_few,nparamsfiles]);
thrchosen = nan([n_nferns_tests_lots,n_fern_depth_tests_few,nparamsfiles]);

nlinesplot = nparamsfiles*n_fern_depth_tests_few;
colors = jet(nlinesplot)*.7;
hfig = 1;
figure(hfig);
clf;
hax = gca;
hold(hax,'on');
set(hax,'XScale','log');

hlines = nan(2,n_fern_depth_tests_few,nparamsfiles);
legends = cell([2,n_fern_depth_tests_few,nparamsfiles]);

for filei = 1:nparamsfiles,
  for fern_depth_i = 1:n_fern_depth_tests_few,
    hlines(1,fern_depth_i,filei) = plot(hax,nan,nan,'x-','color',colors(sub2ind([n_fern_depth_tests_few,nparamsfiles],fern_depth_i,filei),:));
    hlines(2,fern_depth_i,filei) = plot(hax,nan,nan,'o-','color',colors(sub2ind([n_fern_depth_tests_few,nparamsfiles],fern_depth_i,filei),:));
    legends{1,fern_depth_i,filei} = sprintf('File %d, depth=%d, neg',filei,fern_depth_tests_few(fern_depth_i));
    legends{2,fern_depth_i,filei} = sprintf('File %d, depth=%d, pos',filei,fern_depth_tests_few(fern_depth_i));
  end
end
legend(hax,hlines(:)',legends(:)');
xlabel(hax,'N. ferns');
ylabel('Frac wrong');

for filei = 1:nparamsfiles,

  data = all_data{filei};
  hs = all_hs{filei};
  Npos = nnz(hs==2);
  Nneg = nnz(hs==1);
  
  for fern_depth_i = 1:n_fern_depth_tests_few,
    
    fern_depth = fern_depth_tests_few(fern_depth_i);
    
    % train fern model with maximum number of ferns
    max_nferns = max(nferns_tests_lots);
    fern_params = struct('S',fern_depth,'M',max_nferns,'thrr',threshold_range);
    
    % reseed the random number generator for repeatable behavior
    rng(0);
    
    fern_model0 = fernsClfTrain(data,hs,fern_params);
    
    for nferns_i = 1:n_nferns_tests_lots,
      
      nferns = nferns_tests_lots(nferns_i);
      
      % train using all training data
      if nferns == max_nferns,
        fern_model = fern_model0;
      else
        
        % reseed the random number generator for repeatable behavior
        rng(0);
        
        fern_model = fernsClfChangeNFerns(data,hs,fern_model0,nferns,struct('thrr',threshold_range));
      end
      
      % compute training error
      [hsPr,probs] = fernsClfApply([],fern_model,fern_model.inds);
      % choose the threshold so that the fraction of false positives =
      % fraction of false negatives
      dprobs = diff(probs,1,2);
      [thrs,order] = sort(dprobs);
      hsorder = hs(order);
      fracwrong_pos = cumsum(double(hsorder==2))/Npos;
      fracwrong_neg = 1-cumsum(double(hsorder==1))/Nneg;
      tmp = find(fracwrong_pos >= fracwrong_neg,1);
      if isempty(tmp),
        thr_dlogprob = thrs(end);
      elseif tmp == 1,
        thr_dlogprob = thrs(1);
      else
        thr_dlogprob = (thrs(tmp)+thrs(tmp-1))/2;
      end
      thrchosen(nferns_i,fern_depth_i,filei) = thr_dlogprob;
      confusion_matrix_training(:,:,nferns_i,fern_depth_i,filei) = ...
        ComputeConfusionMatrix(hs,dprobs,thr_dlogprob,thr_dlogprob);
      
      % compute cross validation error
      fern_params_cv = fern_params;
      fern_params_cv.cv_sets = all_expidx{filei};
      fern_params_cv.M = nferns;
      fern_params_cv.ferns = fern_model;
      [~,~,~,probs_cv] = fernsClfCrossValidation( data, hs, fern_params_cv);
      confusion_matrix_cv(:,:,nferns_i,fern_depth_i,filei) = ...
        ComputeConfusionMatrix(hs,diff(probs_cv,1,2),thr_dlogprob,thr_dlogprob);

      fprintf('Nferns = %d, fern depth=%d, file %d, thr = %f, frac. negative wrong: %f(%f), frac. positive wrong: %f(%f).\n',...
        size(fern_model.fids,1),size(fern_model.fids,2),filei,...
        thr_dlogprob,...
        confusion_matrix_cv(1,2,nferns_i,fern_depth_i,filei),confusion_matrix_training(1,2,nferns_i,fern_depth_i,filei),...
        confusion_matrix_cv(2,1,nferns_i,fern_depth_i,filei),confusion_matrix_training(2,1,nferns_i,fern_depth_i,filei));
      
      xline = get(hlines(1,fern_depth_i,filei),'XData');
      yline_neg = get(hlines(1,fern_depth_i,filei),'YData');
      yline_pos = get(hlines(2,fern_depth_i,filei),'YData');
      xline(end+1) = nferns; %#ok<SAGROW>
      yline_neg(end+1) = confusion_matrix_cv(1,2,nferns_i,fern_depth_i,filei); %#ok<SAGROW>
      yline_pos(end+1) = confusion_matrix_cv(2,1,nferns_i,fern_depth_i,filei); %#ok<SAGROW>
      set(hlines(1,fern_depth_i,filei),'XData',xline,'YData',yline_neg);
      set(hlines(2,fern_depth_i,filei),'XData',xline,'YData',yline_pos);
      set(hax,'XLim',[min(nferns_tests_lots)-1,max(nferns_tests_lots)+1],'YLim',[0,1]);
      drawnow;
      
    end
    
    save EvaluateRandomFerns_NFerns20110617.mat confusion_matrix_cv confusion_matrix_training fern_depth_i filei thrchosen
    
  end
  
end

cost = reshape(confusion_matrix_cv(1,2,:,:,:)+confusion_matrix_cv(2,1,:,:,:)+...
  .5*(confusion_matrix_cv(1,3,:,:,:)+confusion_matrix_cv(2,3,:,:,:)),...
  [n_nferns_tests_lots,n_fern_depth_tests_few,nparamsfiles]);
cost_nfernstest = min(min(cost,[],2),[],3);

min_cost_nfernstest = min(cost_nfernstest(:));
best_nferns = nferns_tests_lots(find(cost_nfernstest <= min_cost_nfernstest*(1+epsilon),1));
% best_nferns = 54

hfig = 3;
figure(hfig);
clf;
plot(nferns_tests_lots,cost_nfernstest,'k.-');
set(gca,'xscale','log');
xlabel('Number of ferns');
ylabel('Cost');

%% test the effect of fern depth

confusion_matrix_training = nan([2,3,n_nferns_tests_few,n_fern_depth_tests_lots,nparamsfiles]);
confusion_matrix_cv = nan([2,3,n_nferns_tests_few,n_fern_depth_tests_lots,nparamsfiles]);
thrchosen = nan([n_nferns_tests_lots,n_fern_depth_tests_few,nparamsfiles]);

nlinesplot = nparamsfiles*n_nferns_tests_few;
colors = jet(nlinesplot)*.7;
hfig = 2;
figure(hfig);
clf;
hax = gca;
hold(hax,'on');
%set(hax,'XScale','log');

hlines = nan(2,n_nferns_tests_few,nparamsfiles);
legends = cell([2,n_nferns_tests_few,nparamsfiles]);

for filei = 1:nparamsfiles,
  for nferns_i = 1:n_nferns_tests_few,
    hlines(1,nferns_i,filei) = plot(hax,nan,nan,'x-','color',colors(sub2ind([n_nferns_tests_few,nparamsfiles],nferns_i,filei),:));
    hlines(2,nferns_i,filei) = plot(hax,nan,nan,'o-','color',colors(sub2ind([n_nferns_tests_few,nparamsfiles],nferns_i,filei),:));
    legends{1,nferns_i,filei} = sprintf('File %d, nferns=%d, neg',filei,nferns_tests_few(nferns_i));
    legends{2,nferns_i,filei} = sprintf('File %d, nferns=%d, pos',filei,nferns_tests_few(nferns_i));
  end
end
legend(hax,hlines(:)',legends(:)');
xlabel(hax,'Fern depth');
ylabel('Frac wrong');

for filei = 1:nparamsfiles,

  data = all_data{filei};
  hs = all_hs{filei};
  Npos = nnz(hs==2);
  Nneg = nnz(hs==1);
  
  for fern_depth_i = 16:n_fern_depth_tests_lots,
    
    fern_depth = fern_depth_tests_lots(fern_depth_i);
    
    % train fern model with maximum number of ferns
    max_nferns = max(nferns_tests_few);
    fern_params = struct('S',fern_depth,'M',max_nferns,'thrr',threshold_range);
    
    % reseed the random number generator for repeatable behavior
    rng(0);
    
    fern_model0 = fernsClfTrain(data,hs,fern_params);
    
    for nferns_i = 1:n_nferns_tests_few,
      
      nferns = nferns_tests_few(nferns_i);
      
      % train using all training data
      if nferns == max_nferns,
        fern_model = fern_model0;
      else
        
        % reseed the random number generator for repeatable behavior
        rng(0);
        
        fern_model = fernsClfChangeNFerns(data,hs,fern_model0,nferns,struct('thrr',threshold_range));
      end
      
      % compute training error
      [hsPr,probs] = fernsClfApply([],fern_model,fern_model.inds);
      
      
      % choose the threshold so that the fraction of false positives =
      % fraction of false negatives
      dprobs = diff(probs,1,2);
      [thrs,order] = sort(dprobs);
      hsorder = hs(order);
      fracwrong_pos = cumsum(double(hsorder==2))/Npos;
      fracwrong_neg = 1-cumsum(double(hsorder==1))/Nneg;
      tmp = find(fracwrong_pos >= fracwrong_neg,1);
      if isempty(tmp),
        thr_dlogprob = thrs(end);
      elseif tmp == 1,
        thr_dlogprob = thrs(1);
      else
        thr_dlogprob = (thrs(tmp)+thrs(tmp-1))/2;
      end
      thrchosen(nferns_i,fern_depth_i,filei) = thr_dlogprob;
      
      confusion_matrix_training(:,:,nferns_i,fern_depth_i,filei) = ...
        ComputeConfusionMatrix(hs,diff(probs,1,2),thr_dlogprob,thr_dlogprob);
      
      % compute cross validation error
      fern_params_cv = fern_params;
      fern_params_cv.cv_sets = all_expidx{filei};
      fern_params_cv.M = nferns;
      fern_params_cv.ferns = fern_model;
      [~,~,~,probs_cv] = fernsClfCrossValidation( data, hs, fern_params_cv);
      confusion_matrix_cv(:,:,nferns_i,fern_depth_i,filei) = ...
        ComputeConfusionMatrix(hs,diff(probs_cv,1,2),thr_dlogprob,thr_dlogprob);

      fprintf('Fern depth=%d, nferns = %d, file %d, thr = %f, frac. negative wrong: %f(%f), frac. positive wrong: %f(%f).\n',...
        size(fern_model.fids,2),size(fern_model.fids,1),filei,...
        thr_dlogprob,...
        confusion_matrix_cv(1,2,nferns_i,fern_depth_i,filei),confusion_matrix_training(1,2,nferns_i,fern_depth_i,filei),...
        confusion_matrix_cv(2,1,nferns_i,fern_depth_i,filei),confusion_matrix_training(2,1,nferns_i,fern_depth_i,filei));
      
      xline = get(hlines(1,nferns_i,filei),'XData');
      yline_neg = get(hlines(1,nferns_i,filei),'YData');
      yline_pos = get(hlines(2,nferns_i,filei),'YData');
      xline(end+1) = fern_depth; %#ok<SAGROW>
      yline_neg(end+1) = confusion_matrix_cv(1,2,nferns_i,fern_depth_i,filei); %#ok<SAGROW>
      yline_pos(end+1) = confusion_matrix_cv(2,1,nferns_i,fern_depth_i,filei); %#ok<SAGROW>
      set(hlines(1,nferns_i,filei),'XData',xline,'YData',yline_neg);
      set(hlines(2,nferns_i,filei),'XData',xline,'YData',yline_pos);
      set(hax,'XLim',[min(fern_depth_tests_lots)-1,max(fern_depth_tests_lots)+1],'YLim',[0,1]);
      drawnow;
      
    end
    
    save EvaluateRandomFerns_FernDepth20110617.mat confusion_matrix_cv confusion_matrix_training nferns_i filei thrchosen
    
  end
  
end

cost = reshape(confusion_matrix_cv(1,2,:,:,:)+confusion_matrix_cv(2,1,:,:,:)+...
  .5*(confusion_matrix_cv(1,3,:,:,:)+confusion_matrix_cv(2,3,:,:,:)),...
  [n_nferns_tests_few,n_fern_depth_tests_lots,nparamsfiles]);
cost_depthtest = min(min(cost,[],1),[],3);

hfig = 4;
figure(hfig);
clf;
plot(fern_depth_tests_lots,cost_depthtest,'k.-');
set(gca,'xscale','linear');
xlabel('Fern depth');
ylabel('Cost');

min_cost_depthtest = min(cost_depthtest(:));
best_fern_depth = fern_depth_tests_lots(find(cost_depthtest <= min_cost_depthtest*(1+epsilon),1));
% best_fern_depth = 10

%% feature set selection

fern_depth = 10;
nferns = 54;