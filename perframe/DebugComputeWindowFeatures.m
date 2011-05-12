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
labeldata = load(fullfile(rootdatadir,labelmatnames{1}));
paramsfilename = fullfile(rootdatadir,'SharpturnPerFrameParams.xml');


%% set parameters
[params,cellparams] = ReadPerFrameParams(paramsfilename);

%% do it

fns = fieldnames(cellparams);
y = cell(1,numel(fns));
feature_names = cell(1,numel(fns));
for i = 1:numel(fns),
  fn = fns{i};
  j = find(strcmp(labeldata.perframefns,fn));
  if isempty(j),
    error('Unknown perframe field %s',fn);
  end
  [y{i},feature_names{i}] = ...
    ComputeWindowFeatures(labeldata.perframedata{j},cellparams.(fn){:});  
end
%% plot stuff

figure(1);
clf;
plot(1:N,x,'ko-','markerfacecolor','k');
hold on;

fns_plot = {{'stat','mean','trans','none','radius',9,'offset',0},...
  {'stat','min','trans','flip','radius',11,'offset',0},...
  {'stat','max','trans','none','radius',3,'offset',-3},...
  {'stat','hist','trans','abs','radius',3,'offset',0,'bin',nbins},...
  {'stat','prctile','trans','abs','radius',20,'offset',0,'prctile',97.5},...
  {'stat','change','trans','none','radius',20,'offset',-20,'change_window_radius',2},...
  {'stat','std','trans','none','radius',20,'offset',0},...
  {'stat','harmonic','trans','none','radius',20,'offset',0,'num_harmonic',1},...
  {'stat','diff_neighbor_mean','trans','abs','radius',20,'offset',0},...
  {'stat','diff_neighbor_min','trans','abs','radius',20,'offset',0},...
  {'stat','diff_neighbor_max','trans','none','radius',20,'offset',0},...
  {'stat','zscore_neighbors','trans','none','radius',20,'offset',0}};

for i1 = 1:numel(fns_plot),
  fn1 = fns_plot{i1};
  disp(fn1);
  ismatch = false;
  % find the matching feature
  for i2 = 1:numel(feature_names),
    ismatch = true;
    fn2 = feature_names{i2};
    for j1 = 1:2:numel(fn1)-1,
      j2 = find(strcmp(fn1{j1},fn2(1:2:end)),1);
      if isempty(j2),
        ismatch = false;
        break;
      else
        j2 = j2*2-1;
        if ischar(fn1{j1+1}) && ~strcmpi(fn1{j1+1},fn2{j2+1}),
          ismatch = false;
          break;
        elseif ~ischar(fn1{j1+1}) && abs(fn1{j1+1}-fn2{j2+1}>.0001),
          ismatch = false;
          break;
        end
      end
    end
    if ismatch,
      break;
    end
  end
  if ~ismatch,
    fprintf('\n');
    warning('No match found for fns_plot %d',i1);
    disp(fns_plot{i1});
    fprintf('\n');
    input('Press enter to continue: ');    
    continue;
  end
  hplot = plot(1:N,y(i2,:),'r.-');
  disp(feature_names{i2});
  input('Press enter to continue: ');
  if ishandle(hplot) && i1 ~= numel(fns_plot),
    delete(hplot);
  end
end

