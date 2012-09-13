%% set up path

addpath ../misc;
addpath ../filehandling

%rootdatadir = '../experiments/wildtype/results';
if ispc,
  rootdatadir = 'C:\Data\JAABA\groundtruth_pBDPGAL4U_data';
else
  rootdatadir = '/groups/branson/bransonlab/projects/JAABA/data/groundtruth_pBDPGAL4U_data';
end

%% parameters

% experiment_names = {
%   'EXT_CSMH_None_Rig1Plate15BowlA_20120519T170213'
%   'FCF_attP2_1500062_None_Rig1Plate15BowlA_20120519T172815'
%   'FCF_cantons_1500002_None_Rig1Plate15BowlA_20120519T160453'
%   };
experiment_names = {
  'pBDPGAL4U_TrpA_Rig1Plate10BowlB_20110609T091905'
  'pBDPGAL4U_TrpA_Rig1Plate15BowlD_20120210T152130'
  'pBDPGAL4U_TrpA_Rig2Plate17BowlA_20120202T145723'
  'pBDPGAL4U_TrpA_Rig2Plate17BowlC_20110811T133805'
  'pBDPGAL4U_TrpA_Rig2Plate17BowlD_20111130T094443'
};
expdirs = cellfun(@(x) fullfile(rootdatadir,x),experiment_names,'UniformOutput',false);

classifierparamsfiles = {'../experiments/JAABA_classifier_params1.txt'
  '../experiments/JAABA_classifier_params2.txt'
  '../experiments/JAABA_classifier_params_wing.txt'};
nintervals = 10;
nframesperseg = 500;
keyword = 'GroundTruthBDP_2K';
timestamp = datestr(now,'yyyymmddTHHMMSS');
figpos = [100,100,1280,720];
forcecompute = false;

behaviorname_dict = {'pivot_tail','TailPivotTurn'
  'pivot_head','CenterPivotTurn'};

%% make sure all the scores are computed, up-to-date

classifierparams = [];
for filei = 1:numel(classifierparamsfiles),
  
  classifierparamscurr = ReadClassifierParamsFile(classifierparamsfiles{filei});
  if filei == 1,
    classifierparams = classifierparamscurr;
  else
    classifierparams = structappend(classifierparams,classifierparamscurr);
  end
end

nbehaviors = numel(classifierparams);
dobreak = false;
for i = 1:nbehaviors,
  if dobreak,
    break;
  end
  classifierfilecurr = classifierparams(i).classifierfile;
  configfilecurr = classifierparams(i).configfile;
  behavior = classifierparams(i).behaviors.names;
  if iscell(behavior),
    behavior = sprintf('%s_',behavior{:});
    behavior = behavior(1:end-1);
  end
  tmp = load(classifierfilecurr);
  classifier_timestamp = tmp.classifierTS;
  for j = 1:numel(expdirs),
    if dobreak,
      break;
    end
    scoresfile_curr = fullfile(expdirs{j},classifierparams(i).file.scorefilename);
    docompute = false;
    if ~exist(scoresfile_curr,'file'),
      fprintf('Scores file %s does not exist.\n',scoresfile_curr);
      docompute = true;
    else
      tmp = load(scoresfile_curr);
      scores_timestamp = tmp.timestamp;
      if scores_timestamp ~= classifier_timestamp,
        if forcecompute,
          docompute = true;
        else
          if scores_timestamp > classifier_timestamp,
            question = sprintf('Scores are newer than classifier for %s, %s. Recompute scores?',behavior,experiment_names{j});
            def = 'No';
          else
            question = sprintf('Scores are older than classifier for %s, %s. Recompute scores?',behavior,experiment_names{j});
            def = 'Yes';
          end
          res = questdlg(question,'Recompute?','Yes','No','Cancel',def);
          if strcmpi(res,'Yes'),
            docompute = true;
          elseif strcmpi(res,'Cancel'),
            dobreak = true;
            break;
          end
        end
      end
    end
    
    if docompute,
      fprintf('Computing %s...\n',scoresfile_curr);
      JLabelBatch(expdirs{j},'classifierfiles',classifierfilecurr,'configfiles',configfilecurr);
    end
    
  end
  
end



%% choose what to show

behaviorsuse = {'Chase','Jump','WingGrooming','Righting','Backup','PivotTail'};
weightingscheme = 'isbout';

[targets,expdirs_chosen,framestarts,nframesperseg] = ...
  ChooseExampleIntervals(expdirs,classifierparamsfile,...
  'nintervals',nintervals,'nframesperseg',nframesperseg,...
  'weightingscheme',weightingscheme);

% sort based on expdir, then framestart
maxframestart = max(framestarts);
[~,~,expi] = unique(expdirs_chosen);
tmp = expi'*10*maxframestart + framestarts;
[~,order] = sort(tmp);
targets = targets(order);
expdirs_chosen = expdirs_chosen(order);
framestarts = framestarts(order);
[unique_expdirs,~,expi] = unique(expdirs_chosen);

%% make the movies

avibasename = sprintf('JAABAResultsMovie_%s_%s',keyword,timestamp);

for i = 1:nintervals,
  [~,experiment_name] = fileparts(expdirs_chosen{i});
  aviname = sprintf('%s_%02d_%s_%02d.avi',avibasename,i,experiment_name,targets(i));
  MakeJAABAResultsMovie(expdirs_chosen{i},classifierparamsfiles,...
    'nframesperseg',nframesperseg,...
    'targets',targets(i),...
    'framestarts',framestarts(i),...
    'aviname',aviname,...
    'figpos',figpos);
end

%%

expi = 2;
expdir = expdirs{expi};
avibasename = sprintf('JAABAResultsMovie_%s_%s',keyword,timestamp);
targets = [9,18,3,17];
framestarts = [9646,10500,14201,14650];

[~,experiment_name] = fileparts(expdirs_chosen{i});
aviname = sprintf('%s_%02d_%s.avi',avibasename,i,experiment_name);
MakeJAABAResultsMovie(expdir,classifierparamsfiles,...
  'nframesperseg',nframesperseg,...
  'targets',targets,...
  'framestarts',framestarts,...
  'aviname',aviname,...
  'figpos',figpos,...
  'printnone',false);
