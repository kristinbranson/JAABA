%% set up path

addpath ../misc;
addpath ../filehandling

%%

expdirs = {
  '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/GMR_26E12_AE_01_TrpA_Rig1Plate15BowlA_20111202T160314'
  '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/GMR_64D11_AE_01_TrpA_Rig1Plate10BowlB_20110401T144120'
  '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/GMR_21G01_AE_01_TrpA_Rig1Plate15BowlD_20110812T111328'
  '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/GMR_55C09_AE_01_TrpA_Rig2Plate14BowlC_20110330T151605'
  '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/GMR_69A12_AE_01_TrpA_Rig1Plate10BowlD_20110408T110922'
  '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/GMR_68B12_AE_01_TrpA_Rig1Plate15BowlD_20120208T161323'
  '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/GMR_92B03_AE_01_TrpA_Rig1Plate15BowlA_20111104T095956'
  '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/GMR_20C03_AE_01_TrpA_Rig1Plate10BowlA_20110720T133402'
  '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/GMR_16B01_AE_01_TrpA_Rig1Plate15BowlC_20120307T153820'
  '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/GMR_50B07_AE_01_TrpA_Rig2Plate14BowlA_20110331T092823'
  '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/specific_behaviors/GMR_67B01_AE_01_TrpA_Rig2Plate17BowlB_20111014T091451'
  '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/specific_behaviors/GMR_22E02_AE_01_TrpA_Rig2Plate17BowlD_20111201T093552'
  '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/specific_behaviors/GMR_13A08_AE_01_TrpA_Rig1Plate15BowlC_20111111T104418'
  '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/specific_behaviors/GMR_73E12_AE_01_TrpA_Rig1Plate15BowlD_20120517T140351'
  '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/specific_behaviors/GMR_52H10_AE_01_TrpA_Rig2Plate17BowlA_20110914T145037'
  '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/specific_behaviors/GMR_17F12_AE_01_TrpA_Rig1Plate15BowlD_20120411T092053'
  };



%% parameters

% experiment_names = {
%   'EXT_CSMH_None_Rig1Plate15BowlA_20120519T170213'
%   'FCF_attP2_1500062_None_Rig1Plate15BowlA_20120519T172815'
%   'FCF_cantons_1500002_None_Rig1Plate15BowlA_20120519T160453'
%   };

expdirs = {
  '/groups/branson/bransonlab/mayank/myFlyBowl/GMR_26E12_AE_01_TrpA_Rig1Plate15BowlA_20111202T160314'
  '/groups/branson/bransonlab/mayank/myFlyBowl/GMR_64D11_AE_01_TrpA_Rig1Plate10BowlB_20110401T144120'
  '/groups/branson/bransonlab/mayank/myFlyBowl/GMR_21G01_AE_01_TrpA_Rig1Plate15BowlD_20110812T111328'
  '/groups/branson/bransonlab/mayank/myFlyBowl/GMR_55C09_AE_01_TrpA_Rig2Plate14BowlC_20110330T151605'
  '/groups/branson/bransonlab/mayank/myFlyBowl/GMR_69A12_AE_01_TrpA_Rig1Plate10BowlD_20110408T110922'
  '/groups/branson/bransonlab/mayank/myFlyBowl/GMR_68B12_AE_01_TrpA_Rig1Plate15BowlD_20120208T161323'
  '/groups/branson/bransonlab/mayank/myFlyBowl/GMR_92B03_AE_01_TrpA_Rig1Plate15BowlA_20111104T095956'
  '/groups/branson/bransonlab/mayank/myFlyBowl/GMR_20C03_AE_01_TrpA_Rig1Plate10BowlA_20110720T133402'
  '/groups/branson/bransonlab/mayank/myFlyBowl/GMR_16B01_AE_01_TrpA_Rig1Plate15BowlC_20120307T153820'
  '/groups/branson/bransonlab/mayank/myFlyBowl/GMR_50B07_AE_01_TrpA_Rig2Plate14BowlA_20110331T092823'
  '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/specific_behaviors/GMR_67B01_AE_01_TrpA_Rig2Plate17BowlB_20111014T091451'
  '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/specific_behaviors/GMR_22E02_AE_01_TrpA_Rig2Plate17BowlD_20111201T093552'
  '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/specific_behaviors/GMR_13A08_AE_01_TrpA_Rig1Plate15BowlC_20111111T104418'
  '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/specific_behaviors/GMR_73E12_AE_01_TrpA_Rig1Plate15BowlD_20120517T140351'
  '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/specific_behaviors/GMR_52H10_AE_01_TrpA_Rig2Plate17BowlA_20110914T145037'
  '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/specific_behaviors/GMR_17F12_AE_01_TrpA_Rig1Plate15BowlD_20120411T092053'
  };

lineinfo = {
  'Mean velocity = 1.2 mm/s'
  'Mean velocity = 5.2 mm/s'
  'Mean velocity = 7.2 mm/s'
  'Mean velocity = 8.8 mm/s'
  'Mean velocity = 11.0 mm/s'
  'Mean distance to closest fly = 7.8 mm'
  'Mean distance to closest fly = 9.6 mm'
  'Mean distance to closest fly = 10.2 mm'
  'Mean distance to closest fly = 10.7 mm'
  'Mean distance to closest fly = 11.8 mm'
  'Fraction of time chasing = 0.17'
  'Fraction of time chasing = 0.08'
  'Fraction of time jumping = 0.01'
  'Fraction of time jumping = 0.07'
  'Fraction of time righting = 0.02'
  'Fraction of time righting = 0.09'
  };

experiment_names = cell(1,numel(expdirs));
for i = 1:numel(expdirs),
  [~,experiment_names{i}] = fileparts(expdirs{i});
end

classifierparamsfile = '../experiments/JAABA_classifier_params_gal4.txt';

nintervalsperexp = 2;
nframesperseg = 250;
keyword = 'GAL4Screen';
timestamp = datestr(now,'yyyymmddTHHMMSS');
figpos = [100,100,1280,720];
forcecompute = false;

%% make sure all the scores are computed, up-to-date
% 
% classifierparams = ReadClassifierParamsFile(classifierparamsfile);
% nbehaviors = numel(classifierparams);
% dobreak = false;
% for i = 1:nbehaviors,
%   if dobreak,
%     break;
%   end
%   classifierfilecurr = classifierparams(i).classifierfile;
%   configfilecurr = classifierparams(i).configfile;
%   behavior = classifierparams(i).behaviors.names;
%   if iscell(behavior),
%     behavior = sprintf('%s_',behavior{:});
%     behavior = behavior(1:end-1);
%   end
%   tmp = load(classifierfilecurr);
%   classifier_timestamp = tmp.classifierTS;
%   for j = 1:numel(expdirs),
%     if dobreak,
%       break;
%     end
%     scoresfile_curr = fullfile(expdirs{j},classifierparams(i).file.scorefilename);
%     docompute = false;
%     if ~exist(scoresfile_curr,'file'),
%       fprintf('Scores file %s does not exist.\n',scoresfile_curr);
%       docompute = true;
%     else
%       tmp = load(scoresfile_curr);
%       scores_timestamp = tmp.timestamp;
%       if scores_timestamp ~= classifier_timestamp,
%         if forcecompute,
%           docompute = true;
%         else
%           if scores_timestamp > classifier_timestamp,
%             question = sprintf('Scores are newer than classifier for %s, %s. Recompute scores?',behavior,experiment_names{j});
%             def = 'No';
%           else
%             question = sprintf('Scores are older than classifier for %s, %s. Recompute scores?',behavior,experiment_names{j});
%             def = 'Yes';
%           end
%           res = questdlg(question,'Recompute?','Yes','No','Cancel',def);
%           if strcmpi(res,'Yes'),
%             docompute = true;
%           elseif strcmpi(res,'Cancel'),
%             dobreak = true;
%             break;
%           end
%         end
%       end
%     end
%     
%     if docompute,
%       fprintf('Computing %s...\n',scoresfile_curr);
%       JLabelBatch(expdirs{j},classifierfilecurr,configfilecurr);
%     end
%     
%   end
%   
% end

% %% choose one male and female randomly
% 
% targets = nan(2,nexps);
% for i = 1:nexps,
%   
%   % load the trx
%   trxfile = fullfile(expdirs{i},'registered_trx.mat');
%   load(trxfile,'trx');
%   
%   % find one male and one female that is alive the entire interval
%   nflies = numel(trx);
%   isallowedmale = false(1,nflies);
%   isallowedfemale = false(1,nflies);
%   
%   t0 = T0 + min([trx.firstframe]) - 1;
%   t1 = t0+nframesplot-1;
%   
%   for fly = 1:nflies,
%     if trx(fly).firstframe > t0-buffert || trx(fly).endframe < t1 + buffert,
%       continue;
%     end
%     i0 = t0+trx(fly).off;
%     i1 = t1+trx(fly).off;
%     if all(strcmpi(trx(fly).sex(i0:i1),'M')),
%       isallowedmale(fly) = true;
%     elseif all(strcmpi(trx(fly).sex(i0:i1),'F')),
%       isallowedfemale(fly) = true;
%     end
%   end
%   flymale = randsample(find(isallowedmale),1);
%   flyfemale = randsample(find(isallowedfemale),1);
%   
%   % load the scores
%   targets(1,i) = flymale;
%   targets(2,i) = flyfemale;
%   
% end

%% choose what to show

behaviorsuse = {'Chase','Jump','Righting'};
weightingscheme = 'random';
downweight_samesex = 100000;
downweight_samefly = 100000;

nexps = numel(expdirs);

targets = nan(nintervalsperexp,nexps);
framestarts = nan(nintervalsperexp,nexps);
expdirs_chosen = cell(nintervalsperexp,nexps);

for i = 1:nexps,
  [targets(:,i),~,framestarts(:,i),~] = ...
    ChooseExampleIntervals(expdirs(i),classifierparamsfile,...
    'nintervals',nintervalsperexp,'nframesperseg',nframesperseg,...
    'weightingscheme',weightingscheme,...
    'downweight_samesex',downweight_samesex,...
    'downweight_samefly',downweight_samefly);
  [expdirs_chosen{:,i}] = deal(expdirs{i});
end

% sort based on expdir, then framestart
% maxframestart = max(framestarts(:));
% [~,~,expi] = unique(expdirs_chosen(:));
% tmp = expi'*10*maxframestart + framestarts;
% [~,order] = sort(tmp);
% targets = targets(order);
% expdirs_chosen = expdirs_chosen(order);
% framestarts = framestarts(order);
% [unique_expdirs,~,expi] = unique(expdirs_chosen);

%% make the movies

avibasename = sprintf('JAABAResultsMovie_%s_%s',keyword,timestamp);
pxwidthradius = 175;
pxheightradius = 175;
behavior2color = {'Chase',[0,.45,.7]
  'Jump',[.8 .4 0]
  'Righting',[0 .6 .5]};

for i = 1:nexps,
  
  [~,name] = fileparts(expdirs{i});
  shortname = regexprep(name,'^GMR_','');
  shortname = regexprep(shortname,'_A[ED]_01.*$','');
  avinamecurr = [avibasename,'_',shortname,'.avi'];
  
  [succeeded,aviname,figpos,height,width] = ...
    make_jaaba_result_movie_plotallflies(expdirs{i},...
    'aviname',avinamecurr,'nzoomr',5,'nzoomc',4,...
    'maxnframes',480,'firstframes',10000,...
    'figpos',[100,100,1280,720],'movietitle',lineinfo{i},...
    'doplotwings',false,...
    'classifierparamsfiles',{classifierparamsfile},...
    'debug',false,...
    'behavior2color',behavior2color);
  
end
  
%% compress

for i = 1:nexps,
  
  [~,name] = fileparts(expdirs{i});
  shortname = regexprep(name,'^GMR_','');
  shortname = regexprep(shortname,'_A[ED]_01.*$','');
  aviname = [avibasename,'_',shortname,'.avi'];
  
  
  [path,base,ext] = fileparts(aviname);
  tmpfile = fullfile(path,[base,'_tmp',ext]);
  newheight = 4*ceil(height/4);
  newwidth = 4*ceil(width/4);
  cmd = sprintf('mencoder %s -o %s -ovc xvid -xvidencopts fixed_quant=4 -vf scale=%d:%d',...
    aviname,tmpfile,newwidth,newheight);
  status = system(cmd);
  if status ~= 0,
    error('Failed to compress avi to %s',tmpfile);
  end
  unix(sprintf('mv %s %s',tmpfile,aviname));
  
end
 
%% concatenate

lineinfo = {
  'Mean velocity = 1.2 mm/s'
  'Mean velocity = 5.2 mm/s'
  'Mean velocity = 7.2 mm/s'
  'Mean velocity = 8.8 mm/s'
  'Mean velocity = 11.0 mm/s'
  'Mean distance to closest fly = 7.8 mm'
  'Mean distance to closest fly = 9.6 mm'
  'Mean distance to closest fly = 10.2 mm'
  'Mean distance to closest fly = 10.7 mm'
  'Mean distance to closest fly = 11.8 mm'
  'Fraction of time chasing = 0.17'
  'Fraction of time chasing = 0.08'
  'Fraction of time jumping = 0.01'
  'Fraction of time jumping = 0.07'
  'Fraction of time righting = 0.02'
  'Fraction of time righting = 0.09'
  };

lineorder = [11:16,1:10];
cmd = 'mencoder -oac copy -ovc copy -forceidx';
for i = 1:nexps,
  
  [~,name] = fileparts(expdirs{lineorder(i)});
  shortname = regexprep(name,'^GMR_','');
  shortname = regexprep(shortname,'_A[ED]_01.*$','');
  aviname = [avibasename,'_',shortname,'.avi'];
  cmd = [cmd,' ',aviname];
end
cmd = [cmd,sprintf(' -o %s.avi',avibasename)];

%%

%   
%   for jj = 1:nintervalsperexp,
%     i = sub2ind([nintervalsperexp,nexps],jj,ii);
%     
%     [~,experiment_name] = fileparts(expdirs_chosen{i});
%     aviname = sprintf('%s_%02d_%s_%02d.avi',avibasename,i,experiment_name,targets(i));
%     
%     load(fullfile(expdirs_chosen{i},'registered_trx.mat'),'trx');
%     fly = targets(i);
%     i0 = framestarts(i)+trx(fly).off;
%     i1 = i0 + nframesperseg - 1;
%     nmaleframes = nnz(strcmpi(trx(fly).sex(i0:i1),'M'));
%     if nmaleframes >= nframesperseg/2,
%       sex = 'male';
%     else
%       sex = 'female';
%     end
%     
%     titlepagetext = cell(1,2);
%     titlepagetext{1} = sprintf('Ground-truth test line %d, %s fly',ii,sex);
%     titlepagetext{2} = lineinfo{ii};
%     MakeJAABAResultsMovie(expdirs_chosen{i},classifierparamsfile,...
%       'nframesperseg',nframesperseg,...
%       'targets',targets(i),...
%       'framestarts',framestarts(i),...
%       'aviname',aviname,...
%       'figpos',figpos,...
%       'pxwidthradius',pxwidthradius,...
%       'pxheightradius',pxheightradius,...
%       'printnone',false,...
%       'titlepagetext',titlepagetext);
%   end
% end

%%

expi = 2;
expdir = expdirs{expi};
avibasename = sprintf('JAABAResultsMovie_%s_%s',keyword,timestamp);
targets = [9,18,3,17];
framestarts = [9646,10500,14201,14650];

[~,experiment_name] = fileparts(expdirs_chosen{i});
aviname = sprintf('%s_%02d_%s.avi',avibasename,i,experiment_name);
MakeJAABAResultsMovie(expdir,classifierparamsfile,...
  'nframesperseg',nframesperseg,...
  'targets',targets,...
  'framestarts',framestarts,...
  'aviname',aviname,...
  'figpos',figpos,...
  'printnone',false);
