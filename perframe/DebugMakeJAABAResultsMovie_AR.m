%% set up path
%run from perframe dir

addpath ../misc;
addpath ../filehandling



%% parameters

rootdatadir = '/groups/branson/home/robiea/Projects_data/courtship';

experiment_names = {'FCF_pBDPGAL4U_1500437_TrpA_Rig2Plate17BowlD_20121121T152832'};

% list of jab files
% make sure experiment has current scores for all jab files listed -
% not checking for this anymore
classifierparamsfiles = {'/groups/branson/home/robiea/Projects_data/courtship/classifier_list_20140114.txt'};

behaviorname_dict = {'attempted_copulation','AttemptedCop'
   'chase','Chase'
   'copulation_cb','Copulation'
   'orientation_cb', 'Orientation'
   'wing_extension','WingExtension'};

%% set up manually


keyword = 'Courtship_152832';
timestamp = datestr(now,'yyyymmddTHHMMSS');

% for this set up target should be flies in the same roi, multiple
% segements per fly ok. length of targets, framestarts, nframesperseg should
% all be the same
targets = [5,5];
framestarts = [6483,28975];
nframesperseg = [368,219];
% size of movie
figpos = [100,100,800,720];


%% make the movies
expdirs = cellfun(@(x) fullfile(rootdatadir,x),experiment_names,'UniformOutput',false);
avibasename = sprintf('JAABAResultsMovie_%s_%s',keyword,timestamp);

% 1 movie at a time, 1 roi at a time
i = 1;
j = 1;
[~,experiment_name] = fileparts(expdirs{i});
% aviname = sprintf('%s_%02d_%s_%02d.avi',avibasename,experiment_name);
aviname = sprintf('%s_%s_%02d.avi',avibasename,experiment_name,targets(i));
trxfile = fullfile(expdirs{i},'registered_trx.mat');
td = load(trxfile);
roifile = fullfile(expdirs{i},'roidata.mat');
rd = load(roifile);

% not need if not using fixedaxis param
roi = td.trx(targets(j)).roi;
cx = rd.centerx(roi);
cy = rd.centery(roi);
r = rd.radii(roi);

fixedaxis = [cx-r-20,cx+r+20,cy-r-20,cy+r+20];

MakeJAABAResultsMovie(expdirs{i},classifierparamsfiles,...
   'nframesperseg',nframesperseg,...
   'targets',targets,...
   'framestarts',framestarts,...
   'aviname',aviname,...
   'figpos',figpos,...
   'behaviorname_dict',behaviorname_dict,...
   'fixedaxis',fixedaxis);




%% make sure all the scores are computed, up-to-date

% classifierparams = ReadClassifierParamsFile(classifierparamsfiles);
%
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
%         keyboard;
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
%       JAABADetect(expdirs{j},'classifierfiles',{classifierfilecurr},'configfiles',{configfilecurr});
%     end
%
%   end
%
% end
%


%% choose what to show
%
% %behaviorsuse = {'Chase','Jump','WingGrooming','Righting','Backup','pivot_tail','pivot_head','WingExtension','WingFlick','AttenptedCopulation'};
% behaviorsuse = {'WingExtension','WingFlick','AttemptedCopulation','Copulation','pivot_head'};
% weightingscheme = 'sumlogsumscore';
%
% [targets,expdirs_chosen,framestarts,nframesperseg] = ...
%   ChooseExampleIntervals(expdirs,classifierparamsfiles,...
%   'nintervals',nintervals,'nframesperseg',nframesperseg,...
%   'weightingscheme',weightingscheme);
%
% % sort based on expdir, then framestart
% maxframestart = max(framestarts);
% [~,~,expi] = unique(expdirs_chosen);
% tmp = expi'*10*maxframestart + framestarts;
% [~,order] = sort(tmp);
% targets = targets(order);
% expdirs_chosen = expdirs_chosen(order);
% framestarts = framestarts(order);
% [unique_expdirs,~,expi] = unique(expdirs_chosen);
%
% for i = 1:nintervals,
%   fprintf('%d: %s, fly %d, frame %d\n',i,expdirs_chosen{i},targets(i),framestarts(i));
% end






%% BDP

% keyword = 'GroundTruthBDP_2K';
% expi = 2;
% expdir = expdirs{expi};
% avibasename = sprintf('JAABAResultsMovie_%s_%s',keyword,timestamp);
% targets = [9,18,3,17];
% framestarts = [9646,10500,14201,14650];
% 
% [~,experiment_name] = fileparts(expdirs_chosen{expi});
% aviname = sprintf('%s_%02d_%s.avi',avibasename,expi,experiment_name);
% MakeJAABAResultsMovie(expdir,classifierparamsfiles,...
%    'nframesperseg',nframesperseg,...
%    'targets',targets,...
%    'framestarts',framestarts,...
%    'aviname',aviname,...
%    'figpos',figpos,...
%    'printnone',false,...
%    'behaviorname_dict',behaviorname_dict,...
%    'behaviororder',behaviororder,...
%    'textprefix','pBDPGAL4U',...
%    'textsuffix','1/3rd speed');

%%
% 
% keyword = 'GroundTruth71G01_1K';
% expi = 2;
% expi = 5;
% expdir = expdirs{expi};
% avibasename = sprintf('JAABAResultsMovie_%s_%s',keyword,timestamp);
% targets = [9,18,3,17];
% targets = [10,20];
% framestarts = [9646,10500,14201,14650];
% framestarts = [4396,7345];
% 
% [~,experiment_name] = fileparts(expdir);
% aviname = sprintf('%s_%02d_%s.avi',avibasename,expi,experiment_name);
% MakeJAABAResultsMovie(expdir,classifierparamsfiles,...
%    'nframesperseg',nframesperseg,...
%    'targets',targets,...
%    'framestarts',framestarts,...
%    'aviname',aviname,...
%    'figpos',figpos,...
%    'printnone',false,...
%    'behaviorname_dict',behaviorname_dict,...
%    'behaviororder',behaviororder,...
%    'textprefix','R71G01',...
%    'textsuffix','1/6th speed');