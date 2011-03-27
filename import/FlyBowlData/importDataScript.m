%% set up path

addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis;


%% creating sharpturn training data with 3 features 20110326:
% absdtheta
% dtheta
% velmag

rootdir = '../../../data/fly_data.sharpTurns.absdtheta_dtheta_velmag';
matfiledir = fullfile(rootdir,'TrainingData','MatFiles');
matfiles = dir(fullfile(matfiledir,'labeledsharpturns_*.mat'));
matfiles = cellfun(@(s)fullfile(matfiledir,s),{matfiles.name},'UniformOutput',false);
featureparamsfile = fullfile(rootdir,'Params','SharpturnSVMFeatureParams.txt');

% create experiment directories

expdirs = cell(1,numel(matfiles));
for i = 1:numel(matfiles),
  tmp = load(matfiles{i});
  inexpdir = fileparts(tmp.moviename);
  [~,basename] = fileparts(inexpdir);
  outexpdir = fullfile(matfiledir,basename);
%   if ~exist(outexpdir,'file'),
%     [success,msg] = copyfile(inexpdir,outexpdir);
%     if ~success,
%       error('Could not copy %s to %s:\n%s',inexpdir,outexpdir,msg);
%     end
%   end
  expdirs{i} = outexpdir;
end

% convert files
for i = 1:numel(matfiles),
  importData(matfiles{i},expdirs{i},featureparamsfile);
end
