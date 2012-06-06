function MakeJAABAResultsMovie(expdir,classifierparamsfile,varargin)

classifierparams = ReadClassifierParamsFile(classifierparamsfile);
nbehaviors = numel(classifierparams);

%% create trx structure

trx = Trx('trxfilestr',classifierparams(1).file.trxfilename,...
  'moviefilestr',classifierparams(1).file.moviefilename,...
  'perframedir',classifierparams(1).file.perframedir);
if isfield(classifierparams,'perframe') && isfield(classifierparams(1).perframe,'params'),
  trx.SetPerFrameParams(classifierparams(1).perframe.params);
end  
if isfield(classifierparams,'perframe') && isfield(classifierparams(1).perframe,'landmark_params'),
  trx.SetLandmarkParams(classifierparams(1).perframe.landmark_params);
end  

trx.AddExpDir(expdir);

%% 

%% create JData structures

JData = cell(1,nbehaviors);
for i = 1:nbehaviors,
  JData{i} = JLabelData(classifierparams(i).configfile,...
    'setstatusfn',@fprintf_wrapper,'clearstatusfn',@() fprintf('Done.\n'));
  JData{i}.SetClassifierFileName(classifierparams(i).classifierfile,false);
  JData{i}.AddExpDirNoPreload(expdir);
end

%%  