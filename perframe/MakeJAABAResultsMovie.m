function MakeJAABAResultsMovie(expdir,classifierparamsfile,varargin)

%% parse parameters

[framestarts,nframesperseg,fracstarts,fracperseg] = myparse(varargin,'framestarts',[],...
  'nframesperseg',[],'fracstarts',0,'fracperseg',1/60);

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

%% for video reading

[readframe,nframes,fid,headerinfo] = get_readframe_fcn(trx.movienames{1});

%% 

%% clean up

fclose(fid);