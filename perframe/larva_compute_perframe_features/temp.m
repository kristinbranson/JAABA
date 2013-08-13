trx = Trx('trxfilestr','trx.mat',...
                        'moviefilestr','movie.mov',...
                        'perframedir','perframe',...
                        'default_landmark_params',landmark_params,...
                        'perframe_params',perframe_params);
trx.AddExpDir(expdir);

% clean this data to force computation
if forcecompute,
  %deletefns = setdiff(perframefns,Trx.TrajectoryFieldNames());
  trx.CleanPerFrameData();
end

% compute each
for i = 1:nfns,
  fn = perframefns{i};
  fprintf(logfid,'Computing %s...\n',fn);
  trx.(fn); 
end