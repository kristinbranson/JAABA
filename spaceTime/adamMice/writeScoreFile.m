function writeScoreFile(fname,allScores,behName,jab,clsTimestamp)
% Write a score file, with some checking for filesystem weirdness.
%
% jab: jabname or Macguffin

% XXXAL TODO 'classifierfilename' field should be renamed.
PtoSave = struct(...
  'behaviorName',behName,...
  'allScores',allScores,...
  'classifierfilename',jab,...
  'timestamp',clsTimestamp); %#ok<NASGU> Saved to mat file 

if exist(fname,'file'),
  [p,n,e] = fileparts(fname);
  bakfilename = fullfile(p,[n,'.bak',e]);
  if exist(bakfilename,'file'),
    try
      delete(bakfilename);
    end
  end
  try
    [success,msg] = movefile(fname,bakfilename);
    if ~success,
      error(msg);
    end
  catch ME,
    warning('Could not move existing scores file %s to backup location %s: %s',fname,bakfilename,getReport(ME));
  end
end
%fprintf('Saving scores %s at %s...\n',fname,datestr(now));
save(fname,'-struct','PtoSave');

for tryi = 1:5,
  try
    tmp = load(fname);
  catch ME,
    warning('Try %d: error testing whether we could reload the saved file %s: %s',tryi,fname,getReport(ME));
    continue;
  end
  break;
end