function [issuccess,msgs] = CheckScores(expdirs,classifierparamsfiles,classifierparams)

if nargin < 3,
  classifierparams = ReadClassifierParamsFile(classifierparamsfiles);
end
nbehaviors = numel(classifierparams);
nexps = numel(expdirs);

issuccess = true(nbehaviors,nexps);
msgs = cell(1,nexps);

for i = 1:nbehaviors,

  classifierfilecurr = classifierparams(i).classifierfile;
  behavior = classifierparams(i).behaviors.names;
  if iscell(behavior),
    behavior = sprintf('%s_',behavior{:});
    behavior = behavior(1:end-1);
  end
  try
    tmp = load(classifierfilecurr,'classifierTS');
    classifier_timestamp = tmp.classifierTS;
  catch ME,
    warning('Could not load timestamp for %s classifier from %s\n%s',behavior,classifierfilecurr,getReport(ME));
    continue;
  end
  for j = 1:nexps,
    scoresfile_curr = fullfile(expdirs{j},classifierparams(i).file.scorefilename);
    if ~exist(scoresfile_curr,'file'),
      msgs{j}{end+1} = sprintf('Scores file for behavior %s %s does not exist.\n',behavior,scoresfile_curr);
      issuccess(i,j) = false;
      continue;
    end
    try
      tmp = load(scoresfile_curr,'timestamp');
      scores_timestamp = tmp.timestamp;
    catch ME,
      msgs{j}{end+1} = sprintf('Could not load timestamp for %s from %s\n%s',behavior,scoresfile_curr,getReport(ME));
      issuccess(i,j) = false;
      continue;
    end
    if scores_timestamp ~= classifier_timestamp,
      if scores_timestamp > classifier_timestamp,
        msgs{j}{end+1} = sprintf('Timestamp for %s scores = %s more recent than classifier = %s',...
          behavior,datestr(scores_timestamp),datestr(classifier_timestamp));
      else
        msgs{j}{end+1} = sprintf('Timestamp for %s scores = %s less recent than classifier = %s',...
          behavior,datestr(scores_timestamp),datestr(classifier_timestamp));
      end
      issuccess(i,j) = false;
      continue;
    end
    [~,name] = fileparts(expdirs{j});
    %fprintf('%s, %s OK\n',name,behavior);
  end
end