rootoutputdir = '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth';

tmp = dir(rootoutputdir);
experiment_names = {tmp.name};
idx = cellfun(@isempty,regexp(experiment_names,'_\d{8}T\d{6}','once'));
experiment_names(idx) = [];
expdirs = cellfun(@(s) fullfile(rootoutputdir,s),experiment_names,'UniformOutput',false);

for i = 1:numel(expdirs),
  expdir = expdirs{i};
  fprintf('%s\n',expdirs{i});
  perframedir = fullfile(expdir,'perframe');
  if ~exist(perframedir,'dir'),
    continue;
  end
  matfiles = dir(fullfile(perframedir,'*.mat'));
  for j = 1:numel(matfiles),
    filename = fullfile(perframedir,matfiles(j).name);
    [~,link] = unix(sprintf('readlink %s',filename));
    link = strtrim(link);
    if isempty(link),
      %fprintf('Keeping %s\n',filename);
      continue;
    end
    fprintf('Deleting %s\n',filename);
    delete(filename);
  end
end
