function [classifiersets,featureconfigperset] = GroupClassifiers(configparams,scorefilenames)

nclassifiers = numel(configparams);
if nargin < 2,
  scorefilenames = cell(1,nclassifiers);
  for i = 1:nclassifiers,
    scorefilenames{i} = configparams{i}.file.scorefilename;
  end
end

% group into compatible config files
unique_featureconfigfiles = {};
repidx = [];
featureconfigidx = zeros(1,nclassifiers);

for i = 1:nclassifiers,
  didmatch = false;
  for j = 1:numel(unique_featureconfigfiles),
    if ~strcmp(unique_featureconfigfiles{j},configparams{i}.file.featureconfigfile),
      continue;
    end
    configparams_global = configparams{repidx(j)};
    if ~strcmp(configparams{i}.file.perframedir,configparams_global.file.perframedir),
      continue;
    end
    if ~(isfield(configparams_global,'perframe')&& isfield(configparams_global.perframe,'params')),
      didmatch = true;
      break;
    end
    fns1 = fieldnames(configparams_global.perframe.params);
    fns2 = fieldnames(configparams{i}.perframe.params);
    fns = intersect(fns1,fns2);
    docontinue = false;
    for k = 1:numel(fns),
      fn = fns{k};
      if ndims(configparams_global.perframe.params.(fn)) ~= ...
          ndims(configparams{i}.perframe.params.(fn)) || ...
          ~all(size(configparams_global.perframe.params.(fn)) == ...
          size(configparams{i}.perframe.params.(fn))) || ...
          ~all(configparams_global.perframe.params.(fn)(:) == ...
          configparams{i}.perframe.params.(fn)(:)),
        docontinue = true;
        break;
      end
    end
    if docontinue,
      continue;
    end
    didmatch = true;
    break;
  end
  if didmatch,
    featureconfigidx(i) = j;
  else
    featureconfigidx(i) = numel(unique_featureconfigfiles)+1;
    unique_featureconfigfiles{end+1} = configparams{i}.file.featureconfigfile; %#ok<AGROW>
    repidx(end+1) = i; %#ok<AGROW>
  end
  
end

featureconfigfiles = cell(1,nclassifiers);
for i = 1:nclassifiers,
  featureconfigfiles{i} = configparams{i}.file.featureconfigfile;
end
[unique_featureconfigfiles,~,featureconfigidx] = unique(featureconfigfiles);
featureconfigs = cell(1,numel(unique_featureconfigfiles));
for i = 1:numel(unique_featureconfigfiles),
  featureconfigs{i} = ReadXMLParams(unique_featureconfigfiles{i});
end


dependencies = cell(1,nclassifiers);
for i = 1:nclassifiers,
  if isfield(configparams{i},'scoresinput'),
    scorefilesinput = {configparams{i}.scoresinput.scorefilename};
    for k = 1:numel(scorefilesinput),
      if isempty(regexp(scorefilesinput{k},'\.mat$','once')),
        scorefilesinput{k} = [scorefilesinput{k},'.mat'];
      end
    end
    
    dependencies{i} = find(ismember(scorefilenames,scorefilesinput));
    
  else
    j = featureconfigidx(i);
    fnscurr = fieldnames(featureconfigs{j}.perframe);
    isscorefn = ~cellfun(@isempty,regexp(fnscurr,'^[sS]core','once'));
    if any(isscorefn),
      fnscurr = fnscurr(isscorefn);
      for k = 1:numel(fnscurr),
        if isempty(regexp(fnscurr{k},'\.mat$','once')),
          fnscurr{k} = [fnscurr{k},'.mat'];
        end
      end
      dependencies{i} = find(ismember(scorefilenames,fnscurr));
    end
  end
end

% remove groups that have no dependencies
tmp = cellfun(@isempty,dependencies);
idxnodependency = find(tmp);
idxdependency = find(~tmp);
classifiersets = {};
setdependencies = {};
featureconfigperset = {};
if ~isempty(idxnodependency),
  [fc,~,setidx] = unique(featureconfigfiles(idxnodependency));
  for i = 1:numel(fc),
    classifiersets{end+1} = idxnodependency(setidx==i); %#ok<AGROW>
    featureconfigperset{end+1} = fc{i}; %#ok<AGROW>
    setdependencies{end+1} = []; %#ok<AGROW>
  end
end

% this is not necessarily the most efficient way of ordering, but it should
% be safe
E = false(numel(idxdependency));
for i = 1:numel(idxdependency),
  [ism,idx] = ismember(dependencies{idxdependency(i)},idxdependency);
  if any(ism),
    E(idx,i) = true;
  end
end
order = graphtopoorder(sparse(E));
for i = idxdependency(order),
  fc = featureconfigfiles{i};
  if ~isempty(classifiersets) && strcmp(featureconfigperset{end},fc) && ...
      ~any(ismember(classifiersets{end},dependencies{i})) && ...
      isempty(setdiff(dependencies{i},setdependencies{end})) && ...
      isempty(setdiff(setdependencies{end},dependencies{i})),
    classifiersets{end}(end+1) = i;
  else
    classifiersets{end+1} = i; %#ok<AGROW>
    featureconfigperset{end+1} = fc; %#ok<AGROW>
    setdependencies{end+1} = dependencies{i}; %#ok<AGROW>
  end
end