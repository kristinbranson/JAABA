% scores = JAABADetect(expdir,'classifierfiles',classifierfiles,'configfiles',configfiles)
% scores = JAABADetect(expdir,'classifierparamsfile',classifierparamsfile)
function scores = JAABADetect(expdir,varargin)

[blockSize,classifierfiles,configfiles,classifierparamsfile,configparams] = ...
  myparse(varargin,'blockSize',10000,...
  'classifierfiles',{},'configfiles',{},...
  'classifierparamsfile',0,...
  'configparams',[]);

if ischar(classifierparamsfile),
  if ~exist(classifierparamsfile,'file'),
    error('File %s does not exist',classifierparamsfile);
  end
  fid = fopen(classifierparamsfile,'r');
  classifierfiles = {};
  configfiles = {};
  while true,
    l = fgetl(fid);
    if ~ischar(l),
      break;
    end
    if isempty(l),
      continue;
    end
    if strcmp(l(1),'%'), continue; end
    ws = regexp(l,',','split');
    classifierfiles{end+1} = ws{1}; %#ok<AGROW>
    configfiles{end+1} = ws{2}; %#ok<AGROW>
  end
  fclose(fid);
end

nclassifiers = numel(classifierfiles);
if nclassifiers == 0,
  error('No classifiers input');
end
% read parameters
fprintf('Reading classifier parameters...\n');
if isempty(configparams),
  if nclassifiers ~= numel(configfiles),
    error('Number of classifier files and number of config files do not match');
  end
  configparams = cell(1,nclassifiers);
  for i = 1:nclassifiers,
    configparams{i} = ReadXMLParams(configfiles{i});
  end
else
  if nclassifiers ~= numel(configparams),
    error('Number of classifier files and number of config params do not match');
  end
  if ~iscell(configparams),
    configparams = num2cell(configparams);
  end
end
classifiers = cell(1,nclassifiers);
for i = 1:nclassifiers,
  classifiers{i} = load(classifierfiles{i});
end

% check that the config files are compatible
configparams_global = configparams{1};
for i = 2:nclassifiers,
  if ~strcmp(configparams{i}.file.perframedir,configparams_global.file.perframedir),
    error('Incompatible perframedirs');
  end
  fns1 = fieldnames(configparams_global.perframe.params);
  fns2 = fieldnames(configparams{i}.perframe.params);
  fns = intersect(fns1,fns2);
  for j = 1:numel(fns),
    if ndims(configparams_global.perframe.params.(fns{j})) ~= ...
        ndims(configparams{i}.perframe.params.(fns{j})) || ...
        ~all(size(configparams_global.perframe.params.(fns{j})) == ...
        size(configparams{i}.perframe.params.(fns{j}))) || ...
        ~all(configparams_global.perframe.params.(fns{j})(:) == ...
        configparams{i}.perframe.params.(fns{j})(:)),
      error('Incompatible perframe params');
    end
  end
end

% get a list of all window features
wfidx_per_classifier = cell(1,nclassifiers);
wfs_per_classifier = cell(1,nclassifiers);
wfs = {};
fprintf('Getting list of window features...\n');
for i = 1:nclassifiers,
  
  fprintf('Classifier %d\n',i);
  
  % all per-frame features for this classifier
  pffs_curr = fieldnames(classifiers{i}.windowfeaturesparams);
  
  % list of all window features for this classifier
  feature_names = {};
  for j = 1:numel(pffs_curr),
    fn = pffs_curr{j};
    [~,feature_names_curr] = ComputeWindowFeatures([0,0],...
      classifiers{i}.windowfeaturescellparams.(fn){:});
    feature_names_curr = cellfun(@(x) [{fn},x],feature_names_curr,'UniformOutput',false);
    feature_names = [feature_names,feature_names_curr]; %#ok<AGROW>
  end
  
  % which features are actually used
  dims = [classifiers{i}.classifier.dim];
  feature_names = feature_names(dims);
  wfs_per_classifier{i} = feature_names;
  
  % put these in with the rest of the classifiers' window features
  wfidx_per_classifier{i} = nan(1,numel(dims));
  for j = 1:numel(feature_names),
    idxcurr = find(WindowFeatureNameCompare(feature_names{j},wfs),1);
    if isempty(idxcurr),
      idxcurr = numel(wfs)+1;
      wfs{idxcurr} = feature_names{j}; %#ok<AGROW>
    end
    wfidx_per_classifier{i}(j) = idxcurr;
  end
  
end
wf2pff = cellfun(@(x)x{1},wfs,'UniformOutput',false);
[pffs,~,wf2pffidx] = unique(wf2pff);
pfidx_per_classifier = cell(1,nclassifiers);
for i = 1:nclassifiers,
  pfidx_per_classifier{i} = wf2pffidx(wfidx_per_classifier{i});
end

% create parameters for ComputeWindowFeatures which will compute all these
% features
fprintf('Converting window feature names to ComputeWindowFeatures parameters...\n');
windowfeaturescellparams = struct;
for pfi = 1:numel(pffs),
  pf = pffs{pfi};
  wfidx = wf2pffidx==pfi;
  windowfeaturescellparams.(pf) = WindowFeatureName2Params(wfs(wfidx));
end

% make a temporary config file containing all features
configparams_global.behaviors.names = 'JAABADetect_DummyBehavior';
configparams_global.featureparamlist = struct;
for i = 1:numel(pffs),
  configparams_global.featureparamlist.(pffs{i}) = struct;
end
configfile_global = tempname();
SaveXMLParams(configparams_global,configfile_global,'params');

% initialize with temporary configfile
data = JLabelData(configfile_global,'openmovie',false);
[success,msg] = data.AddExpDirNoPreload(expdir);
if ~success,
  error(msg);
end
expi = 1;

% check if all the data inecessary is there
if ~data.filesfixable,
  fprintf('Experiment %s is missing required files that cannot be generated within this interface.',expdir);
  return;
end

% generate missing files if possible
if data.filesfixable && ~data.allfilesexist,
  [success,msg] = data.GenerateMissingFiles(data.nexps,false);
  if ~success,
    fprintf('Error generating missing required files for experiment %s: %s. Removing...',expdir,msg);
    return;
  end
end

isfirst = true;
scores = cell(1,nclassifiers);
for i = 1:nclassifiers,
  scores{i} = cell(1,data.nflies_per_exp(expi));
end

for flies = 1:data.nflies_per_exp(expi),
  
  fprintf('Classifying flies %s...\n',mat2str(flies));
  
  % load per-frame data
  if flies == 1,
    perframedata = data.perframedata;
  else
    file = data.GetPerframeFiles(expi);
    for j = 1:numel(file),
      if ~exist(file{j},'file'),
        error('Per-frame data file %s does not exist',file{j});
      end
      tmp = load(file{j});
      perframedata{j} = tmp.data{flies(1)};
    end
  end

  % frames for this fly
  tStart = 1;
  off = 1-data.firstframes_per_exp{expi}(flies);
  tEnd = data.endframes_per_exp{expi}(flies)+off;
  for i = 1:nclassifiers,
    scores{i}{flies} = nan(1,tEnd-off);
  end
  
  % one block of frames at a time
  for t0 = tStart:blockSize:tEnd,

    t1 = min(t0+blockSize-1,tEnd);
    
    % compute all window features
    X = single([]);
    feature_names = {};
    % loop through per-frame features
    for j = 1:numel(pffs),
      fn = pffs{j};

      t11 = min(t1,numel(perframedata{j}));
      [x_curr,feature_names_curr] = ...
        ComputeWindowFeatures(perframedata{j},windowfeaturescellparams.(fn){:},'t0',t0,'t1',t11);
      
      if t11 < t1,
        x_curr(:,end+1:end+t1-t11) = nan;
      end
      
      % add the window data for this per-frame feature to X
      nold = size(X,1);
      nnew = size(x_curr,2);
      if nold > nnew,
        warning('Number of examples for per-frame feature %s does not match number of examples for previous features',fn);
        x_curr(:,end+1:end+nold-nnew) = nan;
      elseif nnew > nold && ~isempty(X),
        warning('Number of examples for per-frame feature %s does not match number of examples for previous features',fn);
        X(end+1:end+nnew-nold,:) = nan;
      end
      X = [X,single(x_curr)']; %#ok<AGROW>
      % add the feature names
      feature_names = [feature_names,cellfun(@(s) [{fn},s],feature_names_curr,'UniformOutput',false)]; %#ok<AGROW>
    end
    
    % find indexes for each classifier into feature_names
    if isfirst,
      
      classifiers_indexed = cell(1,nclassifiers);

      for i = 1:nclassifiers,
        wfidx_per_classifier{i} = nan(1,numel(wfs_per_classifier{i}));
        for j = 1:numel(wfs_per_classifier{i}),
          idxcurr = find(WindowFeatureNameCompare(wfs_per_classifier{i}{j},feature_names));
          if numel(idxcurr) ~= 1,
            error('Error matching wfs for classifier %d with window features computed',i);
          end
          wfidx_per_classifier{i}(j) = idxcurr;
        end
        classifiers_indexed{i} = classifiers{i}.classifier;
        for j = 1:numel(classifiers_indexed{i}),
          classifiers_indexed{i}(j).dim = j;
        end
      end
      
      isfirst = false;
      
    end
    
    % apply the classifiers
    for i = 1:nclassifiers,
      Xcurr = X(:,wfidx_per_classifier{i});
      % there's some weird offseting in the stored scores files
      scores{i}{flies}(t0-off:t1-off) = myBoostClassify(Xcurr,classifiers_indexed{i});
    end

    fprintf('%d%% done.\n',round(100*(t1-tStart)/(tEnd-tStart)));
    
  end
  
end

% save scores
for i = 1:nclassifiers,
  allScores = struct;
  allScores.scores = scores{i};
  allScores.tStart = data.firstframes_per_exp{expi};
  allScores.tEnd = data.endframes_per_exp{expi};
  for flies = 1:data.nflies_per_exp(expi),
    [i0s i1s] = get_interval_ends(allScores.scores{flies}>0);
    allScores.t0s{flies} = i0s;
    allScores.t1s{flies} = i1s;
  end
  
  if ~isfield(configparams{i}.file,'scorefilename'),
    if ~isfield(configparams{i},'behaviors'),
      error('configparams %d does not have field behaviors',i);
    elseif ~isfield(configparams{i}.behaviors,'names'),
      error('configparams{%d}.behaviors does not have field names',i);
    else
      scorefilename = ['scores_',configparams{i}.behaviors.names,'.mat'];
    end
  else
    scorefilename = configparams{i}.file.scorefilename;
  end
  %scorefilename = ['test_',scorefilename];
  sfn = fullfile(expdir,scorefilename);
      
  if exist(sfn,'file'),
    [didbak,msg] = copyfile(sfn,[sfn,'~']);
    if ~didbak,
      warning('Could not create backup of %s: %s',sfn,msg);
    end
  end
  timestamp = classifiers{i}.classifierTS; %#ok<NASGU>
  classifierfilename = classifierfiles{i}; %#ok<NASGU>
  save(sfn,'allScores','timestamp','classifierfilename');  
end