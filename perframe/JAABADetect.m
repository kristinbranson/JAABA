function JAABADetect(expdir,classifierfiles,configfiles,varargin)

[blockSize] = myparse(varargin,'blockSize',5000);

nclassifiers = numel(classifierfiles);
if nclassifiers ~= numel(configfiles),
  error('Number of classifier files and number of config files do not match');
end

% read parameters
configparams = cell(1,nclassifiers);
classifiers = cell(1,nclassifiers);
for i = 1:nclassifiers,
  configparams{i} = ReadXMLParams(configfiles{i});
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

% get a list of all per-frame features
fprintf('Getting list of all per-frame features...\n');
pffs_per_classifier = cell(1,nclassifiers);
for i = 1:nclassifiers,
  fprintf('Classifier %d\n',i);
  pffs_per_classifier{i} = fieldnames(classifiers{i}.windowfeaturesparams);
end
pffs = unique([pffs_per_classifier{:}]);
npffs = numel(pffs);
pffidx_per_classifier = cell(1,nclassifiers);
for i = 1:nclassifiers,
  [~,pffidx_per_classifier{i}] = ismember(pffs_per_classifier{i},pffs);
end

% for each per-frame feature, get a list of window features
wfnames_per_classifier = cell(1,npffs);
wfs_per_classifier = cell(1,npffs);
wfnames = cell(1,npffs);
wfs = cell(1,npffs);
wfidx_per_classifier = cell(1,npffs);
for fni = 1:npffs,
  fn = pffs{fni};
  
  wfnames_per_classifier{fni} = cell(1,nclassifiers);
  wfs_per_classifier{fni} = cell(1,nclassifiers);
  
  for i = 1:nclassifiers,
    j = find(pffidx_per_classifier{i}==fni);
    if isempty(j),
      continue;
    end
    [~,feature_names_curr] = ComputeWindowFeatures([0,0],...
      classifiers{i}.windowfeaturescellparams.(fn){:});
    print_names_curr = WindowFeatureParams2String(fn,feature_names_curr);
    wfnames_per_classifier{fni}{i} = [wfnames_per_classifier{fni}{i},print_names_curr];
    wfs_per_classifier{fni}{i} = [wfs_per_classifier{fni}{i},feature_names_curr];
  end

  [wfnames{fni},idx] = unique([wfnames_per_classifier{fni}{:}]);
  tmp = [wfs_per_classifier{fni}{:}];
  wfs{fni} = tmp(idx);

  wfidx_per_classifier{fni} = cell(1,nclassifiers);
  for i = 1:nclassifiers,
    [~,wfidx_per_classifier{fni}{i}] = ismember(wfnames_per_classifier{fni}{i},wfnames);
  end
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
data = JLabelData(configfile_global);
[success,msg] = data.AddExpDir(expdir);
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

for flies = 1:data.nflies_per_exp(expi),
  
  % load per-frame data
  if flies == 1,
    perframedata = data.perframedata;
  else
    fprintf('Loading per-frame data for %s, flies %s',expdir,mat2str(flies));
    file = data.GetPerframeFiles(expi);
    for j = 1:numel(obj.allperframefns),
      if ~exist(file{j},'file'),
        error('Per-frame data file %s does not exist',file{j});
      end
    end
    tmp = load(file{j});
    perframedata{j} = tmp.data{flies(1)};
  end

  % frame for this fly
  tStart = data.firstframes_per_exp{expi}(flies);
  tEnd = data.endframes_per_exp{expi}(flies);
  scores{flies} = nan(nclassifiers,tEnd);
  
  % one block of frames at a time
  for t0 = tStart:blockSize:tEnd

    t1 = min(t0+blockSize-1,tEnd);
    
    % compute window features
    X = [];    
    % loop through per-frame features
    for j = 1:numel(pffs),
      fn = pffs{j};

      t11 = min(t1,numel(perframedata));
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
      X = [X,x_curr']; %#ok<AGROW>
      % add the feature names
      feature_names = [feature_names,cellfun(@(s) [{fn},s],feature_names_curr,'UniformOutput',false)]; %#ok<AGROW>
    end
    X = single(X);
    
    end
      
      
      
      X = JLabelData.ComputeWindowDataChunkStatic(pffs,pffs,
      perframefile,flies,windowfeaturescellparams,t0-tStart+1,curt1-tStart+1);
        
          scores(t0:t1) = myBoostClassify(X,classifier);
        end
        scoresA{flies} = scores;
        fprintf('Prediction done for %d fly, total number of flies:%d\n',flies,numFlies);
      end
  
  
  
  
  success = false; msg = '';
      
      if ~exist('mode','var'), mode = 'center'; end
      if ~exist('forceCalc','var'), forceCalc = false; end
      
      % Check if the features have been configured.
      if isempty(fieldnames(obj.windowfeaturesparams))
        obj.ShowSelectFeatures();
        if isempty(fieldnames(obj.windowfeaturesparams))
          error('No features selected!');
        end
      end
      
      % choose frames to compute:
      
      % bound at start and end frame of these flies
      T0 = max(obj.GetTrxFirstFrame(expi,flies));
      T1 = min(obj.GetTrxEndFrame(expi,flies));
      
      switch lower(mode),
        case 'center',
          % go forward r to find the end of the chunk
          t1 = min(t+obj.windowdatachunk_radius,T1);
          % go backward 2*r to find the start of the chunk
          t0 = max(t1-2*obj.windowdatachunk_radius,T0);
          % go forward 2*r again to find the end of the chunk
          t1 = min(t0+2*obj.windowdatachunk_radius,T1);
        case 'start',
          t0 = max(t,T0);
          t1 = min(t0+2*obj.windowdatachunk_radius,T1);
        case 'end',
          t1 = min(t,T1);
          t0 = max(t1-2*obj.windowdatachunk_radius,T0);
        otherwise
          error('Unknown mode %s',mode);
      end
      
      % find a continuous interval that covers all uncomputed ts between t0
      % and t1
      off = 1-t0;
      n = t1-t0+1;
      docompute = true(1,n);
      if ~isempty(obj.windowdata.exp) && ~forceCalc,
        tscomputed = obj.windowdata.t(obj.FlyNdx(expi,flies));
        tscomputed = tscomputed(tscomputed >= t0 & tscomputed <= t1);
        docompute(tscomputed+off) = false;
      end
      
      X = [];
      feature_names = {};
      if ~any(docompute),
        t1 = t0-1;
        success = true;
        return;
      end
      
      t0 = find(docompute,1,'first') - off;
      t1 = find(docompute,1,'last') - off;
      i0 = t0 - obj.GetTrxFirstFrame(expi,flies) + 1;
      i1 = t1 - obj.GetTrxFirstFrame(expi,flies) + 1;
      
      %       try
      
      curperframefns = obj.curperframefns;
      allperframefns = obj.allperframefns;
      perframeInMemory = ~isempty(obj.flies) && obj.IsCurFly(expi,flies);
      perframedata_all = obj.perframedata;
      perframefile = obj.GetPerframeFiles(expi);
      x_curr_all = cell(1,numel(curperframefns));
      feature_names_all = cell(1,numel(curperframefns));
      windowfeaturescellparams = obj.windowfeaturescellparams;
      
      % loop through per-frame fields to check that they exist.
      for j = 1:numel(curperframefns),
        fn = curperframefns{j};
        
        % get per-frame data
        ndx = find(strcmp(fn,allperframefns));
        if isempty(ndx),
          success = false;
          msg = 'Window features config file has a perframe feature that is not defined in params file';
          return;
        end
        
        if ~exist(perframefile{ndx},'file'),
          if ~isdeployed
            res = questdlg(sprintf('Experiment %s is missing some perframe files (%s, possibly more). Generate now?',obj.expnames{expi},perframefile{ndx}),'Generate missing files?','Yes','Cancel','Yes');
          else
            res = 'Yes';
          end
          
          if strcmpi(res,'Yes'),
            for ndx = 1:obj.nexps  
              [success1,msg1] = obj.GenerateMissingFiles(ndx);
              if ~success1,
                success = success1; msg = msg1;
                return;
              end
            end
            
          else
            success = false;
            msg = sprintf('Cannot compute window data for %s ',expdir);
            return;
          end

          
        end
        
      end
      
      parfor j = 1:numel(curperframefns),
        fn = curperframefns{j};
        
        % get per-frame data
        ndx = find(strcmp(fn,allperframefns));
        if perframeInMemory,
          perframedata = perframedata_all{ndx};
        else
          perframedata = load(perframefile{ndx});
          perframedata = perframedata.data{flies(1)};
        end
        
        i11 = min(i1,numel(perframedata));
        [x_curr,feature_names_curr] = ...
          ComputeWindowFeatures(perframedata,windowfeaturescellparams.(fn){:},'t0',i0,'t1',i11);
        if any(imag(x_curr(:)))
          fprintf('Feature values are complex, check input\n');
        end
        
        if i11 < i1,
          x_curr(:,end+1:end+i1-i11) = nan;
        end
        
        x_curr_all{j} = x_curr;
        feature_names_all{j} = feature_names_curr;
        
      end
      
      for j = 1:numel(curperframefns),
        fn = curperframefns{j};
        x_curr = x_curr_all{j};
        feature_names_curr = feature_names_all{j};
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
        X = [X,x_curr']; %#ok<AGROW>
        % add the feature names
        feature_names = [feature_names,cellfun(@(s) [{fn},s],feature_names_curr,'UniformOutput',false)]; %#ok<AGROW>
      end
      %       catch ME,
      %         msg = getReport(ME);
      %         return;
      %       end
      X = single(X);
      success = true;
  
  
end


data = JLabelData(configfilename);
data.SetClassifierFileNameWoExp(classifierfilename);

existingNdx = find(strcmp(expdir,data.expdirs));
if ~isempty(existingNdx)
  data.PredictWholeMovie(existingNdx);
  ndx = existingNdx;
else
  [success,msg] = data.AddExpDir(expdir);
  if ~success, fprintf(msg), return, end
  
  if ~data.filesfixable,
    fprintf('Experiment %s is missing required files that cannot be generated within this interface.',expdir); 
    return;
  end

  if data.filesfixable && ~data.allfilesexist,
    [success,msg] = data.GenerateMissingFiles(data.nexps,false);
    if ~success,
      fprintf('Error generating missing required files for experiment %s: %s. Removing...',expdir,msg);
      return;
    end
    
    [success,msg] = data.PreLoadLabeledData();
    if ~success,
      fprintf('Error computing window data for experiment %s: %s. Removing...',expdir,msg);
      return;
    end
  end
  
  ndx = find(strcmp(expdir,data.expdirs));
end

data.PredictWholeMovie(ndx);
sfn = data.GetFile('scores',ndx);
save(sfn,'classifierfilename','-append');