classdef JLabelData < handle
  
  properties (Access=public)
    
    % config parameters: locations of files
    moviefilename = 0;
    trxfilename = 0;
    labelfilename = 0;
    %windowfilename = 0;
    perframedir = 0;
    featureparamsfilename = 0;
    classifierfilename = '';
    rootoutputdir = 0;
    configfilename = '';
    
    % experiment directories
    expdirs = {};
    expnames = {};
    outexpdirs = {};
    nflies_per_exp = [];
    firstframes_per_exp = {};
    endframes_per_exp = {};
    nexps = 0;
    
    % files per experiment directory
    filetypes = {'movie','trx','label','perframedir'};%,'window'};
    
    % stuff stored in classifier mat file
    classifiervars = {'expdirs','outexpdirs','expnames','nflies_per_exp',...
      'firstframes_per_exp','endframes_per_exp',...
      'moviefilename','trxfilename','labelfilename','featureparamsfilename',...
      'configfilename','rootoutputdir','classifiertype','classifier','classifier_params'};%'windowfilename',
    
    % last used paths
    defaultpath = '';
    
    % window feature params
    windowfeaturesparams = struct;
    windowfeaturescellparams = {};
    perframefns = {};
    
    % whether files are missing
    fileexists = [];
    filetimestamps = [];
    allfilesexist = true;
    filesfixable = true;
    
    % current experiment
    expi = 0;
    % current flies
    flies = [];
    
    % last-used trajectories (one experiment, all flies)
    trx = {};

    % last-used per-frame data (one fly)
    perframedata = {};
    
    % computed and cached window features
    windowdata = struct('X',[],'exp',[],'flies',[],'t',[],...
      'labelidx_old',[],'labelidx_new',[],'featurenames',{{}},...
      'predicted',[],'predicted_probs',[],'isvalidprediction',[]);
    windowdatachunk_radius = 500;
    
    % labels structure
    labels = [];
    labelidx = [];
    labelidx_off = 0;
    t0_curr = 0;
    t1_curr = 0;
    predictedidx = [];
    erroridx = [];
    suggestedidx = [];
        
    % behavior names
    labelnames = {};
    nbehaviors = 0;
    
    % classifier
    classifiertype = 'ferns';
    classifier = [];
    classifier_params = struct;
    
    % functions for setting a status bar
    setstatusfn = '';
    clearstatusfn = '';
    
  end
  
  methods (Access=private)
    
    
  end
    
  methods (Access=public,Static=true)

    function res = IsRequiredFile(file)
      res = ismember(file,{'movie','trx','perframedir'});%,'window'});
    end    
    function res = CanGenerateFile(file)
      res = ismember(file,{'perframedir'});%,'window'});
    end
    function res = IsOutputFile(file)
      res = ismember(file,{'label'});%,'window'});
    end
    function res = IsPerFlyFile(file)
      res = ismember(file,{});%{'window'});
    end
    
  end
  
  methods (Access=public)
  
    function obj = JLabelData(configfilename,varargin)
      
      if mod(numel(varargin),2) ~= 0,
        error('Number of inputs to JLabelData constructor must be even');
      end
      
      % config file
      [success,msg] = obj.SetConfigFileName(configfilename);
      if ~success,
        error(msg);
      end
      
      % parse optional arguments in order
      s = varargin(1:2:end);
      v = varargin(2:2:end);
      
      % movie
      i = find(strcmpi(s,'moviefilename'),1);
      if ~isempty(i),
        [success,msg] = obj.SetMovieFileName(v{i});
        if ~success,
          error(msg);
        end

      end
      
      % trx
      i = find(strcmpi(s,'trxfilename'),1);
      if ~isempty(i),
        [success,msg] = obj.SetTrxFileName(v{i});
        if ~success,
          error(msg);
        end
      end
      
      % label
      i = find(strcmpi(s,'labelfilename'),1);
      if ~isempty(i),
        [success,msg] = obj.SetLabelFileName(v{i});
        if ~success,
          error(msg);
        end
      end
      
      % perframedir
      i = find(strcmpi(s,'perframedir'),1);
      if ~isempty(i),
        [success,msg] = obj.SetPerFrameDir(v{i});
        if ~success,
          error(msg);
        end
      end
      
      % featureparamsfilename
      i = find(strcmpi(s,'featureparamsfilename'),1);
      if ~isempty(i),
        [success,msg] = obj.SetFeatureParamsFileName(v{i});
        if ~success,
          error(msg);
        end
      end
      
      % rootoutputdir
      i = find(strcmpi(s,'rootoutputdir'),1);
      if ~isempty(i),
        [success,msg] = obj.SetRootOutputDir(v{i});
        if ~success,
          error(msg);
        end
      end
      
%       % window
%       i = find(strcmpi(s,'windowfilename'),1);
%       if ~isempty(i),
%         [success,msg] = obj.SetWindowFileName(v{i});
%         if ~success,
%           error(msg);
%         end
%       end
      
      % classifier
      i = find(strcmpi(s,'classifierfilename'),1);
      if ~isempty(i),
        [success,msg] = obj.SetClassifierFileName(v{i});
        if ~success,
          error(msg);
        end
      end
      
      % default path
      i = find(strcmpi(s,'defaultpath'),1);
      if ~isempty(i),
        [success,msg] = obj.SetDefaultPath(v{i});
        if ~success,
          error(msg);
        end
      end
      
      i = find(strcmpi(s,'setstatusfn'),1);
      if ~isempty(i),
        obj.setstatusfn = v{i};
      end

      i = find(strcmpi(s,'clearstatusfn'),1);
      if ~isempty(i),
        obj.clearstatusfn = v{i};
      end
      
      % make sure everything gets set, one way or another
      requiredfns = {'moviefilename','trxfilename','labelfilename','rootoutputdir'}; % 'windowfilename',
      for i = 1:numel(requiredfns),
        fn = requiredfns{i};
        if isnumeric(obj.(fn)),
          error('%s did not get initialized',fn);
        end
      end
      
      [success,msg] = obj.UpdateStatusTable();
      if ~success,
        error(msg);
      end
      
    end

    function [success,msg] = SetConfigFileName(obj,configfilename)
      
      success = false;
      msg = '';
      if ~ischar(configfilename),
        return;
      end
      try
        configparams = ReadXMLParams(configfilename);
      catch ME,
        msg = sprintf('Error reading config file %s: %s',configfilename,getReport(ME));
        return;
      end
      obj.configfilename = configfilename;
      if isfield(configparams,'file'),
        if isfield(configparams.file,'moviefilename'),
          [success1,msg] = obj.SetMovieFileName(configparams.file.moviefilename);
          if ~success1,
            return;
          end
        end
        if isfield(configparams.file,'trxfilename'),
          [success1,msg] = obj.SetTrxFileName(configparams.file.trxfilename);
          if ~success1,
            return;
          end
        end
        if isfield(configparams.file,'labelfilename'),
          [success1,msg] = obj.SetLabelFileName(configparams.file.labelfilename);
          if ~success1,
            return;
          end
        end
        if isfield(configparams.file,'perframedir'),
          [success1,msg] = obj.SetPerFrameDir(configparams.file.perframedir);
          if ~success1,
            return;
          end
        end
%         if isfield(configparams.file,'windowfilename'),
%           [success1,msg] = obj.SetWindowFileName(configparams.file.windowfilename);
%           if ~success1,
%             return;
%           end
%         end
        if isfield(configparams.file,'rootoutputdir'),
          [success1,msg] = obj.SetRootOutputDir(configparams.file.rootoutputdir);
          if ~success1,
            return;
          end
        end
        if isfield(configparams.file,'featureparamfilename'),
          [success1,msg] = obj.SetFeatureParamsFileName(configparams.file.featureparamfilename);
          if ~success1,
            return;
          end
        end
      end
      
      if isfield(configparams,'behaviors'),
        
        % read in behavior names
        if isfield(configparams.behaviors,'names'),
          obj.labelnames = configparams.behaviors.names;
          if ~iscell(obj.labelnames),
            obj.labelnames = {obj.labelnames};
          end
          % add none label
          if ~ismember('none',lower(obj.labelnames)),
            obj.labelnames{end+1} = 'None';
          end
        else
          obj.labelnames = {'Behavior','None'};
        end
          
        obj.nbehaviors = numel(obj.labelnames);
        
        % rearrange so that None is the last label
        nonei = find(strcmpi('None',obj.labelnames),1);
        obj.labelnames = obj.labelnames([1:nonei-1,nonei+1:obj.nbehaviors,nonei]);
        
      end
      
    end

    function [success,msg] = PreLoadWindowData(obj,expi,flies,ts)
      
      success = false;
      msg = '';
      
      if numel(expi) ~= 1,
        error('Usage: expi must be a scalar');
      end
      
      if size(flies,1) ~= 1,
        error('Usage: one set of flies must be selected');
      end
      
      if isempty(obj.windowdata.exp),
        missingts = ts;
        tscurr = [];
      else      
        idxcurr = obj.windowdata.exp == expi & all(bsxfun(@eq,obj.windowdata.flies,flies),2);
        tscurr = obj.windowdata.t(idxcurr);
        missingts = setdiff(ts,tscurr);
      end
        
      if isempty(missingts),
        success = true;
        return;
      end

      % get labels for current flies
      [labelidx,t0_labelidx] = obj.GetLabelIdx(obj.expi,flies);
      
      nts0 = numel(missingts);
      
      while true,
        
        t = missingts(1);
        obj.SetStatus('Computing window data %d%% done...',round(100*(nts0-numel(missingts))/nts0));

        [success1,msg,t0,t1,X,feature_names] = obj.ComputeWindowDataChunk(expi,flies,t);
        if ~success1,
          return;
        end
        
        tsnew = t0:t1;
        idxnew = ~ismember(tsnew,tscurr);
        m = nnz(idxnew);
        
        obj.windowdata.X(end+1:end+m,:) = X(idxnew,:);
        obj.windowdata.exp(end+1:end+m,1) = expi;
        obj.windowdata.flies(end+1:end+m,:) = repmat(flies,[m,1]);
        obj.windowdata.t(end+1:end+m,1) = tsnew(idxnew);
        obj.windowdata.labelidx_old(end+1:end+m,1) = 0;
        obj.windowdata.labelidx_new(end+1:end+m,1) = labelidx(t0-t0_labelidx+1:t1-t0_labelidx+1);
        obj.windowdata.predicted(end+1:end+m,1) = 0;
        obj.windowdata.isvalidprediction(end+1:end+m,1) = false;

        missingts(missingts >= t0 & missingts <= t1) = [];
        
        if isempty(missingts),
          obj.ClearStatus();
          break;
        end
        
      end
      obj.windowdata.featurenames = feature_names;
      
      success = true;
      
    end

    function [success,msg,t0,t1,X,feature_names] = ComputeWindowDataChunk(obj,expi,flies,t)
      
      success = false;
      msg = '';
      
      % choose frames to compute
      T0 = max(obj.GetTrxFirstFrame(expi,flies));
      T1 = min(obj.GetTrxEndFrame(expi,flies));
      t1 = min(t+obj.windowdatachunk_radius,T1);
      t0 = max(t1-2*obj.windowdatachunk_radius,T0);
      t1 = min(t0+2*obj.windowdatachunk_radius,T1);
      off = 1-t0;
      n = t1-t0+1;
      docompute = true(1,n);
      if ~isempty(obj.windowdata.exp),
        tscomputed = obj.windowdata.t(obj.windowdata.exp == expi & all(bsxfun(@eq,obj.windowdata.flies,flies),2));
        tscomputed = tscomputed(tscomputed >= t0 & tscomputed <= t1);
        docompute(tscomputed+off) = false;
      end
      t0 = find(docompute,1,'first') - off;
      t1 = find(docompute,1,'last') - off;

      try
        X = [];
        feature_names = {};

        % loop through per-frame fields
        for j = 1:numel(obj.perframefns),
          fn = obj.perframefns{j};

          % get per-frame data
          if all(flies == obj.flies),
        
            % use pre-loaded per-frame data
            [x_curr,feature_names_curr] = ...
              ComputeWindowFeatures(obj.perframedata{j},obj.windowfeaturescellparams.(fn){:},'t0',t0,'t1',t1);
            
          else
            
            % load in data
            perframedir = obj.GetFile('perframedir',expi);
            perframedata = load(fullfile(perframedir,[fn,'.mat']));
            % TODO: adapt for multiple flies labeled at once
            [x_curr,feature_names_curr] = ...
              ComputeWindowFeatures(perframedata.data{fly},obj.windowfeaturescellparams.(fn){:},'t0',t0,'t1',t1);
            
          end
          
          nold = size(X,1);
          nnew = size(x_curr,2);
          if nold > nnew,
            x_curr(:,end+1:end+nold-nnew) = nan;
          elseif nnew > nold,
            X(end+1:end+nnew-nold,:) = nan;
          end
          X = [X,x_curr']; %#ok<AGROW>
          feature_names = [feature_names,cellfun(@(s) [{fn},s],feature_names_curr,'UniformOutput',false)]; %#ok<AGROW>
        end
      catch ME,
        msg = getReport(ME);
        return;
      end
      
      success = true;
     
    end
    
    function [success,msg] = SetMovieFileName(obj,moviefilename)

      success = false;
      msg = '';

      if ischar(moviefilename),
        if ischar(obj.moviefilename) && strcmp(moviefilename,obj.moviefilename),
          success = true;
          return;
        end
        oldmoviefilename = obj.moviefilename;
        obj.moviefilename = moviefilename;
        [success1,msg] = obj.CheckMovies();
        if ~success1,
          obj.moviefilename = oldmoviefilename;
          return;
        end
        [success,msg] = obj.UpdateStatusTable('movie');
        % TODO: remove bad experiments
      end
      
    end
    
    function [successes,msg] = CheckMovies(obj,expis)
      
      successes = [];
      msg = '';
      
      if nargin < 2,
        expis = 1:obj.nexps;
      end
      
      if isempty(expis),
        return;
      end
      
      successes = true(1,numel(expis));
      for i = 1:numel(expis),
        moviefilename = obj.GetFile('movie',expis(i));
        obj.SetStatus('Checking movie %s...',moviefilename);
        if ~exist(moviefilename,'file'),
          successes(i) = false;
          msg1 = sprintf('File %s missing',moviefilename);
          if isempty(msg),
            msg = msg1;
          else
            msg = sprintf('%s\n%s',msg,msg1);
          end
        else
          
          try
            [readframe,~,movie_fid] = ...
              get_readframe_fcn(moviefilename);
            if movie_fid <= 0,
              error('Could not open movie %s for reading',moviefilename);
            end
            readframe(1);
            fclose(movie_fid);
          catch ME,
            successes(i) = false;
            msg1 = sprintf('Could not parse movie %s: %s',moviefilename,getReport(ME));
            if isempty(msg),
              msg = msg1;
            else
              msg = sprintf('%s\n%s',msg,msg1);
            end
          end
          
        end
      end
      
      obj.ClearStatus();
      
    end
    
    function [success,msg] = SetTrxFileName(obj,trxfilename)
      
      success = false;
      msg = '';
      if ischar(trxfilename),
        if ischar(obj.trxfilename) && strcmp(trxfilename,obj.trxfilename),
          success = true;
          return;
        end
        obj.trxfilename = trxfilename;
        [success,msg] = obj.UpdateStatusTable('trx');        
        % TODO: check that trx are parsable, remove bad experiments
      end
      
    end
    
    function [success,msg] = SetLabelFileName(obj,labelfilename)
      
      success = false;
      msg = '';

      if ischar(labelfilename),
        if ischar(obj.labelfilename) && strcmp(labelfilename,obj.labelfilename),
          success = true;
          return;
        end

        % reload labels from file
        for expi = 1:obj.nexps,
          [success1,msg] = obj.LoadLabelsFromFile(expi);
          if ~success1,
            return;
          end
        end
             
        % TODO: remove bad experiments
        
        obj.labelfilename = labelfilename;
        [success,msg] = obj.UpdateStatusTable('label');   
        
      end
      
    end
    
    function [success,msg] = LoadLabelsFromFile(obj,expi)
      
      success = false;
      msg = '';
      
      labelfilename = obj.GetFile('label',expi);
      if exist(labelfilename,'file'),

        obj.SetStatus('Loading labels for %s',obj.expdirs{expi});

        try
          loadedlabels = load(labelfilename,'t0s','t1s','names','flies','off');
          obj.labels(expi).t0s = loadedlabels.t0s;
          obj.labels(expi).t1s = loadedlabels.t1s;
          obj.labels(expi).names = loadedlabels.names;
          obj.labels(expi).flies = loadedlabels.flies;
          obj.labels(expi).off = loadedlabels.off;
        catch ME,
          msg = getReport(ME);
          obj.ClearStatus();
          return;
        end
        
        obj.ClearStatus();
        
      else
        
        obj.labels(expi).t0s = {};
        obj.labels(expi).t1s = {};
        obj.labels(expi).names = {};
        obj.labels(expi).flies = [];
        obj.labels(expi).off = [];
              
      end
      
      % TODO: update windowdata
      
      success = true;
      
    end
    
%     function AddWindowDataLabeled(obj,expi,flies)
%       
%       if numel(expi) ~= 1,
%         error('expi must be a scalar');
%       end
%       
%       if nargin < 3,
%         flies = (1:obj.nflies_per_exp(expi))';
%       end
%       
%       for i = 1:size(flies,1),
%         flies_curr = flies(i,:);
%         
%         % labels for this experiment, fly
%         labels_curr = obj.GetLabels(expi,flies_curr);
%         
%         % no labels?
%         if isempty(labels_curr.t0s),
%           continue;
%         end
%         
%         % loop through all behaviors
%         n = max(labels_curr.t1s);
%         idx = zeros(1,n);
%         for j = 1:obj.nbehaviors,
%           for k = find(strcmp(labels_curr.names,obj.labelnames{j})),
%             t0 = labels_curr.t0s(k);
%             t1 = labels_curr.t1s(k);
%             idx(t0+obj.labels_bin_off:t1+obj.labels_bin_off,j) = j;
%           end
%         end
%         m = nnz(idx);
%         if obj.expi == expi && all(flies_curr == obj.flies),
%           obj.windowdata_labeled(:,end+1:end+m) = obj.windowdata_curr(:,idx~=0);
%         else
%           windowfilenames = obj.GetFile('window',expi);
%           % TODO: make this work for multiple flies
%           windowdata_curr = load(windowfilenames{flies_curr(1)});
%           obj.windowdata_labeled(:,end+1:end+m) = windowdata_curr(:,idx~=0);
%         end
%         obj.exp_labeled(end+1:end+m) = expi;
%         obj.flies_labeled(end+1:end+m,:) = repmat(flies_curr,[m,1]);
%         obj.isintrainingset(end+1:end+m) = false;
%         obj.labelidx_labeled(end+1:end+m) = idx;
%         
%       end
%       
%     end


%     function UpdateWindowDataLabeled(obj)
% 
%       % indices into cached data for current experiment and flies
%       idxcurr = obj.exp_labeled' == obj.expi & all(bsxfun(@eq,obj.flies_labeled,obj.flies),2);
%       
%       % frames of current experiment and flies that have old labeled data
%       
%       % indices into cached data that are for this exp, these flies, have
%       % labelidx_old
%       idxcurr1 = idxcurr & obj.labelidx_old_labeled ~= 0;
%       % which frames for expi, flies that have labelidx_old ~= 0
%       tsold = obj.ts_labeled(idxcurr1);
%       % indices into labels_bin that have labelidx_old ~= 0
%       idxold = tsold+obj.labels_bin_off;
%             
%       % keep/add windowdata if labelidx_old ~= 0 or if labels_bin ~= 0
%       % indices into labels_bin for which new labelidx_new ~= 0
%       cacheidx = any(obj.labels_bin,2);
%       % or labelidx_old ~= 0
%       cacheidx(idxold) = true;
%       m = nnz(cacheidx);
% 
%       % labelidx_old for these frames
%       labelidx_old = zeros(1,m);
%       labelidx_old(idxold) = obj.labelidx_old_labeled(idxcurr1);
%       
%       % remove all data from the cache for the current exp, flies
%       obj.windowdata_labeled(:,idxcurr) = [];
%       obj.exp_labeled(idxcurr) = [];
%       obj.flies_labeled(idxcurr,:) = [];
%       obj.labelidx_old_labeled(idx) = [];
%       obj.labelidx_new_labeled(idx) = [];
%       obj.ts_labeled(idx) = [];
%       
%       % convert labels_bin to integer
%       n = size(obj.labels_bin,1);
%       labelidx = zeros(1,n);
%       for i = 1:obj.nbehaviors,
%         labelidx(obj.labels_bin(:,i)) = i;
%       end
%       
%       % add this data      
%       obj.windowdata_labeled(:,end+1:end+m) = obj.windowdata_curr(:,cacheidx);
%       obj.exp_labeled(end+1:end+m) = obj.expi;
%       obj.flies_labeled(end+1:end+m,:) = repmat(obj.flies,[m,1]);
%       obj.labelidx_old_labeled(end+1:end+m) = labelidx_old;
%       obj.labelidx_new_labeled(end+1:end+m) = labelidx;
%       obj.ts_labeled(end+1:end+m) = find(cacheidx) + obj.labels_bin_off;
% 
%     end
    
%     function RemoveFromWindowDataLabeled(obj,expi,flies)
% 
%       idx = ismember(obj.exp_labeled,expi);
%       if nargin >= 3,
%         idx = idx & ismember(obj.flies_labeled,flies,'rows');
%       end
%       obj.windowdata_labeled(:,idx) = [];
%       obj.exp_labeled(idx) = [];
%       obj.flies_labeled(idx,:) = [];
%       obj.isintrainingset(idx) = [];
%       obj.labelidx_labeled(idx) = [];
%       
%     end

    
    function [success,msg] = SetPerFrameDir(obj,perframedir)
      
      success = false;
      msg = '';

      if ischar(perframedir),
        if ischar(obj.perframedir) && strcmp(perframedir,obj.perframedir),
          success = true;
          return;
        end

        obj.perframedir = perframedir;
        
        % TODO: check per-frame directories are okay, remove bad
        % experiments
        
        [success,msg] = obj.UpdateStatusTable('label');
      end
      
    end
    
%     function [success,msg] = SetWindowFileName(obj,windowfilename)
%       
%       success = false;
%       msg = '';
% 
%       if ischar(windowfilename),
%         obj.windowfilename = windowfilename;
%         if ischar(obj.windowfilename) && strcmp(windowfilename,obj.windowfilename),
%           success = true;
%           return;
%         end
% 
%         % TODO: check window data for existing experiments, remove bad experiments, update
%         % windowdata_labeled, etc. 
%         
%         [success,msg] = obj.UpdateStatusTable('window');
%       end
%       
%     end
    
    function [success,msg] = SetDefaultPath(obj,defaultpath)
      
      success = false;
      msg = '';
      
      if ischar(defaultpath),
        
        if ~isempty(defaultpath) && ~exist(defaultpath,'file'),
          msg = sprintf('defaultpath directory %s does not exist',defaultpath);
          return;
        end
          
        obj.defaultpath = defaultpath;
        success = true;
      end

    end
    
    function [success,msg] = SetRootOutputDir(obj,rootoutputdir)
      
      success = true;
      msg = '';
      if ischar(rootoutputdir),
        if ischar(obj.rootoutputdir) && strcmp(obj.rootoutputdir,rootoutputdir),
          success = true;
          return;
        end
        if ~exist(rootoutputdir,'file'),
          msg = sprintf('root output directory %s does not exist',rootoutputdir);
          success = false;
          return;
        end
        obj.rootoutputdir = rootoutputdir;
        for i = 1:obj.nexps,
          obj.outexpdirs{i} = fullfile(rootoutputdir,obj.expnames{i});
        end
        % TODO: check all files are okay, remove bad experiments
        
        [success,msg] = obj.UpdateStatusTable();
      end
      
    end    
    
    function [success,msg] = SetClassifierFileName(obj,classifierfilename)

      success = false;
      msg = '';
      
      if ~isempty(classifierfilename) && exist(classifierfilename,'file'),
        try
          obj.SetStatus('Loading classifier from %s',obj.classifierfilename);
          loadeddata = load(obj.classifierfilename,obj.classifiervars{:});

          % remove all experiments
          obj.RemoveExpDirs(1:obj.nexps);
          
          % set config file
          if ~strcmp(obj.configfilename,'configfilename'),
            obj.SetConfigFileName(loadeddata.configfilename);
          end

          % set movie
          [success,msg] = obj.SetMovieFileName(loadeddata.moviefilename);
          if ~success,
            error(msg);
          end

          % trx
          [success,msg] = obj.SetTrxFileName(loadeddata.trxfilename);
          if ~success,
            error(msg);
          end
      
          % label
          [success,msg] = obj.SetLabelFileName(loadeddata.labelfilename);
          if ~success,
            error(msg);
          end
      
          % perframedir
          [success,msg] = obj.SetPerFrameDir(loadeddata.perframedir);
          if ~success,
            error(msg);
          end

          % featureparamsfilename
          [success,msg] = obj.SetFeatureParamsFileName(loadeddata.featureparamsfile);
          if ~success,
            error(msg);
          end
      
          % rootoutputdir
          [success,msg] = obj.SetRootOutputDir(loadeddata.rootoutputdir);
          if ~success,
            error(msg);
          end
           
          % set experiment directories
          obj.SetExpDirs(loadeddata.expdirs,loadeddata.outexpdirs,loadeddata.expnames,...
            loadeddata.nflies_per_exp,loadeddata.firstframes_per_exp,loadeddata.endframes_per_exp); 

          [success,msg] = obj.UpdateStatusTable();
          if ~success,
            error(msg);
          end
          
          % update cached data
          obj.windowdata = struct('X',[],'exp',[],'flies',[],'t',[],...
            'labelidx_old',[],'labelidx_new',[],'featurenames',{{}},...
            'predicted',[],'predicted_probs',[],'isvalidprediction',[]);
          [success,msg] = obj.PreLoadLabeledData();
          if ~success,
            error(msg);
          end          
          
          obj.classifier = loadeddata.classifier;
          obj.classifiertype = loadeddata.classifiertype;
          obj.classifier_params = loadeddata.classifier_params;
          
        catch ME,
          errordlg(getReport(ME),'Error loading classifier from file');
        end
        
        obj.ClearStatus();
        
        obj.classifierfilename = classifierfilename;
        
      end

    end

    function [success,msg] = PreLoadLabeledData(obj)

      success = false;
      msg = '';
      
      for expi = 1:obj.nexps,
        for i = 1:size(obj.labels(expi).flies,1),
          flies = obj.labels(expi).flies(i,:);
          labels_curr = obj.GetLabels(expi,flies);
          ts = [];
          for j = 1:numel(labels_curr.t0s),
            ts = [ts,labels_curr.t0s(j):labels_curr.t1s(j)]; %#ok<AGROW>
          end
          [success1,msg] = obj.PreLoadWindowData(expi,flies,ts);
          if ~success1,
            return;
          end            
        end
      end
      success = true;
      
    end
    
    function SaveClassifier(obj)
      
      try
        s = struct;
        for i = 1:numel(obj.classifiervars),
          fn = obj.classifiervars{i};
          s.(fn) = obj.(fn);
        end
        save(obj.classifierfilename,'-struct','s');
      catch ME,
        errordlg(getReport(ME),'Error saving classifier to file');
      end      
      
    end
    
    function [success,msg] = SetExpDirs(obj,expdirs,outexpdirs,nflies_per_exp,firstframes_per_exp,endframes_per_exp)

      success = false;
      msg = '';
      
      if isnumeric(expdirs),
        return;
      end
      
      if nargin < 2,
        error('Usage: obj.SetExpDirs(expdirs,[outexpdirs],[nflies_per_exp])');
      end
      
      isoutexpdirs = nargin > 2 && ~isnumeric(outexpdirs);
      isnflies = nargin > 3 && ~isempty(nflies_per_exp);
      isfirstframes = nargin > 4 && ~isempty(firstframes_per_exp);
      isendframes = nargin > 5 && ~isempty(endframes_per_exp);
      
      % check inputs
      
      % sizes must match
      if isoutexpdirs && numel(expdirs) ~= numel(outexpdirs),
        error('expdirs and outexpdirs do not match size');
      end
      if isnflies && numel(expdirs) ~= numel(nflies_per_exp),
        error('expdirs and nflies_per_exp do not match size');
      end
      
      oldexpdirs = obj.expdirs;
      
      % remove oldexpdirs
      [success1,msg] = obj.RemoveExpDirs(setdiff(oldexpdirs,expdirs));
      if ~success1,
        return;
      end

      % add new expdirs
      idx = find(~ismember(expdirs,oldexpdirs));
      success = true;
      for i = idx,
        params = cell(1,nargin-1);
        params{1} = expdirs{i};
        if isoutexpdirs,
          params{2} = outexpdirs{i};
        end
        if isnflies,
          params{3} = nflies_per_exp(i);
        end
        if isfirstframes,
          params{4} = firstframes_per_exp{i};
        end
        if isendframes,
          params{5} = endframes_per_exp{i};
        end
        [success1,msg1] = obj.AddExpDir(params{:});
        success = success && success1;
        if isempty(msg),
          msg = msg1;
        else
          msg = sprintf('%s\n%s',msg,msg1);
        end
      end

    end
    
    function [success,msg] = GetTrxInfo(obj,expi)
      success = true;
      msg = '';
      obj.SetStatus('Reading trx info for experiment %s',obj.expdirs{expi});
      if numel(obj.nflies_per_exp) < expi || ...
          numel(obj.firstframes_per_exp) < expi || ...
          numel(obj.endframes_per_exp) < expi || ...
          isnan(obj.nflies_per_exp(expi)) || ...
          isnan(obj.firstframes_per_exp(expi)) || ...
          isnan(obj.endframes_per_exp(expi)),
        trxfile = fullfile(obj.expdirs{expi},obj.GetFileName('trx'));
        if ~exist(trxfile,'file'),
          msg = sprintf('Trx file %s does not exist, cannot count flies',trxfile);
          success = false;
        else
          try
            % REMOVE THIS
            global CACHED_TRX; %#ok<TLEV>
            if isempty(CACHED_TRX),
              hwait = mywaitbar(0,sprintf('Loading trx to determine number of flies for %s',obj.expnames{expi}));
              trx = load_tracks(trxfile);
              if ishandle(hwait), delete(hwait); end
              CACHED_TRX = trx;
            else
              trx = CACHED_TRX;
            end
          catch ME,
            msg = sprintf('Could not load trx file for experiment %s to count flies: %s',obj.expdirs{expi},getReport(ME));
          end
          obj.nflies_per_exp(expi) = numel(trx);
          obj.firstframes_per_exp{expi} = [trx.firstframe];
          obj.endframes_per_exp{expi} = [trx.endframe];
        end
      end
      
      obj.ClearStatus();
      
    end
    
    function [success,msg] = AddExpDir(obj,expdir,outexpdir,nflies_per_exp,firstframes_per_exp,endframes_per_exp)

      success = false;
      msg = '';
      
      if isnumeric(expdir),
        return;
      end
      
      if nargin < 2,
        error('Usage: obj.AddExpDirs(expdir,[outexpdir],[nflies_per_exp])');
      end

      % make sure directory exists
      if ~exist(expdir,'file'),
        msg = sprintf('expdir %s does not exist',expdir);
        return;
      end
      
      isoutexpdir = nargin > 2 && ~isnumeric(outexpdir);
      istrxinfo = nargin > 5 && ~isempty(nflies_per_exp);

      % base name
      [~,expname] = myfileparts(expdir);
      
      % expnames and rootoutputdir must match
      if isoutexpdir,
        [rootoutputdir,outname] = myfileparts(outexpdir); %#ok<*PROP>
        if ~strcmp(expname,outname),
          msg = sprintf('expdir and outexpdir do not match base names: %s ~= %s',expname,outname);
          return;
        end
        if ischar(obj.rootoutputdir) && ~strcmp(rootoutputdir,obj.rootoutputdir),
          msg = sprintf('Inconsistent root output directory: %s ~= %s',rootoutputdir,obj.rootoutputdir);
          return;
        end
      elseif ~ischar(obj.rootoutputdir),
        msg = 'rootoutputdir not yet set';
        return;
      else
        rootoutputdir = obj.rootoutputdir;        
      end

      if ~isoutexpdir,
        outexpdir = fullfile(rootoutputdir,expname);
      end
      
      % create missing outexpdirs
      if ~exist(outexpdir,'file'),
        [success1,msg1] = mkdir(rootoutputdir,expname);
        if ~success1,
          msg = (sprintf('Could not create output directory %s, failed to set expdirs: %s',outexpdir,msg1));
          return;
        end
      end
            
      % okay, checks succeeded, start storing stuff
      obj.nexps = obj.nexps + 1;
      obj.expdirs{end+1} = expdir;
      obj.expnames{end+1} = expname;
      obj.rootoutputdir = rootoutputdir;
      obj.outexpdirs{end+1} = outexpdir;
      
      % get trxinfo
      if istrxinfo,
        obj.nflies_per_exp(end+1) = nflies_per_exp;
        obj.firstframes_per_exp(end+1) = firstframes_per_exp;
        obj.endframes_per_exp(end+1) = endframes_per_exp;
      else
        obj.nflies_per_exp(end+1) = nan;
        obj.firstframes_per_exp{end+1} = [];
        obj.endframes_per_exp{end+1} = [];
        [success1,msg1] = obj.GetTrxInfo(obj.nexps);
        if ~success1,
          msg = sprintf('Error getting basic trx info: %s',msg1);
          obj.RemoveExpDirs(obj.nexps);
          return;
        end
      end
      
      % load labels for this experiment
      [success1,msg] = obj.LoadLabelsFromFile(obj.nexps);
      if ~success1,
        obj.RemoveExpDirs(obj.nexps);
        return;
      end
      
      [success1,msg1] = obj.UpdateStatusTable('',obj.nexps);
      if ~success1,
        msg = msg1;
        obj.RemoveExpDirs(obj.nexps);
        return;
      end
      
      % save default path
      obj.defaultpath = expdir;
      
      success = true;
      
    end
   
    function [success,msg] = RemoveExpDirs(obj,expi)
      
      success = false;
      msg = '';
      
      if any(obj.nexps < expi) || any(expi < 1),
        msg = sprintf('expi = %s must be in the range 1 < expi < nexps = %d',mat2str(expi),obj.nexps);
        return;
      end

      obj.expdirs(expi) = [];
      obj.expnames(expi) = [];
      obj.outexpdirs(expi) = [];
      obj.nflies_per_exp(expi) = [];
      obj.firstframes_per_exp(expi) = [];
      obj.endframes_per_exp(expi) = [];
      obj.nexps = obj.nexps - numel(expi);
      obj.labels(expi) = [];
      % TODO: exp2labeloff

      % update current exp, flies
      if ~isempty(obj.expi) && obj.expi > 0 && ismember(obj.expi,expi),
        
        % change to different experiment, by default choose fly 1
        % TODO: allow for more than one fly to be selected at once
        obj.expi = 0;
        obj.flies = nan(size(obj.flies));
        % TODO: may want to save labels somewhere before just overwriting
        % labelidx
        if obj.nexps > 0,
          obj.PreLoad(obj.nexps,1);
        end

      end
      
      % TODO: windowdata_labeled, etc
      
      success = true;
      
    end
    
    function res = GetFileName(obj,file)
      switch file,
        case 'movie',
          res = obj.moviefilename;
        case 'trx',
          res = obj.trxfilename;
        case 'label',
          res = obj.labelfilename;
%         case 'window',
%           res = obj.windowfilename;
        case 'perframedir',
          res = obj.perframedir;
        otherwise
          error('Unknown file type %s',file);
      end
    end
    
    function [filename,timestamp] = GetFile(obj,file,expi)
      
      % base name
      fn = obj.GetFileName(file);
      
      % if this is an output file, only look in output experiment directory
      if JLabelData.IsOutputFile(file),
        expdirs_try = obj.outexpdirs(expi);
      else
        % otherwise, first look in output directory, then look in input
        % directory
        expdirs_try = {obj.outexpdirs{expi},obj.expdirs{expi}};
      end
      
      % initialize timestamp = -inf if we never set
      timestamp = -inf;

      % loop through directories to look in
      for j = 1:numel(expdirs_try),
        expdir = expdirs_try{j}; 
        
        % are there per-fly files?
        if JLabelData.IsPerFlyFile(file),
          
          % if per-fly, then there will be one file per fly
          filename = cell(1,obj.nflies_per_exp(expi));
          [~,name,ext] = fileparts(fn);
          file_exists = true;
          for fly = 1:obj.nflies_per_exp(expi),
            filename{fly} = fullfile(expdir,sprintf('%s_fly%02d%s',name,fly,ext));
          end
          
          % check this directory, get timestamp
          timestamp = -inf;
          for fly = 1:obj.nflies_per_exp(expi),
            % doesn't exist? then just return timestamp = -inf
            if ~exist(filename{fly},'file');
              file_exists = false;
              timestamp = -inf;
              break;
            end
            tmp = dir(filename{fly});
            timestamp = max(tmp.datenum,timestamp);
          end
          
          % file exists? then don't search next directory
          if file_exists,
            break;
          end
          
        else
          
          % just one file to look for
          filename = fullfile(expdir,fn);
          if exist(filename,'file'),
            tmp = dir(filename);
            timestamp = tmp.datenum;
            break;
          end
          
        end
        
      end
      
    end
    
    
    function [success,msg] = GenerateMissingFiles(obj,expi)
      
      success = true;
      msg = '';
      
      for i = 1:numel(obj.filetypes),
        file = obj.filetypes{i};
        if obj.IsRequiredFile(file) && obj.CanGenerateFile(file) && ...
            ~obj.FileExists(file,expi),
          fprintf('Generating %s for %s...\n',file,obj.expnames{expi});
          switch file,
%             case 'window',
%               [success1,msg1] = obj.GenerateWindowFeaturesFiles(expi);
%               success = success && success1;
%               if ~success1,
%                 msg = [msg,'\n',msg1]; %#ok<AGROW>
%               end
            case 'perframedir',
          end
        end
      end
      [success1,msg1] = obj.UpdateStatusTable();
      success = success && success1;
      if isempty(msg),
        msg = msg1;
      else
        msg = sprintf('%s\n%s',msg,msg1);
      end
      
    end
    
%     function [success,msg] = GenerateWindowFeaturesFiles(obj,expi,doforce)
%       
%       success = false;
%       msg = '';
%       
%       if ~exist('doforce','var'),
%         doforce = false;
%       end
% 
%       hwait = mywaitbar(0,'Computing window features...');
% 
%       filenames = obj.GetFile('window',expi);
%       perframedir = obj.GetFile('perframedir',expi);
%       for fly = 1:obj.nflies_per_exp(expi),
%         filename = filenames{fly};
%         if ~doforce && exist(filename,'file'),
%           fprintf('File %s exists, skipping\n',filename);
%           continue;
%         end
%         try
%         X = [];
%         feature_names = {};
%         for j = 1:numel(obj.perframefns),
%           fn = obj.perframefns{j};
%           perframedata = load(fullfile(perframedir,[fn,'.mat']));
%           hwait = mywaitbar((fly-1+(j-1)/numel(obj.perframefns))/obj.nflies_per_exp(expi),hwait,...
%             sprintf('Computing %s window features (%d/%d) for fly %d/%d',fn,j,numel(obj.perframefns),fly,numel(perframedata.data)));
%           [x_curr,feature_names_curr] = ...
%             ComputeWindowFeatures(perframedata.data{fly},obj.windowfeaturescellparams.(fn){:});
%           nold = size(X,2);
%           nnew = size(x_curr,2);
%           if nold > nnew,
%             x_curr(:,end+1:end+nold-nnew) = nan;
%           elseif nnew > nold,
%             X(:,end+1:end+nnew-nold) = nan;
%           end
%           X = [X;x_curr]; %#ok<AGROW>
%           feature_names = [feature_names,cellfun(@(s) [{fn},s],feature_names_curr,'UniformOutput',false)]; %#ok<AGROW>
%         end
%         hwait = mywaitbar(fly/obj.nflies_per_exp(expi),hwait,sprintf('Saving window features for fly %d/%d to file...',fly,obj.nflies_per_exp(expi)));
%         save(filename,'X','feature_names');
%         catch ME,
%           msg = getReport(ME);
%           return;
%         end
%       end
%       if exist('hwait','var') && ishandle(hwait),
%         delete(hwait);
%       end
%       
%     end
    
    function [success,msg] = SetFeatureParamsFileName(obj,featureparamsfilename)
      success = false;
      msg = '';
      
      if ischar(obj.featureparamsfilename) && strcmp(featureparamsfilename,obj.featureparamsfilename),
        success = true;
        return;
      end
      
      if obj.nexps > 0,
        msg = 'Currently, feature params file can only be changed when no experiments are loaded';
        return;
      end
      try
        [windowfeaturesparams,windowfeaturescellparams] = ...
          ReadPerFrameParams(featureparamsfilename); %#ok<PROP>
      catch ME,
        msg = sprintf('Error reading feature parameters file %s: %s',...
          params.featureparamsfilename,getReport(ME));
        return;
      end
      obj.windowfeaturesparams = windowfeaturesparams; %#ok<PROP>
      obj.windowfeaturescellparams = windowfeaturescellparams; %#ok<PROP>
      obj.featureparamsfilename = featureparamsfilename;
      obj.perframefns = fieldnames(obj.windowfeaturescellparams);
      if numel(obj.perframedata) ~= numel(obj.perframefns),
        obj.perframedata = cell(1,numel(obj.perframefns));
      end
      success = true;
    end
    
    function [success,msg] = UpdateStatusTable(obj,filetypes,expis)

      msg = '';
      success = false;

      if nargin > 1 && ~isempty(filetypes),
        [ism,fileis] = ismember(filetypes,obj.filetypes);
        if any(~ism),
          msg = 'Unknown filetypes';
          return;
        end
      else
        fileis = 1:numel(obj.filetypes);
      end
      if nargin <= 2 || isempty(expis),
        expis = 1:obj.nexps;
      end
      
      % initialize fileexists table
      obj.fileexists(expis,fileis) = false;
      obj.filetimestamps(expis,fileis) = nan;
      
      
      % loop through all file types
      for filei = fileis,
        file = obj.filetypes{filei};
        % loop through experiments
        for expi = expis,

          % check for existence of current file(s)
          [fn,obj.filetimestamps(expi,filei)] = obj.GetFile(file,expi);
          if iscell(fn),
            obj.fileexists(expi,filei) = all(cellfun(@(s) exist(s,'file'),fn));
          else
            obj.fileexists(expi,filei) = exist(fn,'file');
          end
                    
        end
      end

      % store old values to see if latest change broke something
      old_filesfixable = obj.filesfixable;
      old_allfilesexist = obj.allfilesexist;

      % initialize summaries to true
      obj.filesfixable = true;
      obj.allfilesexist = true;

      for filei = 1:numel(obj.filetypes),
        file = obj.filetypes{filei};
        % loop through experiments
        for expi = 1:obj.nexps,
          
          % if file doesn't exist and is required, then not all files exist
          if ~obj.fileexists(expi,filei),
            if JLabelData.IsRequiredFile(file),
              obj.allfilesexist = false;
              % if furthermore file can't be generated, then not fixable
              if ~JLabelData.CanGenerateFile(file),
                obj.filesfixable = false;                
                msg1 = sprintf('%s missing and cannot be generated.',fn);
                if isempty(msg),
                  msg = msg1;
                else
                  msg = sprintf('%s\n%s',msg,msg1);
                end
              end
            end
          end
          
        end
      end

      % fail if was ok and now not ok
      success = ~(old_allfilesexist || old_filesfixable) || ...
        (obj.allfilesexist || obj.filesfixable);
      
    end
    
    function [fe,ft] = FileExists(obj,file,expi)
      filei = find(strcmpi(file,obj.filetypes),1);
      if isempty(filei),
        error('file type %s does not match any known file type',file);
      end
      if nargin < 3,
        expi = 1:obj.nexps;
      end
      fe = obj.fileexists(expi,filei);
      ft = obj.filetimestamps(expi,filei);
    end
    
    function trx = GetTrx(obj,expi,flies,ts)
      
      if numel(expi) ~= 1,
        error('expi must be a scalar');
      end
      
      if expi ~= obj.expi,
        [success,msg] = obj.LoadTrx(expi);
        if ~success,
          error('Error loading trx for experiment %d: %s',expi,msg);
        end
      end

      if nargin < 3,
        trx = obj.trx;
        return;
      end
      
      nflies = numel(flies);
      c = cell(1,nflies);
      trx = struct('x',c,'y',c,'a',c,'b',c,'theta',c);
      for i = 1:numel(flies),
        fly = flies(i);
        if nargin < 4,
          trx(i).x = obj.trx(fly).x;
          trx(i).y = obj.trx(fly).y;
          trx(i).a = obj.trx(fly).a;
          trx(i).b = obj.trx(fly).b;
          trx(i).theta = obj.trx(fly).theta;
        else
          js = min(obj.trx(fly).nframes,max(1,ts + obj.trx(fly).off));
          trx(i).x = obj.trx(fly).x(js);
          trx(i).y = obj.trx(fly).y(js);
          trx(i).a = obj.trx(fly).a(js);
          trx(i).b = obj.trx(fly).b(js);
          trx(i).theta = obj.trx(fly).theta(js);
        end
      end
    end

    function varargout = GetTrxX(obj,expi,flies,ts)
      
      if numel(expi) ~= 1,
        error('expi must be a scalar');
      end
      
      if expi ~= obj.expi,
        [success,msg] = obj.LoadTrx(expi);
        if ~success,
          error('Error loading trx for experiment %d: %s',expi,msg);
        end
      end

      if nargin < 3,
        varargout = {obj.trx.x};
        return;
      end
      
      if nargin < 4,
        varargout = {obj.trx(flies).x};
        return;
      end
      
      nflies = numel(flies);
      varargout = cell(1,nflies);
      for i = 1:numel(flies),
        fly = flies(i);
        js = min(obj.trx(fly).nframes,max(1,ts + obj.trx(fly).off));
        varargout{i} = obj.trx(fly).x(js);
      end
    end

    function varargout = GetTrxY(obj,expi,flies,ts)
      
      if numel(expi) ~= 1,
        error('expi must be a scalar');
      end
      
      if expi ~= obj.expi,
        [success,msg] = obj.LoadTrx(expi);
        if ~success,
          error('Error loading trx for experiment %d: %s',expi,msg);
        end
      end

      if nargin < 3,
        varargout = {obj.trx.x};
        return;
      end
      
      if nargin < 4,
        varargout = {obj.trx(flies).x};
        return;
      end
      
      nflies = numel(flies);
      varargout = cell(1,nflies);
      for i = 1:numel(flies),
        fly = flies(i);
        js = min(obj.trx(fly).nframes,max(1,ts + obj.trx(fly).off));
        varargout{i} = obj.trx(fly).x(js);
      end
    end

    
    function t0 = GetTrxFirstFrame(obj,expi,flies)

      if numel(expi) ~= 1,
        error('expi must be a scalar');
      end
      
      if expi ~= obj.expi,
        [success,msg] = obj.LoadTrx(expi);
        if ~success,
          error('Error loading trx for experiment %d: %s',expi,msg);
        end
      end

      if nargin < 3,
        t0 = [obj.trx.firstframe];
        return;
      end

      t0 = [obj.trx(flies).firstframe];
      
    end

    function t1 = GetTrxEndFrame(obj,expi,flies)

      if numel(expi) ~= 1,
        error('expi must be a scalar');
      end
      
      if expi ~= obj.expi,
        [success,msg] = obj.LoadTrx(expi);
        if ~success,
          error('Error loading trx for experiment %d: %s',expi,msg);
        end
      end

      if nargin < 3,
        t1 = [obj.trx.endframe];
        return;
      end

      t1 = [obj.trx(flies).endframe];
      
    end
    
    function [success,msg] = LoadTrx(obj,expi)

      success = false;
      msg = '';
      
      if numel(expi) ~= 1,
        error('expi must be a scalar');
      end

      if expi < 1,
        msg = 'expi not yet set';
        return;
      end
      
      % load trx
      try
        trxfilename = obj.GetFile('trx',expi);
  
        hwait = mywaitbar(0,'Loading trx');
  
        % TODO: remove this
        global CACHED_TRX; %#ok<TLEV>
        if isempty(CACHED_TRX),
          obj.trx = load_tracks(trxfilename);
          CACHED_TRX = obj.trx;
        else
          obj.trx = CACHED_TRX;
        end
      catch ME,
        msg = sprintf('Error loading trx from file %s: %s',trxfilename,getReport(ME));
        if ishandle(hwait),
          delete(hwait);
        end
        return;
      end

      if ishandle(hwait),
        delete(hwait);
      end
      success = true;
      
    end
    
    function [success,msg] = PreLoad(obj,expi,flies)
      
      success = false;
      msg = '';
      
      if numel(expi) ~= 1,
        error('expi must be a scalar');
      end

      if numel(unique(flies)) ~= numel(flies),
        msg = 'flies must all be unique';
        return;
      end
      
      diffexpi = expi ~= obj.expi;
      diffflies = diffexpi || ~isempty(setdiff(flies,obj.flies)) || ~isempty(setdiff(obj.flies,flies));

      if diffflies && ~isempty(obj.expi) && obj.expi > 0,
        % store labels currently in labelidx to labels
        obj.StoreLabels();
      end
      
      if diffexpi,
        [success1,msg] = obj.LoadTrx(expi);
        if ~success1,
          return;
        end
      end

      obj.expi = expi;
            
      if diffflies,
        
        % set labelidx from labels
        obj.CacheLabelIdx(flies);

        % load perframedata
        obj.SetStatus('Loading per-frame data for %s, flies %s',obj.expdirs{expi},mat2str(obj.flies));
        perframedir = obj.GetFile('perframedir',obj.expi);
        for j = 1:numel(obj.perframefns),
          fn = obj.perframefns{j};
          file = fullfile(perframedir,[fn,'.mat']);
          if ~exist(file,'file'),
            msg = sprintf('Per-frame data file %s does not exist',file);
            return;
          end
          try
            tmp = load(file);
            obj.perframedata{j} = tmp.data{flies(1)};
          catch ME,
            msg = getReport(ME);
          end
          obj.UpdatePredictedIdx();
        end
        obj.ClearStatus();
        
%         % window data for the currenf flies
%         [success1,msg] = obj.LoadWindowData(expi,flies);
%         if ~success1,
%           return;
%         end
        
      end
      obj.flies = flies;
           
      success = true;
      
    end

%     function [success,msg] = LoadWindowData(obj,expi,flies)
% 
%       success = false;
%       msg = '';
%       
%       windowfilenames = obj.GetFile('window',expi);
%       % TODO: make this work for multiple flies
%       try
%         obj.windowdata_curr = load(windowfilenames{flies(1)});
%       catch ME,
%         msg = getReport(ME);
%         return;
%       end
%       
%       success = true;
%       
%     end

    
    function CacheLabelIdx(obj,flies)
      
      if nargin < 2,
        flies = obj.flies;
      end
      
      obj.t0_curr = max(obj.GetTrxFirstFrame(obj.expi,flies));
      obj.t1_curr = min(obj.GetTrxEndFrame(obj.expi,flies));

      n = obj.t1_curr-obj.t0_curr+1;
      obj.labelidx = zeros(1,n);

      % get labels for current flies 
      labels_curr = obj.GetLabels(obj.expi,flies);
      
      % loop through all behaviors
      obj.labelidx_off = labels_curr.off;
      for j = 1:obj.nbehaviors,
        for k = find(strcmp(labels_curr.names,obj.labelnames{j})),
          t0 = labels_curr.t0s(k);
          t1 = labels_curr.t1s(k);
          obj.labelidx(t0+obj.labelidx_off:t1+obj.labelidx_off) = j;
        end
      end
      
    end
    
    function [labelidx,T0,T1] = GetLabelIdx(obj,expi,flies)
      
      T0 = max(obj.GetTrxFirstFrame(expi,flies));
      T1 = min(obj.GetTrxEndFrame(expi,flies));
      n = T1-T0+1;
      off = 1 - T0;
      labels_curr = obj.GetLabels(expi,flies);
      labelidx = zeros(1,n);
      for i = 1:obj.nbehaviors,
        for j = find(strcmp(labels_curr.names,obj.labelnames{i})),
          t0 = labels_curr.t0s(j);
          t1 = labels_curr.t1s(j);
          labelidx(t0+off:t1+off) = i;
        end
      end
      
    end
    
    function labels_curr = GetLabels(obj,expi,flies)

      labels_curr = struct('t0s',[],'t1s',[],'names',{{}},'off',0);
      
      if nargin < 2 || isempty(expi),
        expi = obj.expi;
      end
      
      if nargin < 3 || isempty(flies),
        flies = obj.flies;
      end
      
      [ism,fliesi] = ismember(flies,obj.labels(expi).flies,'rows');
      if ism,
        labels_curr.t0s = obj.labels(expi).t0s{fliesi};
        labels_curr.t1s = obj.labels(expi).t1s{fliesi};
        labels_curr.names = obj.labels(expi).names{fliesi};
        labels_curr.off = obj.labels(expi).off(fliesi);
      else
        if expi ~= obj.expi,
          error('This should never happen -- only should get new labels for current experiment');
        end
        t0_curr = max(obj.GetTrxFirstFrame(expi,flies));
        labels_curr.off = 1-t0_curr;
      end

      
    end
    
    function StoreLabels(obj)
      
      % flies not yet initialized
      if isempty(obj.flies) || all(isnan(obj.flies)) || isempty(obj.labelidx),
        return;
      end
      
      % update labels
      newlabels = struct('t0s',[],'t1s',[],'names',{{}},'flies',[],'off',obj.labelidx_off);
      for j = 1:obj.nbehaviors,
        [i0s,i1s] = get_interval_ends(obj.labelidx==j);
        if ~isempty(i0s),
          n = numel(i0s);
          newlabels.t0s(end+1:end+n) = i0s - obj.labelidx_off;
          newlabels.t1s(end+1:end+n) = i1s - obj.labelidx_off;
          newlabels.names(end+1:end+n) = repmat(obj.labelnames(j),[1,n]);
        end
      end
      [ism,j] = ismember(obj.flies,obj.labels(obj.expi).flies,'rows');
      if ~ism,
        j = size(obj.labels(obj.expi).flies,1)+1;
      end
      obj.labels(obj.expi).t0s{j} = newlabels.t0s;
      obj.labels(obj.expi).t1s{j} = newlabels.t1s;
      obj.labels(obj.expi).names{j} = newlabels.names;
      obj.labels(obj.expi).flies(j,:) = obj.flies;
      obj.labels(obj.expi).off(j) = obj.labelidx_off;
      
      ts = find(obj.labelidx~=0) - obj.labelidx_off;
      [success,msg] = obj.PreLoadWindowData(obj.expi,obj.flies,ts);
      if ~success,
        warning(msg);
      end

      if ~isempty(obj.windowdata.exp),
        idxcurr = obj.windowdata.exp == obj.expi & ...
          all(bsxfun(@eq,obj.windowdata.flies,obj.flies),2);
        obj.windowdata.labelidx_new(idxcurr) = obj.labelidx(obj.windowdata.t(idxcurr)+obj.labelidx_off);
      end
      
      %obj.UpdateWindowDataLabeled(obj.expi,obj.flies);
      
    end
    
    function Train(obj)
      
      % load all labeled data
      success = true;
      for expi = 1:obj.nexps,
        for i = 1:numel(obj.labels(expi).t0s),
          flies = obj.labels(expi).flies(i,:);
          ts = [];
          for j = 1:numel(obj.labels(expi).t0s{i}),
            ts = [ts,obj.labels(expi).t0s{i}(j):obj.labels(expi).t1s{i}(j)]; %#ok<AGROW>
          end
          ts = unique(ts);
          [success,msg] = obj.PreLoadWindowData(expi,flies,ts);
          if ~success,
            break;
          end
        end
      end
      if ~success,
        warning(msg);
        return;
      end

      islabeled = obj.windowdata.labelidx_new ~= 0;
      if ~any(islabeled),
        return;
      end
      
      switch obj.classifiertype,
      
        case 'ferns',
          if isempty(obj.classifier),
            
            % train classifier
            obj.SetStatus('Training fern classifier from %d examples...',numel(islabeled));

            s = struct2paramscell(obj.classifier_params);
            obj.classifier = fernsClfTrain( obj.windowdata.X(islabeled,:), obj.windowdata.labelidx_new(islabeled), s{:} );
            obj.windowdata.labelidx_old = obj.windowdata.labelidx_new;
                        
          else
            
            % new data added to windowdata at the end, so classifier.inds still
            % matches windowdata(:,1:Nprev)
            Nprev = size(obj.windowdata.labelidx_old);
            Ncurr = size(obj.windowdata.labelidx_new);
            waslabeled = obj.windowdata.labelidx_old(1:Nprev) ~= 0;
            islabeled = obj.windowdata.labelidx_new(1:Nprev) ~= 0;
            
            % replace labels for examples that have been relabeled:
            % islabeled & waslabeled will not change
            idx_relabel = islabeled & waslabeled & (obj.windowdata.labelidx_new(1:Nprev) ~= obj.windowdata.labelidx_old(1:Nprev));
            if any(idx_relabel),
              obj.SetStatus('Updating fern classifier for %d relabeled examples...',nnz(idx_relabel));
              [obj.classifier] = fernsClfRelabelTrainingData( obj.windowdata.labelidx_old(waslabeled), ...
                obj.windowdata.labelidx_new(waslabeled), obj.classifier );
              % update labelidx_old
              obj.windowdata.labelidx_old(idx_relabel) = obj.windowdata.labelidx_new(idx_relabel);
            end
            
            % remove training examples that were labeled but now aren't
            idx_remove = waslabeled & ~islabeled(1:Nprev);
            if any(idx_remove),
              obj.SetStatus('Removing %d training examples from fern classifier',nnz(idx_remove));
              [obj.classifier] = fernsClfRemoveTrainingData(obj.windowdata.labelidx_old(waslabeled), idx_remove(waslabeled), obj.classifier );
              % update labelidx_old
              obj.windowdata.labelidx_old(idx_remove) = 0;
            end
            % update islabeled and waslabeled
            islabeled = obj.windowdata.labelidx_new ~= 0;
            waslabeled = [obj.windowdata.labelidx_old ~= 0;false(Ncurr-Nprev,1)];
            % now only examples with islabeled should be in training set
            
            % add training examples that are labeled now but weren't before
            idx_add = ~waslabeled(islabeled);
            if any(idx_add),
              obj.SetStatus('Adding %d new examples to fern classifier...',nnz(idx_add));
              [obj.classifier] = fernsClfAddTrainingData( obj.windowdata.X(islabeled,:), ...
                obj.windowdata.labelidx_new(islabeled), find(idx_add), obj.classifier );
              % update labelidx_old
              obj.windowdata.labelidx_old(~waslabeled&islabeled) = ...
                obj.windowdata.labelidx_new(~waslabeled&islabeled);
            end
            
            % labelidx_old and new should match
            if ~all(obj.windowdata.labelidx_old == obj.windowdata.labelidx_new),
              error('Sanity check: labelidx_old and labelidx_new should match');
            end
          end
      end

      obj.ClearStatus();
      
      % all predictions invalid now
      obj.windowdata.isvalidprediction(:) = false;
      
      % predict for all window data
      obj.PredictLoaded();
      
    end
    
    function PredictLoaded(obj)
      
      if isempty(obj.classifier),
        return;
      end
      
      % apply classifier
      switch obj.classifiertype,
        
        case 'ferns',
          obj.SetStatus('Applying fern classifier to %d windows',size(obj.windowdata.X,1));
          [obj.windowdata.predicted,...
            obj.windowdata.predicted_probs] = ...
            fernsClfApply(obj.windowdata.X,obj.classifier);
          obj.windowdata.isvalidprediction(:) = true;
          obj.ClearStatus();
      end
            
      % transfer to predictidx for current fly
      obj.UpdatePredictedIdx();
      
    end

    function UpdatePredictedIdx(obj)
      n = obj.t1_curr - obj.t0_curr + 1;
      obj.predictedidx = zeros(1,n);
      if isempty(obj.windowdata.exp),
        return;
      end
      idxcurr = obj.windowdata.exp == obj.expi & ...
        all(bsxfun(@eq,obj.windowdata.flies,obj.flies),2) & ...
        obj.windowdata.isvalidprediction;
      obj.predictedidx(obj.windowdata.t(idxcurr)-obj.t0_curr+1) = ...
        obj.windowdata.predicted(idxcurr);

      obj.UpdateErrorIdx();
            
    end
    
    function UpdateErrorIdx(obj)
      n = obj.t1_curr - obj.t0_curr + 1;
      obj.erroridx = zeros(1,n);
      obj.suggestedidx = zeros(1,n);
      idxcurr = obj.predictedidx ~= 0 & obj.labelidx ~= 0;
      obj.erroridx(idxcurr) = double(obj.predictedidx(idxcurr) ~= obj.labelidx(idxcurr))+1;
      
      idxcurr = obj.predictedidx ~= 0 & obj.labelidx == 0;
      obj.suggestedidx(idxcurr) = obj.predictedidx(idxcurr);
    end
        
    function Predict(obj,expi,flies,ts)
      
      % TODO: don't store window data just because predicting. 
      
      if isempty(obj.classifier),
        return;
      end

      if isempty(ts),
        return;
      end
            
      % compute window data
      [success,msg] = obj.PreLoadWindowData(expi,flies,ts);
      if ~success,
        warning(msg);
        return;
      end
      
      % indices into windowdata
      idxcurr = obj.windowdata.exp == expi & all(bsxfun(@eq,obj.windowdata.flies,flies),2) & ...
        ~obj.windowdata.isvalidprediction & ismember(obj.windowdata.t,ts);
      
      % apply classifier
      switch obj.classifiertype,
        
        case 'ferns',
          obj.SetStatus('Applying fern classifier to %d windows',nnz(idxcurr));
          [obj.windowdata.predicted(idxcurr),...
            obj.windowdata.predicted_probs(idxcurr,:)] = ...
            fernsClfApply(obj.windowdata.X(idxcurr,:),obj.classifier);
          obj.windowdata.isvalidprediction(idxcurr) = true;
        obj.ClearStatus();
      end
           
      obj.UpdatePredictedIdx();
      
    end

    function SetStatus(obj,varargin)

      if isempty(obj.setstatusfn),
        fprintf(varargin{:});
        fprintf('\n');
      else
        obj.setstatusfn(sprintf(varargin{:}));
        drawnow;
      end
      
    end
    
    function ClearStatus(obj)
      
      if ~isempty(obj.clearstatusfn),
        obj.clearstatusfn();
        drawnow;
      end
    
    end
  end
    
end