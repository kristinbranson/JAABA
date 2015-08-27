classdef JExperiments < handles
  
  properties (Access = public)
    expname = '';
    expdir = '';
    outexpdir = '';
    
    fileexists = [];
    filetimestamps = [];
    
    windowdata = [];
    scores = [];
    nflies = [];
    rootoutdir = '';
    clipsdir = '';
    
    labels = struct('t0s',{},'t1s',{},'names',{},'flies',{},'off',{},'timestamp',{},'imp_t0s',{},'imp_t1s',{});
    gt_labels = struct('t0s',{},'t1s',{},'names',{},'flies',{},'off',{},'timestamp',{},'imp_t0s',{},'imp_t1s',{});
    
    % cell array of arrays of first frame of each trajectory for each
    % experiment: firstframes_per_exp{expi}(fly) is the first frame of the
    % trajectory of fly for experiment expi.
    firstframes = [];
    % cell array of arrays of end frame of each trajectory for each
    % experiment: endframes_per_exp{expi}(fly) is the last frame of the
    % trajectory of fly for experiment expi.
    endframes = [];
    
    % sex per experiment, fly
    frac_sex=[];
    sex=[];
    hassex=[];
    hasperframesex=[];
    trxInfoExists = false;
    
    randomGTSuggestions = [];
    loadedGTSuggestions = [];
  end
  
  methods (Access = public)
    
    %%%% FIXED %%%%
    function [success, msg] = SetExpdir(expdir,projconf,rootoutdir)
      
      success = true; msg = '';
      
      if nargin<3, rootoutdir = ''; end
      obj.rootoutdir = rootoutdir;
      
      if ~exist(expdir,'file'),
        success = false;
        msg = sprintf('expdir %s does not exist',expdir);
        return;
      end
      
      obj.expdir = expdir;
      [pname,fname] = myfileparts(expdir);
      obj.expname = fname;
      
      if ~isempty(obj.rootoutdir)
        obj.outexpdir = fullfile(obj.rootoutdir,obj.expname);
        if ~exist(obj.outexpdir,'dir'),
          [success,msg1] = mkdir(obj.rootoutputdir,obj.expname);
          if ~success,
            msg = (sprintf('Could not create output directory %s, failed to set expdirs: %s',outexpdir,msg1));
            return;
          end
          
        end
      else
        obj.outexpdir = obj.expdir;
      end

      outclipsdir = fullfile(obj.outexpdir,projconf.GetFileName('clips',obj));
      if ~exist(outclipsdir,'dir'),
        [success,msg1] = mkdir(obj.outexpdir,clipsdir);
        if ~success, obj.clipsdir = ''; end
      end
      
      obj.windowdata = JWindowdata();
      [success,msg] = obj.LoadLabelsFromFile(projconf);
      if ~success, return; end;
      [success,msg] = obj.LoadGTLabelsFromFile(projconf);
      if ~success, return; end;
      
    end
    
    %%%% FIXED %%%%
    function name = GetName(obj)
      name = obj.expname;
    end
    
    %%%% FIXED %%%%
    function nflies = GetNumFlies(obj)
      nflies = obj.numflies;
    end
    
    %%%% FIXED %%%%
    function outexpdir = GetOutexpdir(obj)
      outexpdir = obj.outexpdir;
    end
    
    %%%% FIXED %%%%
    function expdir = GetExpdir(obj)
      expdir = obj.expdir;
    end
    
    %%%% FIXED %%%%
    function SetFilestatus(obj,fileis,fileexists,filetimestamps)
      obj.fileexists(fileis) = fileexists;
      obj.filetimestamps(fileis) = filetimestamps;
    end
    
    %%%% FIXED %%%%
    function [fileexists,filetimestamps] = GetFilestatus(obj)
      fileexists = obj.fileexists;
      filetimestamps = obj.filetimestamps;
    end      
    
    %%%% FIXED %%%%
    function [success,msg] = LoadLabelsFromFile(obj,projconf)
      % [success,msg] = LoadLabelsFromFile(obj,expi)
      % If the label file exists, this function loads labels for experiment
      % expi into obj.labels. Otherwise, it sets the labels to be empty. This
      % does not currently update the windowdata and labelidx (TODO).
      
      success = false; msg = '';
      labelfilename = projconf.GetFile('label',obj);
      
      if exist(labelfilename,'file'),
        
        obj.SetStatus('Loading labels for %s',obj.expname);
        
        try
          loadedlabels = load(labelfilename,'t0s','t1s','names','flies','off','timestamp');
          obj.labels.t0s = loadedlabels.t0s;
          obj.labels.t1s = loadedlabels.t1s;
          obj.labels.names = loadedlabels.names;
          obj.labels.flies = loadedlabels.flies;
          obj.labels.off = loadedlabels.off;
          obj.labelstats.nflies_labeled = size(loadedlabels.flies,1);
          obj.labelstats.nbouts_labeled = numel([loadedlabels.t0s{:}]);

          % Timestamps
          if iscell(loadedlabels.timestamp)
            obj.labels.timestamp = loadedlabels.timestamp;
          else
            for ndx = 1:numel(loadedlabels.flies)
              nBouts = numel(loadedlabels.t0s{ndx});
              if isempty(loadedlabels.timestamp)
                obj.labels.timestamp{ndx}(1:nBouts) = now;
              else
                obj.labels.timestamp{ndx}(1:nBouts) = loadedlabels.timestamp;
              end
            end
          end

          if ~isempty(whos('-file',labelfilename,'imp_t0s'))
            loadedimp = load(labelfilename,'imp_t0s','imp_t1s');
            obj.labels.imp_t0s = loadedimp.imp_t0s;
            obj.labels.imp_t1s = loadedimp.imp_t1s;
          else
            obj.labels.imp_t0s = cell(1,numel(loadedlabels.flies));
            obj.labels.imp_t1s = cell(1,numel(loadedlabels.flies));
          end
        catch ME,
          msg = getReport(ME);
          obj.ClearStatus();
          return;
        end
        
        obj.ClearStatus();
        
      else
        
        obj.labels.t0s = {};
        obj.labels.t1s = {};
        obj.labels.names = {};
        obj.labels.flies = [];
        obj.labels.off = [];
        obj.labels.timestamp = {};
        obj.labels.imp_t0s = {};
        obj.labels.imp_t1s = {};
        obj.labelstats.nflies_labeled = 0;
        obj.labelstats.nbouts_labeled = 0;
        
      end
      
      success = true;
      
    end
    
    %%%% FIXED %%%%
    function [success,msg] = LoadGTLabelsFromFile(obj,projconf)
      % [success,msg] = LoadGTLabelsFromFile(obj,expi)
      % If the label file exists, this function loads labels for experiment
      % expi into obj.gt_labels. Otherwise, it sets the gt_labels to be empty. This
      % does not currently update the windowdata and labelidx (TODO).
      
      success = false; msg = '';
      
      labelfilename = projconf.GetFile('gt_label',obj);
      
      if exist(labelfilename,'file'),
        
        obj.SetStatus('Loading labels for %s',obj.expdirs{expi});
        
        try
          loadedlabels = load(labelfilename,'t0s','t1s','names','flies','off','timestamp');
          obj.gt_labels.t0s = loadedlabels.t0s;
          obj.gt_labels.t1s = loadedlabels.t1s;
          obj.gt_labels.names = loadedlabels.names;
          obj.gt_labels.flies = loadedlabels.flies;
          obj.gt_labels.off = loadedlabels.off;
          obj.gt_labelstats.nflies_labeled = size(loadedlabels.flies,1);
          obj.gt_labelstats.nbouts_labeled = numel([loadedlabels.t0s{:}]);
          
          if iscell(loadedlabels.timestamp)
            obj.gt_labels.timestamp = loadedlabels.timestamp;
          else
            for ndx = 1:numel(loadedlabels.flies)
              nBouts = numel(loadedlabels.t0s{ndx});
              if isempty(loadedlabels.timestamp)
                obj.gt_labels.timestamp{ndx}(1:nBouts) = now;
              else
                obj.gt_labels.timestamp{ndx}(1:nBouts) = loadedlabels.timestamp;
              end
            end
          end
          
          if ~isempty(whos('-file',labelfilename,'imp_t0s'))
            loadedimp = load(labelfilename,'imp_t0s','imp_t1s');
            obj.gt_labels.imp_t0s = loadedimp.imp_t0s;
            obj.gt_labels.imp_t1s = loadedimp.imp_t1s;
          else
            obj.gt_labels.imp_t0s = cell(1,numel(loadedlabels.flies));
            obj.gt_labels.imp_t1s = cell(1,numel(loadedlabels.flies));
          end
          
        catch ME,
          msg = getReport(ME);
          obj.ClearStatus();
          return;
        end
        
        obj.ClearStatus();
        
      else
        
        obj.gt_labels.t0s = {};
        obj.gt_labels.t1s = {};
        obj.gt_labels.names = {};
        obj.gt_labels.flies = [];
        obj.gt_labels.off = [];
        obj.gt_labels.timestamp = {};
        obj.gt_labels.imp_t0s = {};
        obj.gt_labels.imp_t1s = {};
        obj.gt_labelstats.nflies_labeled = 0;
        obj.gt_labelstats.nbouts_labeled = 0;
      end
      
      success = true;
      
    end
    
    %%%% FIXED %%%%
    function [labelidx,T0,T1] = GetLabelIdx(obj,flies,T0,T1)
      % [labelidx,T0,T1] = GetLabelIdx(obj,expi,flies)
      % Returns the labelidx for the input experiment and flies read from
      % labels.
      
      if nargin < 4,
        T0 = max(obj.GetTrxFirstFrame(flies));
        T1 = min(obj.GetTrxEndFrame(flies));
      end
      n = T1-T0+1;
      off = 1 - T0;
      labels_curr = obj.GetLabels(flies);
      labelidx.vals = zeros(1,n);
      labelidx.imp = zeros(1,n);
      labelidx.timestamp = zeros(1,n);
      
      labelnames = obj.data.projconf.GetLabelNames();
      
      for i = 1:numel(labelnames),
        labelname = labelnames{i};
        for j = find(strcmp(labels_curr.names,labelname)),
          t0 = labels_curr.t0s(j);
          t1 = labels_curr.t1s(j);
          if t0>T1 || t1<T0; continue;end
          t0 = max(T0,t0);
          t1 = min(T1,t1);
          labelidx.vals(t0+off:t1-1+off) = i;
          labelidx.timestamp(t0+off:t1-1+off) = labels_curr.timestamp(j);
        end
      end
      for j = 1:numel(labels_curr.imp_t0s)
        t0 = labels_curr.imp_t0s(j); t1 = labels_curr.imp_t1s(j);
        labelidx.imp(t0+off:t1-1+off) = 1;
      end
      
    end
    
    %%%% FIXED %%%%
    function perframefiles = GetPerFrameFiles(obj)
      perframefiles = obj.data.projconf.GetPerFrameFiles(obj);
    end
    
    function [perframedata,T0,T1] = GetPerFrameData(obj,expi,flies,prop,T0,T1)
      % [perframedata,T0,T1] = GetPerFrameData(obj,expi,flies,prop,T0,T1)
      % Returns the per-frame data for the input experiment, flies, and
      % property.
      
      if ischar(prop),
        prop = find(strcmp(prop,handles.perframefn),1);
        if isempty(prop),
          error('Property %s is not a per-frame property');
        end
      end
      
      if obj.IsCurFly(expi,flies)
        if nargin < 5,
          perframedata = obj.perframedata{prop};
          T0 = obj.t0_curr;
          T1 = obj.t0_curr + numel(perframedata) - 1;
        else
          T0 = max(T0,obj.t0_curr);
          T1 = min(T1,obj.t0_curr+numel(obj.perframedata{prop})-1);
          i0 = T0 - obj.t0_curr + 1;
          i1 = T1 - obj.t0_curr + 1;
          perframedata = obj.perframedata{prop}(i0:i1);
        end
        return;
      end
      
      perframedir = obj.GetFile('perframedir',expi);
      tmp = load(fullfile(perframedir,[obj.allperframefns{prop},'.mat']));
      if nargin < 5,
        T0 = max(obj.GetTrxFirstFrame(expi,flies));
        % TODO: generalize to multi-fly
        perframedata = tmp.data{flies(1)};
        T1 = T0 + numel(perframedata) - 1;
        return;
      end
      off = 1 - GetTrxFirstFrame(expi,flies);
      i0 = T0 + off;
      i1 = T1 + off;
      perframedata = tmp.data{flies(1)}(i0:i1);
      
    end
    
    function perframedata = GetPerFrameData1(obj,expi,flies,prop,t)
      % perframedata = GetPerFrameData1(obj,expi,flies,prop,t)
      % Returns the per-frame data for the input experiment, flies, and
      % property.
      
      %       if ischar(prop),
      %         prop = find(strcmp(prop,handles.perframefn),1);
      %         if isempty(prop),
      %           error('Property %s is not a per-frame property');
      %         end
      %       end
      
      if ~isempty(obj.expi) && expi == obj.expi && numel(flies) == numel(obj.flies) && all(flies == obj.flies),
        is = t-obj.t0_curr+1;
        badidx = is > numel(obj.perframedata{prop});
        if any(badidx),
          perframedata = nan(size(is));
          perframedata(~badidx) = obj.perframedata{prop}(is(~badidx));
        else
          perframedata = obj.perframedata{prop}(is);
        end
        return;
      end
      
      perframedir = obj.GetFile('perframedir',expi);
      tmp = load(fullfile(perframedir,[obj.allperframefns{prop},'.mat']));
      off = 1 - obj.GetTrxFirstFrame(expi,flies);
      perframedata = tmp.data{flies(1)}(t+off);
      
    end
    
    function [prediction,T0,T1] = GetPredictedIdx(obj,expi,flies,T0,T1)
      
      if ~isempty(obj.expi) && numel(flies) == numel(obj.flies) && obj.IsCurFly(expi,flies),
        if nargin < 4,
          prediction = struct('predictedidx',obj.predictedidx,...
            'scoresidx', obj.scoresidx);
          T0 = obj.t0_curr;
          T1 = obj.t1_curr;
        else
          prediction = struct(...
            'predictedidx', obj.predictedidx(T0+obj.labelidx_off:T1+obj.labelidx_off),...
            'scoresidx',  obj.scoresidx(T0+obj.labelidx_off:T1+obj.labelidx_off));
        end
        return;
      end
      
      if nargin < 4,
        T0 = max(obj.GetTrxFirstFrame(expi,flies));
        T1 = min(obj.GetTrxEndFrame(expi,flies));
      end
      
      n = T1-T0+1;
      off = 1 - T0;
      prediction = struct('predictedidx', zeros(1,n),...
        'scoresidx', zeros(1,n));
      
      
      if ~isempty(obj.windowdata.exp)
        idxcurr = obj.FlyNdx(expi,flies) & ...
          obj.windowdata.t >= T0 & obj.windowdata.t <= T1 & ...
          obj.windowdata.isvalidprediction;
        prediction.predictedidx(obj.windowdata.t(idxcurr)+off) = ...
          obj.windowdata.predicted(idxcurr);
        prediction.scoresidx(obj.windowdata.t(idxcurr)+off) = ...
          obj.windowdata.scores(idxcurr);
      end
    end
    
    %%%% WHAT? %%%%
    function [idx,T0,T1] = IsBehavior(obj,behaviori,expi,flies,T0,T1)
      % [idx,T0,T1] = IsBehavior(obj,behaviori,expi,flies,T0,T1)
      % Returns whether the behavior is labeled as behaviori for experiment
      % expi, flies from frames T0 to T1. If T0 and T1 are not input, then
      % firstframe to endframe are used.
      
      if ~isempty(obj.expi) && expi == obj.expi && numel(flies) == numel(obj.flies) && all(flies == obj.flies),
        if nargin < 4,
          idx = obj.labelidx.vals == behaviori;
          T0 = obj.t0_curr;
          T1 = obj.t1_curr;
        else
          idx = obj.labelidx.vals(T0+obj.labelidx_off:T1+obj.labelidx_off) == behaviori;
        end
        return;
      end
      
      if nargin < 4,
        T0 = max(obj.GetTrxFirstFrame(expi,flies));
        T1 = min(obj.GetTrxEndFrame(expi,flies));
      end
      n = T1-T0+1;
      off = 1 - T0;
      labels_curr = obj.GetLabels(expi,flies);
      idx = false(1,n);
      for j = find(strcmp(labels_curr.names,obj.labelnames{behaviori})),
        t0 = labels_curr.t0s(j);
        t1 = labels_curr.t1s(j);
        idx(t0+off:t1-1+off) = true;
      end
      
    end
    
    %%%% FIXED %%%%
    function labels_curr = GetLabels(obj,flies)
      % labels_curr = GetLabels(obj,expi,flies)
      % Returns the labels for the input
      
      labels_curr = struct('t0s',[],'t1s',[],'names',{{}},'timestamp',[],'off',0,'imp_t0s',[],'imp_t1s',[]);
      
      if obj.IsGTMode()
        labelsToUse = obj.gt_labels;
      else
        labelsToUse = obj.labels;
      end
      
      [ism,fliesi] = ismember(flies,labelsToUse.flies,'rows');
      if ism,
        labels_curr.t0s = labelsToUse.t0s{fliesi};
        labels_curr.t1s = labelsToUse.t1s{fliesi};
        labels_curr.names = labelsToUse.names{fliesi};
        labels_curr.off = labelsToUse.off(fliesi);
        labels_curr.imp_t0s = labelsToUse.imp_t0s{fliesi};
        labels_curr.imp_t1s = labelsToUse.imp_t1s{fliesi};
        labels_curr.timestamp = labelsToUse.timestamp{fliesi};
      else
        t0_curr = max(obj.GetTrxFirstFrame(expi,flies));
        labels_curr.off = 1-t0_curr;
      end
      
    end
    
    %%%% FIXED %%%%
    function StoreLabels(obj,cache)
      % Store labels cached in labelidx for the current experiment and flies
      % to labels structure. This is when the timestamp on labels gets
      % updated.
      
      % flies not yet initialized
      cachedTarget = cache.getTarget();
      [cachedLabels,cachedLabels_off] = cache.getLabels();
      
      if isempty(cachedTarget) || isempty(cachedLabels),
        return;
      end
      
      obj.StoreLabels1(cachedTarget,cachedLabels,cachedLabels_off);
      
%       % preload labeled window data while we have the per-frame data loaded
%       ts = find(cachedLabels.vals~=0) - cachedLabels_off;
%       [success,msg] = obj.PreLoadWindowData(cache,cachedTarget,ts);
%       if ~success,
%         warning(msg);
%       end
      
    end
    
    %%%% FIXED %%%%
    function StoreLabels1(obj,target,labels,labels_off)
      
      % update labels
      newlabels = struct('t0s',[],'t1s',[],'names',{{}},'flies',[],'timestamp',[],'imp_t0s',[],'imp_t1s',[]);
      labelnames = obj.data.projconf.GetLabelNames();
      nbehaviors = numel(labelnames);
      for j = 1:nbehaviors,
        [i0s,i1s] = get_interval_ends(labels.vals==j);
        
        if ~isempty(i0s),
          n = numel(i0s);
          newlabels.t0s(end+1:end+n) = i0s - labels_off;
          newlabels.t1s(end+1:end+n) = i1s - labels_off;
          newlabels.names(end+1:end+n) = repmat(labelnames(j),[1,n]);
          newlabels.timestamp(end+1:end+n) = labels.timestamp(i0s);
        end
      end
      
      [i0s,i1s] = get_interval_ends(labelidx.imp);
      if ~isempty(i0s),
        newlabels.imp_t0s = i0s - labelidx_off;
        newlabels.imp_t1s = i1s - labelidx_off;
      end
      
      % Store labels according to the mode
      if obj.data.IsGTMode(),
        labelsToUse = 'gt_labels';
      else
        labelsToUse = 'labels';
      end
      
      [ism,j] = ismember(target,obj.(labelsToUse).target,'rows');
      if ~ism,
        j = size(obj.(labelsToUse).target,1)+1;
      end
      
      obj.(labelsToUse).t0s{j} = newlabels.t0s;
      obj.(labelsToUse).t1s{j} = newlabels.t1s;
      obj.(labelsToUse).names{j} = newlabels.names;
      obj.(labelsToUse).target(j,:) = target;
      obj.(labelsToUse).off(j) = labelidx_off;
      obj.(labelsToUse).timestamp{j} = newlabels.timestamp;
      obj.(labelsToUse).imp_t0s{j} = newlabels.imp_t0s;
      obj.(labelsToUse).imp_t1s{j} = newlabels.imp_t1s;
      
    end
    
    %%%% WHAT? %%%%
    function isstart = IsLabelStart(obj,expi,flies,ts)
      
      if obj.expi == expi && all(flies == obj.flies),
        isstart = obj.labelidx.vals(ts+obj.labelidx_off) ~= 0 & ...
          obj.labelidx.vals(ts+obj.labelidx_off-1) ~= obj.labelidx.vals(ts+obj.labelidx_off);
      else
        
        if obj.IsGTMode(),
          labelsToUse = 'gt_labels';
        else
          labelsToUse = 'labels';
        end
        
        [ism,fliesi] = ismember(flies,obj.(labelsToUse).flies,'rows');
        if ism,
          isstart = ismember(ts,obj.(labelsToUse).t0s{fliesi});
        else
          isstart = false(size(ts));
        end
      end
      
    end
    
    function ClearLabels(obj,expi,flies)
      
      if obj.nexps == 0,
        return;
      end
      
      if obj.IsGTMode()
        labelsToUse = 'gt_labels';
        labelstatsToUse = 'gt_labelstats';
      else
        labelsToUse = 'labels';
        labelstatsToUse = 'labelstats';
      end
      
      timestamp = now;
      
      % use all experiments by default
      if nargin < 2,
        expi = 1:obj.nexps;
      end
      
      % delete all flies by default
      if nargin < 3,
        for i = expi(:)',
          obj.(labelsToUse).t0s = {};
          obj.(labelsToUse).t1s = {};
          obj.(labelsToUse).names = {};
          obj.(labelsToUse).flies = [];
          obj.(labelsToUse).off = [];
          obj.(labelsToUse).timestamp = {};
          obj.(labelstatsToUse).nflies_labeled = 0;
          obj.(labelstatsToUse).nbouts_labeled = 0;
          obj.(labelsToUse).imp_t0s = {};
          obj.(labelsToUse).imp_t1s = {};
        end
      else
        if numel(expi) > 1,
          error('If flies input to ClearLabels, expi must be a single experiment');
        end
        % no labels
        if numel(obj.(labelsToUse)) < expi,
          return;
        end
        % which index of labels
        [~,flyis] = ismember(obj.(labelsToUse).flies,flies,'rows');
        for flyi = flyis(:)',
          % keep track of number of bouts so that we can update stats
          ncurr = numel(obj.(labelsToUse).t0s{flyi});
          obj.(labelsToUse).t0s{flyi} = [];
          obj.(labelsToUse).t1s{flyi} = [];
          obj.(labelsToUse).names{flyi} = {};
          obj.(labelsToUse).timestamp{flyi} = [];
          obj.(labelsToUse).imp_t0s{flyi} = [];
          obj.(labelsToUse).imp_t1s{flyi} = [];
          % update stats
          obj.(labelstatsToUse).nflies_labeled = obj.(labelstatsToUse).nflies_labeled - 1;
          obj.(labelstatsToUse).nbouts_labeled = obj.(labelstatsToUse).nbouts_labeled - ncurr;
        end
      end
      
      % clear labelidx if nec
      if ismember(obj.expi,expi) && ((nargin < 3) || ismember(obj.flies,flies,'rows')),
        obj.labelidx.vals(:) = 0;
        obj.labelidx.imp(:) = 0;
        obj.labelidx.timestamp(:) = 0;
      end
      
      % clear windowdata labelidx_new
      for i = expi(:)',
        if nargin < 3,
          idx = obj.windowdata.exp == i;
        else
          idx = obj.windowdata.exp == i & ismember(obj.windowdata.flies,flies,'rows');
        end
        obj.windowdata.labelidx_new(idx) = 0;
        obj.windowdata.labelidx_imp(idx) = 0;
        obj.UpdateErrorIdx();
      end
      
    end
    
    function [success,msg] = PredictWholeMovie(obj,classifierH)
      
      success = true; msg = '';
      
      if ~classifierH.istrained(),
        success = false;
        msg = 'No Classifier have been trained yet';
        return;
      end
      
      numFlies = obj.numflies;
      scoresA = cell(1,numFlies);
      
      tStartAll = obj.GetTrxFirstFrame(expi);
      tEndAll = obj.GetTrxEndFrame(expi);
      perframefile = obj.GetPerframeFiles(expi);
      windowfeaturescellparams = obj.windowfeaturescellparams;
      curperframefns = obj.curperframefns;
      allperframefns = obj.allperframefns;
      classifier = obj.classifier;
      
      obj.SetStatus('Classifying current movie..');
      
      parfor flies = 1:numFlies
        blockSize = 5000;
        tStart = tStartAll(flies);
        tEnd = tEndAll(flies);
        
        scores = nan(1,tEnd);
        
        for curt0 = tStart:blockSize:tEnd
          curt1 = min(curt0+blockSize-1,tEnd);
          X = JLabelData.ComputeWindowDataChunkStatic(curperframefns,...
            allperframefns,perframefile,flies,windowfeaturescellparams,curt0-tStart+1,curt1-tStart+1);
          
          scores(curt0:curt1) = myBoostClassify(X,classifier);
        end
        scoresA{flies} = scores;
        fprintf('Prediction done for flynum:%d, total number of flies:%d\n',flies,numFlies);
      end
      
      allScores = struct;
      allScores.scores = scoresA;
      allScores.tStart = tStartAll;
      allScores.tEnd = tEndAll;
      for flies = 1:numFlies
        [i0s i1s] = get_interval_ends(allScores.scores{flies}>0);
        allScores.t0s{flies} = i0s;
        allScores.t1s{flies} = i1s;
      end
      obj.SaveScores(allScores,expi);
      obj.LoadScores(expi,obj.GetFile('scores',expi));
      obj.ClearStatus();
      
    end
    
    %%%% FIXED %%%%
    function flyStats = GetFlyStats(obj,expi,flyNum)
      % Calculates statistics such as number of labeled bouts, predicted bouts
      % and change in scores.
      
      [ism,j] = ismember(flyNum,obj.labels.flies,'rows');
      labelnames = obj.data.projconf.GetLabelNames();
      if ism,
        flyStats.nbouts = numel(obj.labels.t0s{j});
        posframes = 0; negframes = 0;
        for ndx = 1:numel(obj.labels.t0s{j})
          numFrames = obj.labels.t1s{j}(ndx)-obj.labels.t0s{j}(ndx);
          if strcmp(obj.labels.names{j}{ndx},labelnames{1})
            posframes = posframes + numFrames;
          else
            negframes = negframes + numFrames;
          end
        end
        flyStats.posframes = posframes;
        flyStats.negframes = negframes;
        flyStats.totalframes = posframes + negframes;
      else
        flyStats.nbouts = 0;
        flyStats.posframes = 0;
        flyStats.negframes = 0;
        flyStats.totalframes = 0;
      end
      
      [ism,j] = ismember(flyNum,obj.gt_labels.flies,'rows');
      if ism,
        flyStats.gt_nbouts = numel(obj.gt_labels.t0s{j});
        posframes = 0; negframes = 0;
        for ndx = 1:numel(obj.gt_labels.t0s{j})
          numFrames = obj.gt_labels.t1s{j}(ndx)-obj.gt_labels.t0s{j}(ndx);
          if strcmp(obj.gt_labels.names{j}{ndx},labelnames{1})
            posframes = posframes + numFrames;
          else
            negframes = negframes + numFrames;
          end
        end
        flyStats.gt_posframes = posframes;
        flyStats.gt_negframes = negframes;
        flyStats.gt_totalframes = posframes + negframes;
      else
        flyStats.gt_nbouts = 0;
        flyStats.gt_posframes = 0;
        flyStats.gt_negframes = 0;
        flyStats.gt_totalframes = 0;
      end
      
      flyStats.endframe = obj.GetTrxEndFrame(flyNum);
      flyStats.firstframe = obj.GetTrxFirstframe(flyNum);
      flyStats.trajLength = flyStats.endframe-flyStats.firstframe+1;
      
      if obj.hassex,
        if obj.hasperframesex,
          sexfrac = obj.GetSexFrac(flyNum);
          flyStats.sexfrac = round(100*sexfrac.M);
        else
          flyStats.sexfrac = 100*strcmpi(obj.GetSex(expi,flyNum),'M');
        end
      else
        flyStats.sexfrac = [];
      end
      
      loadedscores = obj.scores.GetScores(flyNum);
      if ~isempty(curscores)
        flyStats.nscorepos_loaded = nnz(loadedscores>0);
        flyStats.nscoreneg_loaded = nnz(loadedscores<0);
      else
        flyStats.nscorepos_loaded = [];
        flyStats.nscoreneg_loaded = [];
      end
      
      curScores = obj.windowdata.GetScores(flyNum);
      curLabels = obj.windowdata.GetLabels(flyNum);
      
      if ~isemtpy(curscores)
        curPosMistakes = nnz( curScores<0 & curLabels ==1 );
        curNegMistakes = nnz( curScores>0 & curLabels >1 );
        
        flyStats.nscoreframes = numel(curScores);
        flyStats.nscorepos = nnz(curScores>0);
        flyStats.nscoreneg = nnz(curScores<0);
        flyStats.errorsPos = curPosMistakes;
        flyStats.errorsNeg = curNegMistakes;
      else
        flyStats.nscoreframes = [];
        flyStats.nscorepos = [];
        flyStats.nscoreneg = [];
        flyStats.errorsPos = [];
        flyStats.errorsNeg = [];
      end
      
      flyStats.one2two = [];
      flyStats.two2one = [];
      if ~isempty(obj.classifier_old) && nnz(curScores),
        oldScores = obj.windowdata.GetScoresOld(flyNum);
        flyStats.one2two = nnz(curScores < 0 & oldScores > 0);
        flyStats.two2one = nnz(curScores > 0 & oldScores < 0);
      end
      
      flyStats.validatedErrorsPos = [];
      flyStats.validatedErrorsNeg = [];
      scores_validated = obj.windowdata.GetScoresValidated(flyNum);
      if ~isempty(scores_validated),
        curNdx = obj.FlyNdx(expi,flyNum);
        if nnz(curNdx);
          curScores = obj.windowdata.scores_validated(curNdx);
          curLabels = obj.windowdata.labelidx_cur(curNdx);
          
          curPosMistakes = nnz( curScores(:)<0 & curLabels(:) ==1 );
          curNegMistakes = nnz( curScores(:)>0 & curLabels(:) >1 );
          
          flyStats.validatedErrorsPos = curPosMistakes;
          flyStats.validatedErrorsNeg = curNegMistakes;
        end
      end
      
      flyStats.gt_suggestion_frames = nnz(obj.data.GetGTSuggestionIdx(expi,flyNum));
      
      %       if ~isempty(obj.windowdata.X)
      %         idxcurr = obj.windowdata.exp==expi & obj.windowdata.flies == flyNum;
      %         flyStats.npredictframes = nnz(idxcurr);
      %         flyStats.npredictfrac = nnz(obj.windowdata.scores(idxcurr)>0)/flyStats.nscoreframes;
      %
      %       else
      %         flyStats.npredictframes = [];
      %         flyStats.npredictfrac = [];
      %       end
      
    end
    
    %%%% FIXED %%%%
    function [success,msg] = LoadTrxInfo(obj,projconf)
      
      success = true; msg = '';
      
      if obj.trxInfoExists, return; end
      trxfile = fullfile(obj.expdir,projconf.GetFileName('trx'));
      
      if ~exist(trxfile,'file'),
        msg = sprintf('Trx file %s does not exist, cannot count flies',trxfile);
        success = false;
        return;
      end
      
      global CACHED_TRX; %#ok<TLEV>
      global CACHED_TRX_EXPNAME; %#ok<TLEV>
      if isempty(CACHED_TRX) || isempty(CACHED_TRX_EXPNAME) || ...
          ~strcmp(obj.expname,CACHED_TRX_EXPNAME),
        hwait = mywaitbar(0,sprintf('Loading trx to determine number of flies for %s',obj.expnames{expi}),'interpreter','none');
        trx = load_tracks(trxfile);
        if ishandle(hwait), delete(hwait); end
        CACHED_TRX = trx;
        CACHED_TRX_EXPNAME = obj.expnames{expi};
      else
        fprintf('DEBUG: Using CACHED_TRX. REMOVE THIS\n');
        trx = CACHED_TRX;
      end
      
      obj.StoreTrxInfo(trx);
      
    end
    
    %%%% FIXED %%%% 
    function StoreTrxInfo(obj,trx)
      
      if obj.trxInfoExists, return; end
      
      obj.nflies = numel(trx);
      obj.firstframes = [trx.firstframe];
      obj.endframes = [trx.endframe];
      
      obj.hassex = obj.hassex || isfield(trx,'sex');
      
      % store sex info
      tmp = repmat({nan},[1,numel(trx)]);
      obj.frac_sex = struct('M',tmp,'F',tmp);
      obj.sex = repmat({'?'},[1,numel(trx)]);
      
      if ~isfield(trx,'sex'), return; end;

      if numel(trx) > 1,
        obj.hasperframesex = iscell(trx(1).sex);
      end
      
      if obj.hasperframesex,
        for fly = 1:numel(trx),
          n = numel(trx(fly).sex);
          nmale = nnz(strcmpi(trx(fly).sex,'M'));
          nfemale = nnz(strcmpi(trx(fly).sex,'F'));
          obj.frac_sex(fly).M = nmale/n;
          obj.frac_sex(fly).F = nfemale/n;
          if nmale > nfemale,
            obj.sex{fly} = 'M';
          elseif nfemale > nmale,
            obj.sex{fly} = 'F';
          else
            obj.sex{fly} = '?';
          end
        end
      else
        for fly = 1:numel(trx),
          obj.sex{fly} = trx(fly).sex;
          if strcmpi(trx(fly).sex,'M'),
            obj.frac_sex(fly).M = 1;
            obj.frac_sex(fly).F = 0;
          elseif strcmpi(trx(fly).sex,'F'),
            obj.frac_sex(fly).M = 0;
            obj.frac_sex(fly).F = 1;
          end
        end
      end
      
      obj.trxInfoExists = true;
      
    end
    
    %%%% WHAT? %%%%
    function pos = GetTrxPos1(varargin)
      % [x,y,theta,a,b] = GetTrxPos1(obj,expi,fly,ts)
      % Returns the position for the input experiment, SINGLE fly, and
      % frames. If ts is not input, then all frames are returned.
      
      % moved to separate file so that function could be easily modified
      pos = JLabelData_GetTrxPos(varargin{:});
      
    end
    
    %%%% FIXED %%%%
    function sexfrac = GetSexFrac(obj,fly)
      % x = GetSexFrac(obj,expi,fly)
      % Returns a struct indicating the fraction of frames for which the sex
      % of the fly is M, F
      
      sexfrac = obj.frac_sex_per_exp(fly);
      
    end
    
    %%%% FIXED %%%%
    function t0 = GetTrxFirstFrame(obj,target)
      % t0 = GetTrxFirstFrame(obj,expi,flies)
      % Returns the firstframes for the input experiment and flies. If flies
      % is not input, then all flies are returned.
      
      if nargin < 2;
        t0 = obj.firstframes;
        return;
      end
      
      t0 = obj.firstframes(target);
      
    end
    
    %%%% FIXED %%%%
    function t1 = GetTrxEndFrame(obj,flies)
      % t1 = GetTrxEndFrame(obj,expi,flies)
      % Returns the endframes for the input experiment and flies. If flies
      % is not input, then all flies are returned.
      
      if nargin < 2,
        t1 = obj.endframes;
        return;
      end
      
      t1 = obj.endframes(flies);
      
    end
    
    function LoadScoresDefault(obj,expi)
      sfn = obj.GetFile('scores',expi);
      if ~exist(sfn,'file')
        warndlg(sprintf('No scores file %s at the default location',...
          scoreFileName));
      end
      obj.LoadScores(expi,sfn);
    end
    
    %%%% FIXED %%%%
    function [success,msg] = CheckFileStatus(obj,projconf)
      msg = ''; success = false;
      
      filetypes = projconf.filetypes;
      fileis = 1:numel(filetypes);
      
      % initialize fileexists table
      
      allfilesexist = true;
      filesfixable = true;
      
      obj.fileexists = false(1,fileis);
      obj.filetimestamps = nan(1,fileis) ;
      
      % loop through all file types
      for filei = fileis,
        file = filetypes{filei};
        % loop through experiments
        
        if strcmpi(file,'perframedir'),
          
          [fn,timestamps] = projconf.GetPerframeFiles(expH);
          if isempty(fn),
            obj.fileexists(filei) = false;
            obj.filetimestamps(filei) = -inf;
          else
            obj.fileexists(filei) = all(cellfun(@(s) exist(s,'file'),fn));
            obj.filetimestamps(filei) = max(timestamps);
          end
          
        else
          
          % check for existence of current file(s)
          [fn,obj.filetimestamps(filei)] = projconf.GetFile(file,obj);
          if iscell(fn),
            obj.fileexists(filei) = all(cellfun(@(s) exist(s,'file'),fn));
          else
            obj.fileexists(filei) = exist(fn,'file');
          end
          
        end
        
        % if file doesn't exist and is required, then not all files exist
        if ~obj.fileexists(filei),
          if Projconf.IsRequiredFile(file),
            allfilesexist = false;
            % if furthermore file can't be generated, then not fixable
            if ~Projconf.CanGenerateFile(file),
              filesfixable = false;
              msg1 = sprintf('%s missing and cannot be generated.',file);
              if isempty(msg),
                msg = msg1;
              else
                msg = sprintf('%s\n%s',msg,msg1);
              end
            end
          end
        end
        
      end
      success = allfilesexist || filesfixable;
      
    end
    
    %%%% FIXED %%%%
    function [success,msg] = UpdateWindowData(obj)
      % [success,msg] = PreLoadWindowData(obj,expi,flies,ts)
      % Compute and store the window data for experiment expi, flies flies,
      % and all frames ts.
      % This function finds all frames that currently do not have window data
      % cached. In a loop, it finds the first frame that is missing window
      % data, and computes window data for all frames in a chunk of size
      % 2*obj.windowdatachunk_radius + 1 after this frame using the function
      % ComputeWindowDataChunk. Then, it updates the frames that are missing
      % window data. It proceeds in this loop until there are not frames
      % in the input ts missing window data.
      
      success = false; msg = '';
      
      obj.windowdata.TrimWindowData();
      
      for flies = 1:obj.nflies
        if ~any(obj.labels.flies==flies); continue; end

        labelidxStruct = obj.GetLabelIdx(flies);
        [success,msg1] = obj.windowdata.AddTs(flies,find(labelidxStruct~=0));
      end
      
      if ~success, msg = msg1; return; end
      
      success = true;
      
    end
    
  end
  
end