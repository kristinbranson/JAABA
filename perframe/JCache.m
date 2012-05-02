classdef JCache < handles
  
  properties (Access = public)

    % first frame that all flies currently selected are tracked
    t0_curr = 0;
    % last frame that all flies currently selected are tracked
    t1_curr = 0;

    expi = [];
    target = [];
    scores = [];
    labels = [];
    labels_gt = [];
    scores_old = [];
    scores_validated = [];
    scores_loaded = [];
    scores_processed = [];
    perframedata = {};
    % last-used trajectories (one experiment, all flies)
    trx = {};
    
    % labels for the current experiment and flies, represented as an array
    % such that labelidx(t+labelidx_off) is the index of the behavior for
    % frame t of the movie. labelidx(i) == 0 corresponds to
    % unlabeled/unknown, otherwise labelidx(i) corresponds to behavior
    % labelnames{labelidx{i})
    labelidx = struct('val',[],'imp',[],'timestamp',[]);
    labelidx_off = 0;
    
    % predicted label for current experiment and flies, with the same type
    % of representation as labelidx
    predictedidx = [];
    scoresidx = [];
    scoresidx_old = [];
    scoreTS = [];

    % whether the predicted label matches the true label. 0 stands for
    % either not predicted or not labeled, 1 for matching, 2 for not
    % matching. this has the same representation as labelidx.
    erroridx = [];
    
    % TODO: remove this
    % predictedidx for unlabeled data, same representation as labelidx
    suggestedidx = [];

  end
  
  methods (Access = public)
    
    function SetCurTarget(obj,target)
      obj.target = target;
      %TODO Read data ..
    end
    
    function val = IsCurFly(obj,expi,flies)
      val = all(flies == obj.flies) && (expi==obj.expi);
    end

    function res = IsCached(obj,expi,target)
      res = obj.expi == expi & obj.target ==target;
    end
    
    function expi = obj.GetExp(obj)
      expi = obj.expi;
    end
    
    function target = obj.GetTarget(obj)
      target = obj.target;
    end
    
    function labels = GetLabels(obj)
      labels = obj.labels;
    end
    
    function scores = GetScores(obj,type)
      
      switch type
        case 'current'
          scores = obj.scores;
        case 'old'
          scores = obj.scores_old;
        case 'loaded'
        case 'validated'
        case 'processed'
          
      end
    end
    
    function [success,msg] = PreLoad(obj,flies)
    % [success,msg] = PreLoad(obj,expi,flies)
    % Preloads data associated with the input experiment and flies. If
    % neither the experiment nor flies are changing, then we do nothing. If
    % there is currently a preloaded experiment, then we store the labels
    % in labelidx into labels using StoreLabels. We then load from labels
    % into labelidx for the new experiment and flies. We load the per-frame
    % data for this experiment and flies. If this is a different
    % experiment, then we load in the trajectories for this experiment.  
      
      success = false;
      msg = '';
      
      if numel(expi) ~= 1,
        error('expi must be a scalar');
      end

      if numel(unique(flies)) ~= numel(flies),
        msg = 'flies must all be unique';
        return;
      end
      
      diffexpi = isempty(obj.expi) || expi ~= obj.expi;
      diffflies = diffexpi || numel(flies) ~= numel(obj.flies) || ~all(flies == obj.flies);
      % nothing to do
      if ~diffflies,
        success = true;
        return;
      end

      if ~isempty(obj.expi) && obj.expi > 0,
        % store labels currently in labelidx to labels
        obj.StoreLabels();
      end
      
      if diffexpi,
        
        % load trx
%         try
          trxfilename = obj.GetFile('trx',expi);
          if ~exist(trxfilename,'file')
            msg = sprintf('Trx file %s does not exist',trxfilename);
            success = false;
            return;
          end
          
          obj.SetStatus('Loading trx for experiment %s',obj.expnames{expi});
                    
          % TODO: remove this
          global CACHED_TRX; %#ok<TLEV>
          global CACHED_TRX_EXPNAME; %#ok<TLEV>
          if isempty(CACHED_TRX) || isempty(CACHED_TRX_EXPNAME) || ...
              ~strcmp(obj.expnames{expi},CACHED_TRX_EXPNAME),
            obj.trx = load_tracks(trxfilename);
            CACHED_TRX = obj.trx;
            CACHED_TRX_EXPNAME = obj.expnames{expi};
          else
            fprintf('DEBUG: Using CACHED_TRX. REMOVE THIS\n');
            obj.trx = CACHED_TRX;
          end
          % store trx_info, in case this is the first time these trx have
          % been loaded
          [success,msg] = obj.GetTrxInfo(expi,true,obj.trx);
          if ~success,
            return;
          end
          
%         catch ME,
%           msg = sprintf('Error loading trx from file %s: %s',trxfilename,getReport(ME));
%           if ishandle(hwait),
%             delete(hwait);
%             drawnow;
%           end
%           return;
%         end
 
      end

      % set labelidx from labels
      obj.SetStatus('Caching labels for experiment %s, flies%s',obj.expnames{expi},sprintf(' %d',flies));
      [obj.labelidx,obj.t0_curr,obj.t1_curr] = obj.GetLabelIdx(expi,flies);
      obj.labelidx_off = 1 - obj.t0_curr;
      
      % load perframedata
      obj.SetStatus('Loading per-frame data for %s, flies %s',obj.expdirs{expi},mat2str(flies));
      file = obj.GetPerframeFiles(expi);
      for j = 1:numel(obj.allperframefns),
        if ~exist(file{j},'file'),
          msg = sprintf('Per-frame data file %s does not exist',file{j});
          return;
        end
%         try
          tmp = load(file{j});
          obj.perframedata{j} = tmp.data{flies(1)};
          obj.perframeunits{j} = tmp.units;
%         catch ME,
%           msg = getReport(ME);
%         end
      end
      
      obj.expi = expi;
      obj.flies = flies;

      obj.UpdatePredictedIdx();
      obj.ClearStatus();
           
      success = true;
      
    end
    
    function ClearCachedPerExpData(obj)
    % ClearCachedPerExpData(obj)
    % Clears all cached data for the currently loaded experiment
      obj.trx = {};
      obj.expi = 0;
      obj.flies = nan(size(obj.flies));
      obj.perframedata = {};
      obj.labelidx = struct('vals',[],'imp',[],'timestamp',[]);
      obj.labelidx_off = 0;
      obj.t0_curr = 0;
      obj.t1_curr = 0;
      obj.predictedidx = [];
      obj.scoresidx = [];
      obj.scoresidx_old = [];
      obj.erroridx = [];
      obj.suggestedidx = [];
    end

    function SetLabel(obj,expi,flies,ts,behaviori,important)
    % SetLabel(obj,expi,flies,ts,behaviori)
    % Set label for experiment expi, flies, and frames ts to behaviori. If
    % expi, flies match current expi, flies, then we only set labelidx.
    % Otherwise, we set labels. 
      
      if obj.IsCurFly(expi,flies),
        obj.labelidx.vals(ts+obj.labelidx_off) = behaviori;
        obj.labelidx.imp(ts+obj.labelidx_off) = important;
        obj.labelidx.timestamp(ts+obj.labelidx_off) = now;
      else
        [labelidx,T0] = obj.GetLabelIdx(expi,flies);
        labelidx.vals(ts+1-T0) = behaviori;
        labelidx.imp(ts+1-T0) = important;
        labelidx.timestamp(ts+1-T0) = now;
        obj.StoreLabels1(expi,flies,labelidx,1-T0);        
      end
      
    end
    
    function StoreLabels(obj)
    % Store labels cached in labelidx for the current experiment and flies
    % to labels structure. This is when the timestamp on labels gets
    % updated. 
      
      % flies not yet initialized
      if isempty(obj.flies) || all(isnan(obj.flies)) || isempty(obj.labelidx.vals),
        return;
      end
      
      obj.StoreLabels1(obj.expi,obj.flies,obj.labelidx,obj.labelidx_off);
            
      % preload labeled window data while we have the per-frame data loaded
      ts = find(obj.labelidx.vals~=0) - obj.labelidx_off;
      [success,msg] = obj.PreLoadWindowData(obj.expi,obj.flies,ts);
      if ~success,
        warning(msg);
      end

      % update windowdata's labelidx_new
      if ~isempty(obj.windowdata.exp),
        idxcurr = obj.windowdata.exp == obj.expi & ...
          all(bsxfun(@eq,obj.windowdata.flies,obj.flies),2);
        obj.windowdata.labelidx_new(idxcurr) = obj.labelidx.vals(obj.windowdata.t(idxcurr)+obj.labelidx_off);
        obj.windowdata.labelidx_imp(idxcurr) = obj.labelidx.imp(obj.windowdata.t(idxcurr)+obj.labelidx_off);
      end
      
      %obj.UpdateWindowDataLabeled(obj.expi,obj.flies);
      
    end

    function StoreLabels1(obj,expi,flies,labelidx,labelidx_off)
      
      % update labels
      newlabels = struct('t0s',[],'t1s',[],'names',{{}},'flies',[],'timestamp',[],'imp_t0s',[],'imp_t1s',[]);
      for j = 1:obj.nbehaviors,
        [i0s,i1s] = get_interval_ends(labelidx.vals==j);
        
        if ~isempty(i0s),
          n = numel(i0s);
          newlabels.t0s(end+1:end+n) = i0s - labelidx_off;
          newlabels.t1s(end+1:end+n) = i1s - labelidx_off;
          newlabels.names(end+1:end+n) = repmat(obj.labelnames(j),[1,n]);
          newlabels.timestamp(end+1:end+n) = labelidx.timestamp(i0s);
        end
      end
      
      [i0s,i1s] = get_interval_ends(labelidx.imp);
      if ~isempty(i0s),
        newlabels.imp_t0s = i0s - labelidx_off;
        newlabels.imp_t1s = i1s - labelidx_off;
      end
      
      % Store labels according to the mode
      if obj.IsGTMode(),
        labelsToUse = 'gt_labels';
        labelstatsToUse = 'gt_labelstats';
      else
        labelsToUse = 'labels';
        labelstatsToUse = 'labelstats';
      end
      
      [ism,j] = ismember(flies,obj.(labelsToUse)(expi).flies,'rows');
      if ~ism,
        j = size(obj.(labelsToUse)(expi).flies,1)+1;
      end

      obj.(labelsToUse)(expi).t0s{j} = newlabels.t0s;
      obj.(labelsToUse)(expi).t1s{j} = newlabels.t1s;
      obj.(labelsToUse)(expi).names{j} = newlabels.names;
      obj.(labelsToUse)(expi).flies(j,:) = flies;
      obj.(labelsToUse)(expi).off(j) = labelidx_off;
      obj.(labelsToUse)(expi).timestamp{j} = newlabels.timestamp;
      obj.(labelsToUse)(expi).imp_t0s{j} = newlabels.imp_t0s;
      obj.(labelsToUse)(expi).imp_t1s{j} = newlabels.imp_t1s;

      % store labelstats
      obj.(labelstatsToUse)(expi).nflies_labeled = numel(unique(obj.(labelsToUse)(expi).flies));
      obj.(labelstatsToUse)(expi).nbouts_labeled = numel(newlabels.t1s);
            
    end

  end
  
end