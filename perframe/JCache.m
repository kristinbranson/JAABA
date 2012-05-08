classdef JCache < handles
  
  properties (Access = public)
    
    % first frame that all flies currently selected are tracked
    t0_curr = 0;
    % last frame that all flies currently selected are tracked
    t1_curr = 0;
    
    expH = [];
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
    
    %%%% FIXED %%%%
    function res = IsCached(obj,expH,target)
      res = obj.expH == expH & obj.target ==target;
    end
    
    %%%% FIXED %%%%
    function expH = obj.GetExp(obj)
      expH = obj.expH;
    end
    
    %%%% FIXED %%%%
    function target = obj.GetTarget(obj)
      target = obj.target;
    end
    
    %%%% FIXED %%%%
    function [labels, labels_off] = GetLabels(obj)
      labels = obj.labels;
      labels_off = obj.labels_off;
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
    
    %%%% FIXED %%%%
    function [success,msg] = PreLoad(obj,expH,target)
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
      
      diffexpi = isempty(obj.expH) || expH ~= obj.expH;
      diffflies = diffexpi || numel(target) ~= numel(obj.target) || ~all(target == obj.target);
      % nothing to do
      if ~diffflies,
        success = true;
        return;
      end
      
      if ~isempty(obj.expH) ,
        % store labels currently in labelidx to labels
        expH.StoreLabels(obj);
      end
      
      if diffexpi,
        
        % load trx
        try
          trxfilename = obj.data.projconf.GetFile('trx',expH);
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
          expH.StoreTrxInfo(obj.trx);

        catch ME,
          msg = sprintf('Error loading trx from file %s: %s',trxfilename,getReport(ME));
          if ishandle(hwait),
            delete(hwait);
            drawnow;
          end
          return;
        end
        
      end
      
      % set labelidx from labels
      obj.SetStatus('Caching labels for experiment %s, target %s',expH.GetName,sprintf(' %d',target));
      [obj.labelidx,obj.t0_curr,obj.t1_curr] = expH.GetLabelIdx(target);
      obj.labelidx_off = 1 - obj.t0_curr;
      
      % load perframedata
      obj.SetStatus('Loading per-frame data for %s, target %s',expH.GetName,mat2str(target));
      file = obj.data.projconf.GetPerframeFiles(expH);
      for j = 1:numel(file),
        if ~exist(file{j},'file'),
          msg = sprintf('Per-frame data file %s does not exist',file{j});
          return;
        end
        try
          tmp = load(file{j});
          obj.perframedata{j} = tmp.data{flies(1)};
          obj.perframeunits{j} = tmp.units;
        catch ME,
          msg = getReport(ME);
          return;
        end
      end
      
      obj.expH = expH;
      obj.target = target;
      
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
    
    %%% FIXED %%%%
    function SetLabel(obj,expH,flies,ts,behaviori,important)
      % SetLabel(obj,expi,flies,ts,behaviori)
      % Set label for experiment expi, flies, and frames ts to behaviori. If
      % expi, flies match current expi, flies, then we only set labelidx.
      % Otherwise, we set labels.
      
      if obj.expH == expH && obj.flies == flies
        obj.labelidx.vals(ts+obj.labelidx_off) = behaviori;
        obj.labelidx.imp(ts+obj.labelidx_off) = important;
        obj.labelidx.timestamp(ts+obj.labelidx_off) = now;
      else
        [labelidx,T0] = expH.GetLabelIdx(flies);
        labelidx.vals(ts+1-T0) = behaviori;
        labelidx.imp(ts+1-T0) = important;
        labelidx.timestamp(ts+1-T0) = now;
        expH.StoreLabels1(flies,labelidx,1-T0);
      end
      
    end
    
    %%%% FIXED %%%%
    function StoreLabels(obj)
      expH.StoreLabels();
    end
    
    %%%% FIXED %%%%
    function out = GetTrxValues(obj,infoType,flies,ts)
      % A generic function that return track info.
      
      if nargin < 3,     % No flies given
        switch infoType
          case 'Trx'
            out = obj.trx;
          case 'X'
            out = {obj.trx.x};
          case 'Y'
            out = {obj.trx.y};
          case 'A'
            out = {obj.trx.a};
          case 'B'
            out = {obj.trx.b};
          case 'Theta'
            out = {obj.trx.theta};
          otherwise
            error('Incorrect infotype requested from GetTrxValues with less than 4 arguments');
        end
        return;
        
        
      elseif nargin < 4, % No ts given
        switch infoType
          case 'Trx'
            out = obj.trx(flies);
          case 'X'
            out = {obj.trx(flies).x};
          case 'Y'
            out = {obj.trx(flies).y};
          case 'A'
            out = {obj.trx(flies).a};
          case 'B'
            out = {obj.trx(flies).b};
          case 'Theta'
            out = {obj.trx(flies).theta};
          case 'X1'
            out = [obj.trx(flies).x];
          case 'Y1'
            out = [obj.trx(flies).y];
          case 'A1'
            out = [obj.trx(flies).a];
          case 'B1'
            out = [obj.trx(flies).b];
          case 'Theta1'
            out = [obj.trx(flies).theta];
          otherwise
            error('Incorrect infotype requested from GetTrxValues');
        end
        return
      else               % Everything is given
        nflies = numel(flies);
        fly = flies(1);
        switch infoType
          case 'Trx'
            c = cell(1,nflies);
            trx = struct('x',c,'y',c,'a',c,'b',c,'theta',c,'ts',c,'firstframe',c,'endframe',c);
            for i = 1:numel(flies),
              fly = flies(i);
              js = min(obj.trx(fly).nframes,max(1,ts + obj.trx(fly).off));
              trx(i).x = obj.trx(fly).x(js);
              trx(i).y = obj.trx(fly).y(js);
              trx(i).a = obj.trx(fly).a(js);
              trx(i).b = obj.trx(fly).b(js);
              trx(i).theta = obj.trx(fly).theta(js);
              trx(i).ts = js-obj.trx(fly).off;
              trx(i).firstframe = trx(i).ts(1);
              trx(i).endframe = trx(i).ts(end);
            end
            out = trx;
          case 'X'
            x = cell(1,nflies);
            for i = 1:numel(flies),
              fly = flies(i);
              js = min(obj.trx(fly).nframes,max(1,ts + obj.trx(fly).off));
              x{i} = obj.trx(fly).x(js);
            end
            out = x;
          case 'Y'
            x = cell(1,nflies);
            for i = 1:numel(flies),
              fly = flies(i);
              js = min(obj.trx(fly).nframes,max(1,ts + obj.trx(fly).off));
              x{i} = obj.trx(fly).y(js);
            end
            out = x;
          case 'A'
            x = cell(1,nflies);
            for i = 1:numel(flies),
              fly = flies(i);
              js = min(obj.trx(fly).nframes,max(1,ts + obj.trx(fly).off));
              x{i} = obj.trx(fly).a(js);
            end
            out = x;
          case 'B'
            x = cell(1,nflies);
            for i = 1:numel(flies),
              fly = flies(i);
              js = min(obj.trx(fly).nframes,max(1,ts + obj.trx(fly).off));
              x{i} = obj.trx(fly).b(js);
            end
            out = x;
          case 'Theta'
            x = cell(1,nflies);
            for i = 1:numel(flies),
              fly = flies(i);
              js = min(obj.trx(fly).nframes,max(1,ts + obj.trx(fly).off));
              x{i} = obj.trx(fly).theta(js);
            end
            out = x;
          case 'X1'
            out = obj.trx(fly).x(ts + obj.trx(fly).off);
          case 'Y1'
            out = obj.trx(fly).y(ts + obj.trx(fly).off);
          case 'A1'
            out = obj.trx(fly).a(ts + obj.trx(fly).off);
          case 'B1'
            out = obj.trx(fly).b(ts + obj.trx(fly).off);
          case 'Theta1'
            out = obj.trx(fly).theta(ts + obj.trx(fly).off);
          otherwise
            error('Incorrect infotype requested from GetTrxValues');
        end
      end
      
    end
    
    %%%% FIXED %%%%
    function sex = GetSex(obj,fly,ts,fast)
      % x = GetSex(obj,expi,fly,ts)
      % Returns the sex for the input experiment, SINGLE fly, and
      % frames. If ts is not input, then all frames are returned.
      
      if ~obj.expH.hassex,
        sex = '?';
        return;
      end
      
      if nargin < 4,
        fast = false;
      end
      
      if ~obj.expH.hasperframesex || fast,
        sex = obj.expH.sex(fly);
        return;
      end
      
      if nargin < 3,
        sex = obj.trx(fly).sex;
        return;
      end
      
      sex = obj.trx(fly).sex(ts + obj.trx(fly).off);
      
    end
    
    %%%% FIXED %%%%
    function sex = GetSex1(obj,fly,t)
      % x = GetSex1(obj,expi,fly,t)
      % Returns the sex for the input experiment, SINGLE fly, and
      % SINGLE frame.
      
      if ~obj.expH.hassex,
        sex = '?';
        return;
      end
      
      if ~obj.expH.hasperframesex,
        sex = obj.sex(fly);
        if iscell(sex),
          sex = sex{1};
        end
        return;
      end
      
      sex = obj.trx(fly).sex{t + obj.trx(fly).off};
      
    end
    
    
  end
  
end