classdef JFeatures < handles
  
  properties (Access = public)
    featureconfigfile = '';
    allpf = {};
    curpf = {};
    projconfH = [];
    perframeGenerate = false;
    
    % constant: radius of window data to compute at a time
    windowdatachunk_radius = 500;
    perframe_params = {};
    landmark_params = {};
    
    % file containing feature parameters
    featureparamsfilename = 0;
    
    % parameters of window features, represented as a struct
    windowfeaturesparams = struct;
    
    % parameters of window features, represented as a cell array of
    % parameter name, parameter value, so that it can be input to
    % ComputeWindowFeatures
    windowfeaturescellparams = {};
    
    % State of the basic/compact feature table.
    basicFeatureTable = {};
    featureWindowSize = [];
    
    % per-frame features that are used
    allperframefns = {};
    curperframefns = {};
    perframeunits = {};
    
    % to overwrite or keep the perframe files.
    perframeOverwrite = [];
    
  end
  
  methods (Access = public,Static = true)
    function [X,feature_names] = ...
        ComputeWindowDataChunkStatic(curperframefns,allperframefns,perframefile,flies,windowfeaturescellparams,t0,t1)
      
      %function [X,feature_names] = ...
      %    ComputeWindowDataChunkStatic(perframefns,perframedir,flies,windowfeaturecellparams,t0,t1)
      %
      % Computes a chunk of windowdata between frames t0 and t1 for flies
      % flies.
      %
      % X is the nframes x nfeatures window data, feature_names is a cell array
      % of length nfeatures containing the names of each feature.
      %
      % It then loops through all the per-frame features, and calls
      % ComputeWindowFeatures to compute all the window data for that
      % per-frame feature.
      
      X = [];
      feature_names = {};
      
      for j = 1:numel(curperframefns),
        fn = curperframefns{j};
        ndx = find(strcmp(fn,allperframefns));
        
        perframedata = load(perframefile{ndx});
        perframedata = perframedata.data{flies(1)};
        
        t11 = min(t1,numel(perframedata));
        [x_curr,feature_names_curr] = ...
          ComputeWindowFeatures(perframedata,windowfeaturescellparams.(fn){:},'t0',t0,'t1',t11);
        
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
    
  end
  
  methods (Access = public)
    
    
    function [success,msg] = SetFeatureConfigFile(obj,configfile)
      success = false;
      msg = '';
      
      obj.featureConfigFile = configfile;
      [settings,~] = ReadPerFrameParams(configfile);
      obj.allperframefns =  fieldnames(settings.perframe);
      if isempty(obj.allperframefns)
        msg = 'No perframefns defined';
        return;
      end
      success = true;
      
    end
    
    function allpf = GetAllperframefns(obj)
      allpf = obj.allpf;
    end
    
    function curpf = GetCurperframefns(obj)
      curpf = obj.curpf;
    end
    
    function [filenames,timestamps] = GetPerframeFiles(obj,expH,dowrite)
      % [filenames,timestamps] = GetPerFrameFiles(obj,file,expi)
      % Get the full path to the per-frame mat files for experiment expi
      
      if nargin < 3,
        dowrite = false;
      end
      
      fn = obj.projconfH.GetFileName('perframedir');
      
      % if this is an output file, only look in output experiment directory
      if dowrite && JLabelData.IsOutputFile('perframedir'),
        expdirs_try = expH.GetOutexpdirs;
      else
        % otherwise, first look in output directory, then look in input
        % directory
        expdirs_try = {expH.GetOutexpdirs(),expH.GetExpdirs()};
      end
      
      filenames = cell(1,numel(obj.allpf));
      timestamps = -inf(1,numel(obj.allpf));
      
      for i = 1:numel(obj.allpf),
        
        % loop through directories to look in
        for j = 1:numel(expdirs_try),
          expdir = expdirs_try{j};
          filename = fullfile(expdir,fn,[obj.allpf{i},'.mat']);
          
          if exist(filename,'file'),
            filenames{i} = filename;
            tmp = dir(filename);
            timestamps(i) = tmp.datenum;
          elseif j == 1,
            filenames{i} = filename;
          end
          
        end
      end
    end
    
    function SetGenerateMissingFiles(obj)
      obj.perframeGenerate = true;
    end
    
    function perframeGenerate = GetGenerateMissingFiles(obj)
      perframeGenerate = obj.perframeGenerate;
    end
    
    function [success,msg] = GenerateMissingFiles(obj,expi,isInteractive)
      % [success,msg] = GenerateMissingFiles(obj,expi)
      % Generate required, missing files for experiments expi.
      % TODO: implement this!
      
      success = true;
      msg = '';
      
      if nargin<3
        isInteractive = true;
      end
      
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
              [success1,msg1] = obj.GeneratePerFrameFiles(expi,isInteractive);
              success = success && success1;
              if ~success1,
                msg = [msg,'\n',msg1]; %#ok<AGROW>
              end
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
    
    function [success,msg] = GeneratePerFrameFiles(obj,expi,isInteractive)
      success = false; %#ok<NASGU>
      msg = '';
      
      perframedir = obj.GetFile('perframedir',expi);
      
      if ~isInteractive,
        dooverwrite = false;
      elseif ~isempty(obj.perframeOverwrite)
        if obj.perframeOverwrite
          dooverwrite = true;
        else
          dooverwrite = false;
        end
      elseif exist(perframedir,'dir'),
        res = questdlg('Do you want to overwrite existing files or keep them?',...
          'Regenerate files?','Overwrite','Keep','Keep');
        dooverwrite = strcmpi(res,'Overwrite');
        obj.perframeOverwrite = dooverwrite;
      else
        dooverwrite = true;
      end
      
      expdir = obj.expdirs{expi};
      
      if isInteractive
        hwait = mywaitbar(0,sprintf('Initializing perframe directory for %s',expdir),'interpreter','none');
      else
        fprintf('Initializing perframe directory for %s\n',expdir);
      end
      
      perframetrx = Trx('trxfilestr',obj.GetFileName('trx'),...
        'moviefilestr',obj.GetFileName('movie'),...
        'perframedir',obj.GetFileName('perframedir'),...
        'landmark_params',obj.landmark_params,...
        'perframe_params',obj.perframe_params,...
        'rootwritedir',obj.rootoutputdir);
      
      perframetrx.AddExpDir(expdir,'dooverwrite',dooverwrite,'openmovie',false);
      
      perframefiles = obj.GetPerframeFiles(expi);
      for i = 1:numel(obj.allperframefns),
        fn = obj.allperframefns{i};
        %ndx = find(strcmp(fn,obj.allperframefns));
        file = perframefiles{i};
        if ~dooverwrite && exist(file,'file'),
          continue;
        end
        if isInteractive
          hwait = mywaitbar(i/numel(obj.allperframefns),hwait,sprintf('Computing %s and saving to file %s',fn,file));
        else
          fprintf('Computing %s and saving to file %s\n',fn,file);
        end
        perframetrx.(fn); %#ok<VUNUS>
        
      end
      
      if isInteractive && ishandle(hwait),
        delete(hwait);
      end
      
      success = true;
      
    end
    
    function [success,msg] = SetFeatureParamsFileName(obj,featureparamsfilename)
      % [success,msg] = SetFeatureParamsFileName(obj,featureparamsfilename)
      % Sets the name of the file describing the features to use to
      % featureparamsfilename. These parameters are read in. Currently, the
      % window data and classifier, predictions are not changed. (TODO)
      
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
      %       try
      [windowfeaturesparams,windowfeaturescellparams,basicFeatureTable,featureWindowSize] = ...
        ReadPerFrameParams(featureparamsfilename); %#ok<PROP>
      %       catch ME,
      %         msg = sprintf('Error reading feature parameters file %s: %s',...
      %           params.featureparamsfilename,getReport(ME));
      %         return;
      %       end
      obj.SetPerframeParams(windowfeaturesparams,windowfeaturescellparams); %#ok<PROP>
      obj.featureparamsfilename = featureparamsfilename;
      obj.basicFeatureTable = basicFeatureTable;
      obj.featureWindowSize = featureWindowSize;
      success = true;
    end
    
    function SetPerframeParams(obj,windowfeaturesparams,windowfeaturescellparams)
      obj.windowfeaturesparams = windowfeaturesparams; %#ok<PROP>
      obj.windowfeaturescellparams = windowfeaturescellparams; %#ok<PROP>
      obj.curperframefns = fieldnames(windowfeaturesparams);
    end
    
    function [windowfeaturesparams,windowfeaturescellparams] = GetPerframeParams(obj)
      windowfeaturesparams = obj.windowfeaturesparams; %#ok<PROP>
      windowfeaturescellparams = obj.windowfeaturescellparams; %#ok<PROP>
    end
    
    function ScoresToPerframe(obj,expH,fn)
      outdir = expH.GetOutexpdirs{expi};
      scoresFileIn = fullfile(outdir,fn);
      scoresFileOut = fullfile(outdir,expH.GetFileName('perframe'),fn);
      Q = load(scoresFileIn);
      OUT = struct();
      OUT.units = struct(); OUT.units.num = {'scores'};
      OUT.units.den = {''};
      for ndx = 1:numel(Q.allScores.scores)
        t0 = Q.allScores.tStart{ndx};
        t1 = Q.allScores.tEnd{ndx}-1;
        OUT.data{ndx} = Q.allScores.scores{ndx}(t0:t1);
      end
      save(scoresFileOut,'-struct','OUT');
    end
    
    function [success,msg,t0,t1,X,feature_names] = ComputeWindowDataChunk(obj,expi,flies,t,mode,forceCalc)
      % [success,msg,t0,t1,X,feature_names] = ComputeWindowDataChunk(obj,expi,flies,t)
      % Computes a chunk of windowdata near frame t for experiment expi and
      % flies flies. if mode is 'start', then the chunk will start at t. if
      % it is 'center', the chunk will be centered at t. if mode is 'end',
      % the chunk will end at t. by default, mode is 'center'.
      % t0 and t1 define the bounds of the chunk of window data computed. X
      % is the nframes x nfeatures window data, feature_names is a cell array
      % of length nfeatures containing the names of each feature.
      %
      % This function first chooses an interval of frames around t, depending
      % on the mode. it then chooses a subinterval of this interval that
      % covers all frames in this interval that do not have window data. This
      % defines t0 and t1.
      %
      % It then loops through all the per-frame features, and calls
      % ComputeWindowFeatures to compute all the window data for that
      % per-frame feature.
      %
      % To predict over the whole movie we use forceCalc which
      % forces the function to recalculate all the features even though they
      % were calculated before.
      
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
          res = questdlg(sprintf('Experiment %s is missing some perframe files. Generate now?',obj.expnames{expi}),'Generate missing files?','Yes','Cancel','Yes');
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
        function UpdatePerframeParams(obj,params,cellParams,basicFeatureTable,featureWindowSize)
    % Updates the feature params. Called by SelectFeatures
      obj.SetPerframeParams(params,cellParams)
      if nargin>2
        obj.basicFeatureTable = basicFeatureTable;
        obj.featureWindowSize = featureWindowSize;
      end
      obj.ClearWindowData();
      obj.classifier = [];
      obj.classifier_old = [];
      obj.PreLoadLabeledData();
      % TODO: remove clearwindow features.
    end

  end
  
end