classdef ScoreFile
  % Value class by design
  
  properties
    allScores 
    behaviorName = '';
    timestamp
    version
    jabFileNameAbs    
  end
  
  methods % setters
    function obj = set.allScores(obj,v)
      assert(isscalar(v) && isstruct(v),'.allScores must be a scalar struct.');
      obj.allScores = v;
    end
    function obj = set.behaviorName(obj,v)
      assert(ischar(v),'.behaviorName must be a string.');
      obj.behaviorName = v;
    end
    function obj = set.timestamp(obj,v)
      assert(isscalar(v) && isnumeric(v),'.timestamp must be a numeric scalar.');
      obj.timestamp = v;
    end
    function obj = set.jabFileNameAbs(obj,v)
      assert(ischar(v),'.jabFileNameAbs must be a string.');
      obj.jabFileNameAbs = v;
    end
  end
  
  methods
    
    function obj = ScoreFile(varargin)
      if nargin==1
        s = varargin{1};
        assert(isstruct(s) && isscalar(s));
        obj = obj.initFromStruct(s);
      end
    end
        
    function save(obj,filename)
      warnst = warning('off','MATLAB:structOnObject');
      s = struct(obj); %#ok<NASGU>
      warning(warnst);

      if exist(filename,'file')
        if exist([filename '~'],'file')
          delete([filename '~']);
        end
        
        [didbak,msg] = copyfile(filename,[filename,'~']);
        if ~didbak
          error('ScoreFile:unableToBackupScores', ...
            'Could not create backup of %s: %s. Aborting.',filename,msg);
        end
      end
      save(filename,'-struct','s');
      
%     AL 20150112: This was in JAABAST
%       for tryi = 1:5,
%         try
%           tmp = load(fname);
%         catch ME,
%           warning('Try %d: error testing whether we could reload the saved file %s: %s',tryi,fname,getReport(ME));
%           continue;
%         end
%         break;
%       end
    end

    function obj = initFromStruct(obj,s)
      fn = fieldnames(s);
      for f = fn(:)', f=f{1}; %#ok<FXSET>
        obj.(f) = s.(f);
      end
    end

  end
  
  methods (Static) % ScoreFile
        
    function obj = load(filename)
      assert(exist(filename,'file')==2,'File ''%s'' not found.',filename);
      sf = load(filename);
      
      if ~ScoreFile.isValidScoreMATFile(sf)
        error('ScoreFile:contents','''%s'' does not appear to be a valid score file.',filename);
      end
        
      obj = ScoreFile;
      obj = initFromStruct(obj,sf);    
    end
    
    function tf = isValidScoreMATFile(sf)
      % Verify MAT-file contents of score file
      REQD_FIELDS = {'allScores' 'timestamp' 'version' 'jabFileNameAbs'};
      tf = all(isfield(sf,REQD_FIELDS));
    end
    
    function n = defaultScoreFilename(classifiername)
      n = sprintf('scores_%s.mat',classifiername);
    end
    
  end
  
  methods (Static) % allScores substructure
    
    function allScrs = allScrs(nfly)
      % allScores contstructor. allScores represents all scores for 1
      % experiment, all flies.
      
      allScrs = struct(...
        'scores',{cell(1,nfly)},... % cell vec of timeseries
        'tStart',nan(1,nfly),... 
        'tEnd',nan(1,nfly),...
        'postprocessed',{cell(1,nfly)},...
        'postprocessparams',[],...
        't0s',{cell(1,nfly)},...
        't1s',{cell(1,nfly)},...
        'scoreNorm',nan);
    end
    
    function allScrs = AllScrsInitFromPredData(allScrs,predDataExp,iCls,firstFrms,endFrms)
      % Sets: .scores, .tStart, .tEnd, .postprocessed, .t0s, .t1s from
      % predDataExp{iFly}(iCls).cur, .cur_pp
      % NOT SET: .postprocessparams, .scoreNorm 
      
      nfly = numel(predDataExp);
      assert(isequal(nfly,numel(allScrs.scores),numel(firstFrms),numel(endFrms)));
      
      for fly = 1:nfly
        pd = predDataExp{fly}(iCls);
        assert(all(pd.cur_valid));
        
        curt = pd.t;
        assert(all(curt(2:end)-curt(1:end-1)==1),'Scores are out of order.');
        
        tStart = firstFrms(fly);
        tEnd = endFrms(fly);
        scores = nan(1,tEnd);
        pp = nan(1,tEnd);
        scores(tStart:tEnd) = pd.cur;
        pp(tStart:tEnd) = pd.cur_pp;
        
        allScrs.scores{fly} = scores;
        allScrs.tStart(fly) = tStart;
        allScrs.tEnd(fly) = tEnd;
        allScrs.postprocessed{fly} = pp;
        [i0s,i1s] = get_interval_ends(allScrs.postprocessed{fly}>0);
        allScrs.t0s{fly} = i0s;
        allScrs.t1s{fly} = i1s;
      end
    end
    
  end
  
end