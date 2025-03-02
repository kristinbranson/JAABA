classdef JLabelData < matlab.mixin.Copyable
  % This is the class that is essentially the JLabel "model" in MVC terms 
  % --- it holds all of the critical data that JLabel allows a user to
  % manipulate, but doesn't deal with visual or presentational issues.
  % (OK, in a few places it does, but those should probably be changed.)
  
  % About everythingParams and basicParams: these are typically instances
  % of class Macguffin.  Macguffins represent all of the stuff that gets 
  % saved to the .jab file.  A variable named basicParams is typically a
  % Macguffin that has no experiments, or labels, or classifier set.  They
  % typically arise when creating a new .jab file.  A variable named
  % everythingParams is generally the contents of a .jab file, which may
  % well have experiments, labels, and/or a classifier.
  
  % About lexicons, sublexicons, dialects, subdialects, vocabularies, and 
  % all that: 
  %   The feature lexicon is the set of all possible features, not 
  %   including the score features (a.k.a. scores-as-inputs).
  %
  %   The feature dialect is the feature lexicon plus the score features.
  %   (In this analogy, the score features are like region-specific words,
  %   or jargon.)
  % 
  %   The feature sublexicon is the subset of the feature lexicon that
  %   actually gets calculated for a particular .jab file.
  % 
  %   The feature subdialect is the sublexicon plus the score features.
  % 
  %   The feature vocabulary is the enabled features of the sublexicon, plus
  %   the enabled score features.  (The idea is that the classifier only "speaks"
  %   this vocabulary.)
  %
  %   Note that in all of the above, "features" can mean "per-frame
  %   features" or "window features", depending on the context.
  
  % -----------------------------------------------------------------------
  properties (SetAccess=private, GetAccess=public)
    expi  % currently selected experiment
    flies  % currently selected flies
  end
  
  % -----------------------------------------------------------------------
  properties (Access=public) 
    % type of target (mainly used for plotting)
    % this should be a type of animal, and should be singular, not plural
    targettype
        
    % last-used trajectories (one experiment, all flies)
    trx
    
    % last-used per-frame data (one fly)
    % This seems to generally be a cell array with as many elements as
    % allperframefns.  (I.e. as many features as are in the subdialect.)
    % Each element contains a double array with approximately as many
    % elements as there are frames in the current track, which gives the
    % value of that per-frame feature for that frame, for the current
    % experiment and target.  --ALT, Apr 10 2013
    % 
    % perframedata is common to all classifiers.
    perframedata
    
    % A feature lexicon is the universe of possible features available.
    % This is common to all classifiers. The features used when training
    % a classifier for a particular behavior are a subset of the feature
    % lexicon, and the features used by a particular trained classifier
    % will be a further subset. Each possible feature lexicon has a name,
    % such as 'flies', or 'vivek_mice'
    featureLexiconName
    
    % The actual feature lexicon itself. Note that we should not assume
    % the feature lexicon exactly matches the feature lexicon named
    % featureLexiconName in the global config files.  That's because the
    % user may have modified the featureLexicon for this everything file.
    % (Granted, there's no way to do this in the UI as of Apr 2, 2013, but
    % that may change in the future.  Plus we want to at least make it
    % possible for the user to hack the everything file in this way.)
    featureLexicon
    
    % Logical scalar, true for ST mode
    isST
    
    % windowdata holds computed and cached window features
    % windowdata(iCls).X
    % ...
    % Future optimization: window data computation for multiple
    % classifiers 
    windowdata
    
    
    % selFeatures keeps track of features selected when the user trains a
    % full classifier.
    selFeatures

    % predictdata stores predictions.  It is cell array, with as many
    % elements as there are experiments.  Each element holds another cell
    % array, with as many elements as there are targets in that experiment.
    % Each element of _that_ cell array is a struct array, with fields:
    %
    %   t           - startframe:endframe for exp, target. ALTODO: This state is duplicated identically across multiple classifiers
    %   cur
    %   cur_valid
    %   cur_pp      - This and other *_pp (postprocessed) fields are logical vectors
    %   old
    %   old_valid
    %   old_pp
    %   loaded
    %   loaded_valid
    %   loaded_pp
    %   timestamp
    %
    % I think the bottom-level structs are either empty or scalar, and are
    % only empty if the predictions have not been calculated for that
    % experiment and animal.
    %
    % In the scalar case, all these fields hold row vectors of length equal
    % to the trajectory length in frames. cur holds the scores generated by
    % the current classifier, old holds the scores generated by the old
    % classifier, loaded holds scores that were loaded from a file. _pp
    % fields hold post-processed scores, _valid fields hold whether that
    % score for that experiment, target, and frame is valid, and fields
    % without a suffix hold the raw scores (I think).  t gives the frame
    % number for each element, and timestamp gives the timestamp of the
    % classifier that generated the score.
    %
    % AL: predictdata{expi}{flyi} = [] OR
    % predictdata{expi}{flyi}(timelinei). Currently timelinei corresponds
    % to behaviori/labeli for 1:numRealBehaviors.
    predictdata
    
    % The predictblocks stores which frames, for which tracks, get
    % predicted when training happens.  (B/c you only really need to
    % predict frames where there are labels.)  Also, if a user is looking
    % at a particular stretch of frames for a particular track, that
    % stretch will get added to predictblocks, regardless of whether there
    % are labelled frames anywhere near it.  -- ALT, Mar 04 2013
    %
    % predictblocks(iCls).t0
    % predictblocks(iCls).t1
    % predictblocks(iCls).expi
    % predictblocks(iCls).flies
    predictblocks

    % See Predict.fastPredict
    % fastPredict(iCls).<etc>
    fastPredict
        
    % constant: radius of window data to compute at a time. Currently
    % shared by all classifiers.
    windowdatachunk_radius
    predictwindowdatachunk_radius
    
    % labels struct array. See also Labels.m
    % labels(expi) is the labeled data for experiment expi
    % labels(expi).t0s are the start frames of all labeled sequences for
    % experiment expi
    % labels(expi).t1s are the corresponding end frames of all labeled
    % sequences for experiment expi (strictly, one past the last frame)
    % labels(expi).names is the cell array of the corresponding behavior
    % names for all labeled sequences for experiment expi
    % labels(expi).flies is the nseq x nflies_labeled matrix of the
    % corresponding flies for all labeled sequences for experiment expi
    % t0s{j}, t1s{j}, names{j}, and flies(j,:) correspond to each other. 
    % labels(expi).off is the offset so that labels(expi).t0s(j) +
    % labels(expi).off corresponds to the frame of the movie (since the
    % first frame for the trajectory(s) may not be 1.
    % labels(expi).timestamp is the Matlab timestamp at which labels(expi)
    % was last set
    labels
    
    % labels for the current experiment and flies, represented as an array
    % such that labelidx(t+labelidx_off) is the index of the behavior for
    % frame t of the movie. labelidx(i) == 0 corresponds to
    % unlabeled/unknown, otherwise labelidx(i) corresponds to behavior
    % labelnames{labelidx{i})
    % See Labels.m
    labelidx
    labelidx_off % AL20150227: can prob make this Dependent
    
    % first frame that all flies currently selected are tracked
    t0_curr
    % last frame that all flies currently selected are tracked
    t1_curr
    
    % predicted label for current experiment and fly, with the same type
    % of representation as labelidx, except that a value of zero means that
    % there is no prediction (as when a classifier has not yet been trained, 
    % or has been cleared).  These variables are essentially a cache of the
    % data in self.predictdata.
    predictedidx % Takes values in {0,1,1.5,2}; see GetPredictedIdx
    scoresidx  % Is this really an index?  Isn't it just the score itself?
    scoresidx_old
    scoreTS  % timestamp for the above. ALTODO: does not appear to be used for anything
        
    % names of behaviors, corresponding to labelidx
    labelnames        
  end  
  
  properties
    ntimelines
    iLbl2iCls % 2*nclassifiers-by-1 array. iCls = iLbl2iCls(iLbl) where 
              % iLbl/iCls reference .labelnames/.classifiernames resp.
    iCls2iLbl % nclassifiers-by-1 cell array, each el is 1-by-2 array.
              % [iLblPos iLblNeg] = iCls2iLbl{iCls} where iLblPos, iLblNeg 
              % index .labelnames.

% AL20150303: No usages of labelstats right now, but more importantly 
% seems like this can be lazily computed              
%     % statistics of labeled data per experiment
%     % labelstats(expi).nflies_labeled is the total number of flies labeled,
%     % labelstats(expi).nbouts_labeled is the total number of bouts of
%     % behaviors labeled, labelstats(expi).
%     labelstats
    
    % computing per-frame properties
    perframe_params
    landmark_params
    
    classifiertype % 1-by-nclassifiers cellstr    
    % currently learned classifier. structure depends on the type of
    % classifier. if empty, then no classifier has been trained yet. 
    classifier % 1-by-nclassifiers cell
    classifier_old
    lastFullClassifierTrainingSize
    classifierTS  % 1-by-nclassifiers vector of default time stamps. Number of days since Jan 0, 0000 (typically non-integer)    
    trainstats % 1-by-nclassifiers cell array
    
    % parameters to learning the classifier. struct fields depend on type of classifier.
    % 1-by-nclassifiers cell
    classifier_params
        
    % constant: files per experiment directory
    %filetypes = {'movie','trx','label','gt_label','perframedir','clipsdir','scores'};
    filetypes
    
    % locations of files within experiment directories
    moviefilename
    movieindexfilename
    trxfilename
    scorefilename % cellstr
    perframedir
    clipsdir
    scores
    stfeatures
    trkfilename
    
    % Properties relating to whether there is a movie to show
    %openmovie  % true iff a movie is one of the required files for each experiment
    %ismovie  
    
    % experiment info: expi indexes the following
    
    % cell array of input experiment directory paths for the current GT
    % mode
    expdirs
    
    % array of number of flies in each experiment
    nflies_per_exp
    
    % cell array of arrays of first frame of each trajectory for each
    % experiment: firstframes_per_exp{expi}(fly) is the first frame of the
    % trajectory of fly for experiment expi. 
    firstframes_per_exp
    % cell array of arrays of end frame of each trajectory for each
    % experiment: endframes_per_exp{expi}(fly) is the last frame of the
    % trajectory of fly for experiment expi. 
    endframes_per_exp
    
    % sex per experiment, fly
    frac_sex_per_exp
    sex_per_exp
    
    % whether sex is computed
    hassex
    % whether sex is computed on a per-frame basis
    hasperframesex
    
    % last-used path for loading jabfiles
    defaultpath

    % last-used path for loading experiments
    expdefaultpath

    % windowfeaturesparams{iCls} (row vec) is:
    % struct from curperframefns->[structs containing window feature parameters]
    % Each field holds the parameters for a single per-frame feature in the
    % feature vocabulary, with the field name being the per-frame feature
    % name. Score features are included in the feature vocabulary.
    windowfeaturesparams
    
    % windowfeaturescellparams{iCls} (row vec) is:
    % struct from curperframefns->[cell array of window feature parameters]
    % parameters of window features, represented as a cell array of
    % parameter name, parameter value, so that it can be input to
    % ComputeWindowFeatures
    windowfeaturescellparams
    
    %savewindowfeatures = false;  % whether unsaved changes have been made
    %                             % to the features
    
    % State of the basic/compact feature table.
    %basicFeatureTable = {};  % now calculate on the fly
    %maxWindowRadiusCommonCached = [];  
      % need to remember between calls to SelectFeatures, because it needs
      % to override the max_window_radius in the window-feature amount
      % presets in the feature lexicon, and these WF amount presets are not
      % retained in the feature vocabulary.  (And we don't want to retain
      % them in JLabelData's feature vocabulary, b/c they're not _really_
      % part of the feature vocabulary.)  (But I suppose we could add them
      % if we wanted to...)
    
    % per-frame features that are used
    allperframefns  % The list of all per-frame feature names in 
                    % the lexicon that are actually calculated for
                    % each frame, plus the scores-as-input feature
                    % names.  I would call this the the 'subdialect'.
                    % --ALT, Apr 5, 2013. Common to all classifiers.
    curperframefns  % The list of all per-frame feature names
                    % that are enabled, i.e. used for used for classifier 
                    % training.  This is a subset of allperframefns, and will 
                    % include a subset of the score features (but the
                    % subset might be the empty set).  I would 
                    % call this the 'vocabulary'.  --ALT, May 18, 2013.
                    % Indexed by classifier, ie curperframefns{iCls} (row vec).
      % Thus curperframefns is the subset of allperframefns that are enabled.                    
      
    % units for perframedata
    perframeunits 

    % the scores from other classifiers that that used as features for this
    % classifier.  A structure array with number of elements equal to the
    % number of score features.  classifierfile holds the absolute path of
    % the .jab file holding the external classifier, ts holds the time step
    % of the classifier, and scorefilename is the local base name of the
    % score file stored in the experiment directory (e.g. "Chase_v7")
    scoreFeatures 
    
    % experiment/file management

    % matrix of size nexps x numel(file_types), where
    % fileexists(expi,filei) indicates whether file filetypes{filei} exists
    % for experiment expi
    fileexists 
    
    % timestamps indicating time the files were last edited, same structure
    % as fileexists
    filetimestamps
    
    % whether all necessary files for all experiments exist
    allfilesexist 

    % whether we can generate any missing files
    filesfixable 
    
    % whether user has given permission to generate the perframe files
    perframeGenerate
    
    % to overwrite or keep the perframe files.
    perframeOverwrite
    
    % warn about removing arena features.
    arenawarn
    hasarenaparams
    
    % functions for writing text to a status bar
    setstatusfn
    clearstatusfn
    
    % Show similar frames
    % These have not been updated for multiple classifiers.
    % SimilarFrames/Bagging functionality should all be asserted to run
    % only on single-classifier projects.
    frameFig 
    distMat 
    bagModels
    fastPredictBag 
    % fastPredictBag.classifier
    % fastPredictBag.windowfeaturescellparams
    % fastPredictBag.wfs
    % fastPredictBag.pffs
    % fastPredictBag.ts
    % fastPredictBag.tempname
    % fastPredictBag.curF
    % fastPredictBag.dist{expi}{fly}
    % fastPredictBag.trainDist

    confThresholds % nclassifiers-by-2 array
    
    predictOnlyCurrentFly % nclassifiers logical vec, predict only for current fly.
    
    % Retrain properly
    doUpdate % ALTODO: does not appear to be used
    
    % Ground truthing or not
    gtMode        % Seems like it would make sense to not have
                  % JLabelData know about this---it's really more
                  % a property of the View, not the Model, it seems to
                  % me.  --ALT, Jan 17 2013
                  % Maybe, but making it that way would be a pain in the
                  % ass.  --ALT, Apr 22 2013
                  % This is set to empty unless a file is open.
    
    % Ground truthing suggestion
    randomGTSuggestions % randomGTSuggestions{iExp}(iFly).start, .end (scalars)
    thresholdGTSuggestions
    loadedGTSuggestions % same format/structure as randomGTSuggestions, except .start/.end can be vectors
    balancedGTSuggestions
    GTSuggestionMode
    
    cacheSize % ALTODO comments suggests obsolete
    savewindowdata % 1-by-nclassifiers double
    loadwindowdata % 1-by-nclassifiers double. 
    
    postprocessparams % 1-by-nclassifiers cell
    version
        
    % A place to store things for the other mode.
    % So if the JLabelData instance is created in GT mode, this stores the
    % normal experiment directory names, labels, and label statisitics.  If
    % JLabelData is created in Normal mode, this store the GT experiment
    % directory names, etc.  This keeps the stuff for the other mode out of
    % our hair, but keeps it around so that we can save it to the
    % everything file.
    otherModeLabelsEtc
    
    % you could argue that these are view-related, and so shouldn't be in
    % here, but they get saved to the everything file, so we'll include
    % them here.
    trxGraphicParams
    labelcolors
    unknowncolor
    
    % A slot to store whether the JLabelData object is being used
    % interactively.  Controls whether the object tries to ask the user
    % questions, etc.  At some point, it might be nice to refactor so that
    % this went away, and the JLabelData object _never_ did UI type stuff.
    isInteractive
    
    % .jab-file handling stuff
    thereIsAnOpenFile
    everythingFileNameAbs      % the name of the everything file, if one
                               % is open.  We need this here b/c a new
                               % everything file doesn't have a JLabelData
                               % object yet.
    userHasSpecifiedEverythingFileName         % true iff the everything
                                               % file name was specified by
                                               % the user, as opposed to
                                               % being chosen by default
                                               % when a new file was
                                               % created
    needsave  % true iff a file is open and there are unsaved changes.  False if no file is open.
    
    usePastOnly % whether to only use past information when predicting the current frame
    
    deterministic = false; % scalar logical or double, for testing purposes. If nonzero, may be used as RNG seed
    
    trainWarnCount = 0; % number of training iterations since we warned the user about increasing the iterations.
    trainWarn = true;
    
    % Details of APT project.
    aptInfo = struct()
    fromAPT = false;
    
    % Details of st
    stFeatures = false;
    stInfo = [];
    %stInfo = getSTParams();
  end

 
  % -----------------------------------------------------------------------
  properties (GetAccess=public,SetAccess=immutable,Dependent=true)
    expnames
    nexps
    nTargetsInCurrentExp
    ismovie    % true iff the movie file name is nonempty.  If movie file name is empty, it means we don't try to open movies.
    nbehaviors % TO BE DEPRECATED in favor of nlabelnames.
    nlabelnames % number of labels, including 'none' or <no-behaviors>; equivalently, numel(obj.labelnames)
    nclassifiers
    isMultiClassifier
    classifiernames % 1-by-nclassifiers cellstr
    nobehaviornames % 1-by-nclassifiers cellstr  
    iCls2LblNames % nclassifiers-by-1 cell array, each el is 1-by-2 cellstr
    
    % AL20150304: Currently has no clients; NextJump computes a similar qty
    % on its own
    %
    % whether the predicted label matches the true label. 0 stands for
    % either not predicted or not labeled, 1 for matching, 2 for not
    % matching. (This is the same representation as labelidx for
    % single classifiers.)
    erroridx
  end

  
  %% Getters/Setters
  
  methods
    
    function v = get.nclassifiers(self)
      v = self.ntimelines;
    end
    
    function v = get.nbehaviors(self)
      v = numel(self.labelnames);
    end
    
    function v = get.nlabelnames(self)
      v = numel(self.labelnames);
    end
    
    function v = get.isMultiClassifier(self)
      v = self.nclassifiers>1;
    end
    
    function v = get.classifiernames(self)
      v = self.labelnames(1:self.nclassifiers);
    end
    
    function v = get.nobehaviornames(self)
      v = self.labelnames(self.nclassifiers+1:end);
    end
    
    function v = get.iCls2LblNames(self)
      lblnames = self.labelnames;
      iCls2iLbl = self.iCls2iLbl;
      v = cellfun(@(x)lblnames(x),iCls2iLbl,'uni',0);
    end
    
    function expnames = get.expnames(self)
      % Get the names of all the experiments for the current GT mode
      expDirNames = self.expdirs;
      expnames = cellfun(@fileBaseName,expDirNames,'UniformOutput',false);
    end    
    
    function nexps = get.nexps(self)
      % Get the number of experiments for the current GT mode
      nexps = length(self.expdirs);
    end
    
    function nTargetsInCurrentExp = get.nTargetsInCurrentExp(self)
      if self.nexps>0 && ~isempty(self.expi) && self.expi~=0 ,
        nTargetsInCurrentExp=self.nflies_per_exp(self.expi);
      else
        nTargetsInCurrentExp=[];
      end
    end    
    
    function ismovie = get.ismovie(self)
      % Get whether the movie file name is set
      ismovie = ~isempty(self.moviefilename);
      % ismovie=~isempty(self.moviefilename) && self.openmovie;
    end
      
    function v = get.erroridx(obj)
      % Calculates erroridx based on obj.predictedidx and obj.labelidx.
      
      %MERGESTREVIEWED

      if obj.expi == 0
        % legacy codepath
        v = [];
        return;
      end
      
      n = obj.t1_curr - obj.t0_curr + 1;
      nTL = obj.labelidx.nTL;      
      
      predIdx = obj.PredictedIdxExpandValues(obj.predictedidx);
      lblIdx = obj.labelidx.vals;
      assert(isequal(size(predIdx),size(lblIdx),[nTL n]));
      idxcurr = predIdx~=0 & lblIdx~=0;
      
      v = zeros(nTL,n);
      v(idxcurr) = double(predIdx(idxcurr)~=lblIdx(idxcurr))+1;      
    end
    
    
  end
 
  
  %% Project Config
  
  methods
    
    % ---------------------------------------------------------------------
    function setWindowFeaturesParams(obj,windowFeaturesParams)
      % Updates the feature params.  Called after user clicks Done in 
      % Select Features...
      
      % MERGESTUPDATED
      
      if iscell(windowFeaturesParams)
        assert(numel(windowFeaturesParams)==obj.nclassifiers);
      else
        % convenience API, 'scalar expansion'
        assert(isstruct(windowFeaturesParams) && isscalar(windowFeaturesParams));
        windowFeaturesParams = repmat({windowFeaturesParams},1,obj.nclassifiers);
      end
            
      obj.windowfeaturesparams = windowFeaturesParams;
      obj.windowfeaturescellparams = cellfun(...
        @JLabelData.convertParams2CellParams,windowFeaturesParams,'uni',0);
      obj.curperframefns = cellfun(@fieldnames,windowFeaturesParams,'uni',0);
      oldScoreNorm = {obj.windowdata.scoreNorm};

      obj.clearClassifierProper();
      obj.initWindowData();
      obj.SetWindowFeatureNames();
%       [success,msg]=obj.PreLoadPeriLabelWindowData();
%       if ~success, 
%         error('JLabelData:unableToLoadPerLabelWindowData',msg);
%       end
      obj.needsave = true;
      
      if numel(oldScoreNorm)==obj.nclassifiers
        for iCls = 1:obj.nclassifiers
          if obj.HasLoadedScores(iCls) 
              obj.windowdata(iCls).scoreNorm = oldScoreNorm{iCls};
          end
        end
      end
    end
    
    
    % ---------------------------------------------------------------------
    function renameBehavior(obj,oldBeh,newBeh)
      % oldBeh/newBeh: char, old/new behavior name
      %
      % This method operates atomically for now, eg if there is a throw it
      % will occur before obj is modified.
    
      % Substance here is to update .labelnames, .labels, .labelidx,
      % .otherModeLabelsEtc.labels.
      
      assert(ischar(oldBeh) && ischar(newBeh));

      iCls = find(strcmp(oldBeh,obj.classifiernames));
      if ~isscalar(iCls)
        error('JLabelData:behaviorNotPresent',...
          'Behavior ''%s'' is not present in the project.',oldBeh);
      end
      if strcmp(newBeh,oldBeh)
        % no change
        return;
      end
      if ~JLabelData.isValidBehaviorName(newBeh)
        error('JLabelData:invalidBehaviorName',...
          '''%s'' is an invalid behavior name.',newBeh);
      end
      if any(strcmp(newBeh,obj.labelnames))
        error('JLabelData:duplicateName',...
          'The project already contains a behavior with name ''%s''.',newBeh);
      end
      
      % update .labelnames, .labels, .labelidx
      nCls = obj.nclassifiers;
      oldNoBeh = obj.labelnames{iCls+nCls};
      newNoBeh = Labels.noneOrNoBehaviorName(newBeh,nCls);
      if isfield(obj.labelidx,'labelnames'),
        assert(isequal(obj.labelidx.labelnames,obj.labelnames));
      end
      obj.labelnames{iCls} = newBeh;
      obj.labelnames{iCls+nCls} = newNoBeh;
      if isfield(obj.labelidx,'labelnames'),
        obj.labelidx.labelnames = obj.labelnames;
      end
      obj.labels = Labels.renameBehavior(obj.labels,oldBeh,newBeh,oldNoBeh,newNoBeh);
      
      if ~isempty(obj.otherModeLabelsEtc.labels)
        warning('JLabelData:untestedCodepath',...
          'Untested functionality: Updating other-mode labels.');
          obj.otherModeLabelsEtc.labels = Labels.renameBehavior(...
            obj.otherModeLabelsEtc.labels,oldBeh,newBeh,oldNoBeh,newNoBeh);
      end
      
      obj.needsave = true;        
    end
        
        
    % ---------------------------------------------------------------------
    function n = getBehaviorNames(obj)
      n = obj.classifiernames;
    end
    
    
    function tf = isBehaviorName(obj,name)
      % as opposed to no-behavior name
      tf = any(strcmp(name,obj.classifiernames));
    end
    
    
    function tf = isNoBehaviorName(obj,name)
      tf = any(strcmp(name,obj.nobehaviornames));      
    end
    
    
    function iCls = classifierIndexForName(obj,name)
      tf = strcmp(name,obj.labelnames);
      assert(nnz(tf)==1);
      iCls = obj.labelidx.idxBeh2idxTL(tf);
    end

    
    function noneToNoBeh(obj)
      % For single-classfier projects: Convert 'None' to 'No_<behavior>'
      
      assert(obj.nclassifiers==1);
      
      oldNoBeh = obj.labelnames{2};
      newNoBeh = Labels.noBehaviorName(obj.labelnames{1});
      assert(strcmp(oldNoBeh,'None'));
      obj.labelnames{2} = newNoBeh;
      
      obj.labels = Labels.renameBehaviorRaw(obj.labels,oldNoBeh,newNoBeh);
    end
    
    
    function addClassifier(obj,classifierName)
      tf = strcmp(classifierName,obj.classifiernames);
      if any(tf)
        error('JLabelData:addClassifier','Classifier with name ''%s'' already exists.',...
          classifierName);
      end
      
      obj.StoreLabelsForCurrentAnimal(); % store label info for current fly/expi from .labelIdx -> .labels
      if obj.nclassifiers==1
        obj.noneToNoBeh();
      end
      
      % Get some state before we start modifying props
      origNCls = obj.nclassifiers;
      origLabelNames = obj.labelnames; % not quite original if origNCls==1
%       behnobeh = obj.iCls2LblNames;
%       iLblsRm = obj.iCls2iLbl;
      
      noBehName = Labels.noBehaviorName(classifierName);
      obj.labelnames = [origLabelNames(1:origNCls) {classifierName} origLabelNames(origNCls+1:end) {noBehName}];
      sfn = ScoreFile.defaultScoreFilename(classifierName);
      obj.scorefilename{end+1} = sfn;
      % mapping from old label idxs -> new label idxs
      [~,oldLblIdx2NewLblIdx] = ismember(origLabelNames,obj.labelnames);

      [obj.ntimelines,obj.iLbl2iCls,obj.iCls2iLbl] = Labels.determineNumTimelines(obj.labelnames);
      obj.labelcolors = Labels.augmentColors(obj.labelcolors,numel(obj.labelnames),'lines');
      
      % Labels
      obj.labels = Labels.addClassifier(obj.labels,classifierName,noBehName);
      obj.reinitLabelIdx();
      
      % data/features 
      obj.windowdata(end+1,1) = WindowData.windowdata(1);
      obj.windowdata = WindowData.windowdataRemapLabelIdxs(obj.windowdata,...
        oldLblIdx2NewLblIdx);
      FLDS = {
        'windowfeaturesparams' 
        'windowfeaturescellparams' 
        'curperframefns' 
        'savewindowdata' 
        'loadwindowdata'
        'predictOnlyCurrentFly'
        }';
      for f = FLDS,f=f{1}; %#ok<FXSET>
        val = obj.(f);
        if iscell(val)
          assert(all(cellfun(@(x)isequaln(x,val{1}),val)));
        else
          assert(all(val(1)==val));
        end
      end
      for f = FLDS,f=f{1}; %#ok<FXSET>
        val = obj.(f);
        obj.(f)(end+1) = val(1); % works for both cells and arrs
      end
      
      obj.SetWindowFeatureNames();
            
      % predictions
      obj.PredictDataAddClassifier();      
      obj.predictblocks(end+1,1) = Predict.predictblocks(1);
      % obj.fastPredict -- see call to FindFastPredictParams below
      obj.UpdatePredictedIdx();
      
      % classifier
      cs = ClassifierStuff;
      obj.classifiertype{1,end+1} = cs.type;
      obj.classifier{1,end+1} = cs.params;
      obj.classifier_old{1,end+1} = [];
      obj.classifier_params{1,end+1} = cs.trainingParams;
      obj.classifierTS(1,end+1) = cs.timeStamp;
      obj.confThresholds(end+1,:) = cs.confThresholds;
      obj.windowdata(end).scoreNorm = cs.scoreNorm;
      obj.postprocessparams{1,end+1} = cs.postProcessParams;
      obj.savewindowdata(1,end+1) = cs.savewindowdata; % set this above
      % obj.loadwindowdata % see above      
      obj.trainstats{1,end+1} = [];
      obj.predictOnlyCurrentFly(end+1) = cs.predictOnlyCurrentFly;
      obj.selFeatures(end+1) = cs.selFeatures;
      FLDS = {
        'do' 
        'use' 
        }';
      for f = FLDS,f=f{1}; %#ok<FXSET>
        obj.selFeatures(end).(f) = obj.selFeatures(1).(f);
      end
      for f = FLDS,f=f{1}; %#ok<FXSET>
        val = [obj.selFeatures.(f)];
        assert(all(val(1)==val));
      end
      if obj.selFeatures(1).use,
        for ndx = 1:numel(obj.selFeatures),
          obj.selFeatures(ndx).do = true;
        end
      end
      
      if ~obj.isST
        obj.FindFastPredictParams();
      end
      
      obj.needsave = true;
    end
    
    
    function rmClassifier(obj,classifierName)
      tfRm = strcmp(classifierName,obj.classifiernames);
      assert(nnz(tfRm)<2,'Duplicate existing classifiernames.');
      if nnz(tfRm)==0
        error('JLabelData:rmClassifier','No classifier with name matching ''%s''.',...
          classifierName);
      end
      
      assert(obj.nclassifiers>1,'Cannot remove only classifier.');

      % Get some state before we start modifying props
      origLabelNames = obj.labelnames;
      behnobeh = obj.iCls2LblNames{tfRm};
      iLblsRm = obj.iCls2iLbl{tfRm};
      obj.StoreLabelsForCurrentAnimal(); % store label info for current fly/expi from .labelIdx -> .labels
      
      obj.scorefilename(tfRm) = [];
      obj.labelnames(iLblsRm) = [];
      % mapping from old label idxs -> new label idxs
      [~,oldLblIdx2NewLblIdx] = ismember(origLabelNames,obj.labelnames);
      if numel(obj.labelnames)==2
        oldNoBehavior = obj.labelnames{2};
        newNoBehavior = 'None';
        obj.labelnames{2} = newNoBehavior;
      end
      [obj.ntimelines,obj.iLbl2iCls,obj.iCls2iLbl] = Labels.determineNumTimelines(obj.labelnames);
      obj.labelcolors(iLblsRm,:) = [];      
      
      % Labels      
      obj.labels = Labels.removeClassifier(obj.labels,behnobeh{:});
      if numel(obj.labelnames)==2
        obj.labels = Labels.renameBehaviorRaw(obj.labels,oldNoBehavior,newNoBehavior);
        
        % ALXXX to be fixed
        warning('JLabelData:rmClassifier',...
          'Going from multi-classifier to single-classifier project. Label importance may not be set correctly for GT mode.');
      end
      obj.reinitLabelIdx();
      
      % data/features 
      obj.windowdata(tfRm,:) = [];
      obj.windowdata = WindowData.windowdataRemapLabelIdxs(obj.windowdata,...
        oldLblIdx2NewLblIdx);
      obj.windowfeaturesparams(tfRm) = [];
      obj.windowfeaturescellparams(tfRm) = [];
      obj.curperframefns(tfRm) = [];
      obj.savewindowdata(tfRm) = [];
      obj.loadwindowdata(tfRm) = [];
      obj.predictOnlyCurrentFly(tfRm) = [];
      
      % predictions
      nExp = numel(obj.predictdata);
      for iExp = 1:nExp
        nFly = numel(obj.predictdata{iExp});
        for iFly = 1:nFly
          obj.predictdata{iExp}{iFly}(tfRm) = [];
        end
      end
      obj.predictblocks(tfRm) = [];
      obj.fastPredict(tfRm) = [];
      obj.predictedidx(tfRm,:) = [];
      obj.scoresidx(tfRm,:) = [];
      obj.scoresidx_old(tfRm,:) = [];
      obj.scoreTS(tfRm,:) = [];
      %obj.erroridx(tfRm,:) = [];
      
      % classifier
      obj.classifiertype(tfRm) = [];
      obj.classifier(tfRm) = [];
      obj.classifier_old(tfRm) = [];
      obj.classifierTS(tfRm) = [];
      obj.trainstats(tfRm) = [];
      obj.classifier_params(tfRm) = [];
      obj.postprocessparams(tfRm) = [];
      obj.confThresholds(tfRm,:) = [];
      obj.selFeatures(tfRm) = [];
      if obj.selFeatures(1).use,
        for ndx = 1:numel(obj.selFeatures),
          obj.selFeatures(ndx).do = true;
        end
      end
      
      obj.needsave = true;
    end
    
    
    function reorderClassifiers(obj,behaviorNames)
      
    end
    
     
  end
  
  methods (Static)
    
    % ---------------------------------------------------------------------
    function result = isValidBehaviorName(behaviorName)
      result = ~isempty(regexp(behaviorName,'^[a-zA-Z_0-9]+$','once')) && ...
        isvarname(behaviorName); 
      % AL: second condition because labels.timelinetimestamp, see Labels.m
    end
    
  end
  
  methods % more private
    
    % ---------------------------------------------------------------------    
    function setFeatureSublexicon(obj,featureLexicon,featureLexiconName,sublexiconPFNames)
      % This sets the feature lexicon to the given one. If the
      % featureLexicon is not one of the named ones, then either no
      % featureLexiconName should be given, or the name should be 'custom'.
            
      obj.featureLexiconName = featureLexiconName;
      obj.featureLexicon = featureLexicon;
      
      % Update obj.perframe_params based on the new feature lexicon
      if isfield(featureLexicon,'perframe_params')
        obj.perframe_params = featureLexicon.perframe_params;
        % pf_fields = fieldnames(featureLexicon.perframe_params);
        % for ndx = 1:numel(pf_fields),
        %   obj.perframe_params.(pf_fields{ndx}) = ...
        %     featureLexicon.perframe_params.(pf_fields{ndx});
        % end
      elseif ~isstruct(obj.perframe_params)
        obj.perframe_params = struct;
      end
      if ~isfield(obj.perframe_params,'nroi')
        obj.perframe_params.nroi = 0;
      end
      if ~isfield(obj.perframe_params,'nflies_close')
        obj.perframe_params.nflies_close = [];
      end
      if ~isfield(obj.perframe_params,'fov'),
        %warning('fov not set, initializing to pi');
        obj.perframe_params.fov = pi;
      end
      if ~isfield(obj.perframe_params,'max_dnose2ell_anglerange')
        %warning('max_dnose2ell_anglerange not set, initializing to 127');
        obj.perframe_params.max_dnose2ell_anglerange = 127;
      end
      if ~isfield(obj.perframe_params,'nbodylengths_near')
        %warning('nbodylengths_near not set, initializing to 2.5');
        obj.perframe_params.nbodylengths_near = 2.5000;
      end
      if ~isfield(obj.perframe_params,'thetafil'),
        obj.perframe_params.thetafil = [0.0625 0.2500 0.3750 0.2500 0.0625];
      end
      
      obj.isST = isfield(featureLexicon,'st');

      % Update obj.allperframefns based on the new feature lexicon
      %obj.allperframefns =  fieldnames(featureLexicon.perframe);
      % assert(all(ismember(sublexiconPFNames,fieldnames(featureLexicon.perframe)))); % AL: true?
      obj.allperframefns = sublexiconPFNames;
            
      for i = 1:obj.perframe_params.nroi
        obj.allperframefns{end+1} = sprintf('dist2roi2_%02d',i);
        obj.allperframefns{end+1} = sprintf('angle2roi2_%02d',i);
      end
      
      for ii = 1:numel(obj.perframe_params.nflies_close)
        i = obj.perframe_params.nflies_close(ii);
        obj.allperframefns{end+1} = sprintf('nflies_close_%02d',i);
      end
      
      % Generate the necessary files now, so that any problems occur now.
      if obj.isST
        % none for now
      else
        generateMissingPerframeFiles = obj.perframeGenerate;
        for iExp=1:obj.nexps
          expName = obj.expnames{iExp};
          allPerframeFilesExist = obj.fileOfGivenTypesExistForGivenExps('perframedir',iExp);
          if ~allPerframeFilesExist
            % Figure out whether to try to generate the missing files
            obj.SetStatus('Some files missing for %s...',expName);
            if isempty(generateMissingPerframeFiles) && obj.isInteractive
              % Would be good to move UI stuff out of JLabelData, which is
              % essentially a model in MVC terms --ALT, Apr 30, 2013
              questionString = sprintf('Experiment %s is missing required per-frame files. Generate now?',expName);
              res = questdlg(questionString, ...
                             'Generate missing per-frame files?', ...
                             'Yes','No', ...
                             'Yes');
              if strcmpi(res,'Yes')
                generateMissingPerframeFiles = true;
                obj.perframeGenerate = true;
              elseif isempty(res) || strcmpi(res,'No')
                generateMissingPerframeFiles = false;
              end
            end
            % Generate the missing files, or generate an error  
            if ~isempty(generateMissingPerframeFiles) && generateMissingPerframeFiles
              obj.SetStatus('Generating missing files for %s...',expName);
              [success,msg] = obj.GeneratePerframeFilesExceptScoreFeatures(iExp);
              if ~success
                error('JLabelData:unableToGeneratePerframeFile',msg);
              end
            else
              obj.SetStatus('Not generating missing files for %s...',expName);
              msg = 'Unable to generate per-frame files because user declined to.';
              error('JLabelData:unableToGeneratePerframeFile',msg);
            end
          end
        end
        
        % Re-load the perframe feature signals, since the PFFs may have changed
        [success,msg] = obj.loadPerframeData(obj.expi,obj.flies);
        if ~success
          error('JLabelData:unableToLoadPerframeData',msg);
        end
      end      

      % Clear the classifier, since the feature lexicon has changed
      % This also clears the features currently in use by the classifier
      % trainer
      obj.clearClassifierProper();

      obj.needsave = true;
    end    
    
    
    % ---------------------------------------------------------------------
    function addSingleScoreFeature(obj,scoreFeature)
      % This is called by setScoreFeatures() to add a single score feature.
      % Rolling back on error is handled by setScoreFeatures.
      % We assume the score feature to be added is not already present.
      % This changes obj.scoreFeatures and obj.allperframefns, but not
      % obj.perframedata.
      
      % MERGEST OK

      % Convert the scores files into perframe files.
      pfName = scoreFeature.scorefilename; % this is the base name of the score file, which also serves as the per-frame feature name
      timeStamp = scoreFeature.ts;
      nExps = obj.nexps;
      for iExp = 1:nExps
        expName = obj.expnames{iExp};
        obj.SetStatus('Generating score-based per-frame feature file %s for %s...',pfName,expName);
        [success,msg,timeStamp,newProjectName] = obj.ScoresToPerframe(iExp, ...
          pfName, ...
          timeStamp,...
          scoreFeature.classifierfile);
        if ~success
          error('JLabelData:errorGeneratingPerframeFileFromScoreFile', ...
                sprintf('Error generating score-based per-frame file %s for %s: %s',pfName,expName,msg)); %#ok
        end
        scoreFeature.classifierfile = newProjectName;
      end
      
      scoreFeature.ts = timeStamp;
      
      % store scoreFeatures in self
      assert(~any(strcmp(pfName,obj.allperframefns)),...
        'Perframe function with name ''%s'' already exists.');
      obj.scoreFeatures = [obj.scoreFeatures;scoreFeature];
      obj.allperframefns = [obj.allperframefns;pfName];
    end

    
    % ---------------------------------------------------------------------
    function deleteSingleScoreFeature(obj,scoreFeature)
      % This is called by setScoreFeatures() to delete a single score feature.
      % Rolling back on error is handled by setScoreFeatures.
      % This changes obj.scoreFeatures and obj.allperframefns, but not
      % obj.perframedata.
      
      % MERGEST OK

      % Delete from scoreFeatures
      toBeDeleted = arrayfun(@(s)isequal(s,scoreFeature),obj.scoreFeatures);
      obj.scoreFeatures = obj.scoreFeatures(~toBeDeleted);
      
      % Delete from allperframefns
      pfName = scoreFeature.scorefilename;
      tf = strcmp(pfName,obj.allperframefns);
      assert(isequal(nnz(tf),nnz(toBeDeleted)));
      notPfName = @(string)(~isequal(string,pfName));
      obj.allperframefns = cellFilter(notPfName,obj.allperframefns);      
    end
        
    
     % ---------------------------------------------------------------------
    function tf = getPerFrameFeatureSetIsNonEmpty(self)
      % tf: nclassifiers-by-1 logical vec, true iff the current set of 
      % per-frame features in use is non-empty, i.e. contains at least one per-frame feature.
      tf = ~cellfun(@isempty,self.curperframefns);
    end   
    
    
    % ---------------------------------------------------------------------
    function setScoreFeatures(obj,varargin)
      % Update obj.scoreFeatures, preserving invariants. If an exception
      % occurs, this function will roll back the object to its original
      % state.
      
      %MERGEST OK

      % Process arguments
      if length(varargin)==1
        scoreFeaturesNew = varargin{1};
      elseif length(varargin)==3
        % collect the args into a scorefeatures struct array
        scoreFeaturesFileNameListNew = varargin{1};
        timeStampListNew = varargin{2};
        scoreFileBaseNameListNew = varargin{3};
        scoreFeaturesNew = ...
          collectScoreFeatures(scoreFeaturesFileNameListNew, ...
                               timeStampListNew, ...
                               scoreFileBaseNameListNew);
      else
        error('JLabelData:internalError', ...
              'Internal error: Wrong number of arguments to JLabelData.setScoreFeatures()');
      end
      
      scoreFeaturesOld = obj.scoreFeatures;
      featureNamesInSubdialectOld = obj.allperframefns;
      % Get some other things, in case we need to roll-back
      perframeDataOld = obj.perframedata;
      perframeUnitsOld = obj.perframeunits;

      % Use try block, so we can easily roll back if anything goes amiss
      try 
        % determine which elements of each are kept, added
        [kept,added] = ...
          setDifferencesScoreFeatures(scoreFeaturesOld,scoreFeaturesNew);
        deleted = ~kept;

        % delete each of the deleted score features
        scoreFeaturesDeleted = scoreFeaturesOld(deleted);
        nDeleted = length(scoreFeaturesDeleted);
        for i = 1:nDeleted
          obj.deleteSingleScoreFeature(scoreFeaturesDeleted(i));
        end

        % add each of the added score features
        scoreFeaturesAdded = scoreFeaturesNew(added);
        nAdded = length(scoreFeaturesAdded);
        for i = 1:nAdded
          obj.addSingleScoreFeature(scoreFeaturesAdded(i));
        end

        % Re-load the perframe feature signals, since the PFFs may have changed
        [success,msg] = obj.loadPerframeData(obj.expi,obj.flies);
        if ~success
          error('JLabelData:unableToSetScoreFeatures',msg);
        end
      catch excp
        if isequal(excp.identifier,'JLabelData:unableToSetScoreFeatures') || ...
           isequal(excp.identifier,'JLabelData:errorGeneratingPerframeFileFromScoreFile')
          % unroll changes
          obj.scoreFeatures = scoreFeaturesOld;
          obj.allperframefns = featureNamesInSubdialectOld;
          obj.perframedata = perframeDataOld;
          obj.perframeunits = perframeUnitsOld;
        end
        % We always rethrow the exception, so the caller can inform the user
        rethrow(excp);
      end
      
      obj.clearClassifierProper();
      
      % note that we now have unsaved changes
      obj.needsave = true;
    end
    
    
    function SetWindowFeatureNames(obj)
      % set windowdata.featurenames based on obj.curperframefns and
      % obj.windowfeaturescellparams.
      
      % MERGESTUPDATED

      obj.SetStatus('Setting window feature names...');
      
      FAKEPFDATA = zeros(1,101);

      nclassifiers = obj.nclassifiers;
      for iCls = 1:nclassifiers
        
        %wfparams = obj.windowfeaturesparams{iCls};
        wfcparams = obj.windowfeaturescellparams{iCls};
        curPFFs = obj.curperframefns{iCls};
        Ncpff = numel(curPFFs);
        
%         if isempty(fieldnames(wfparams))
%           feature_names = {cell(1,0)};
%           % msg = 'No features selected';
%           % obj.ClearStatus();
%           % return;
%         else        
        feature_names = cell(1,Ncpff);
        parfor j = 1:Ncpff
          fn = curPFFs{j};
          [~,tmp_feature_names_all] = ...
            ComputeWindowFeatures(FAKEPFDATA,wfcparams.(fn){:},'t0',51,'t1',51);
          feature_names{j} = cellfun(@(s) [{fn},s],tmp_feature_names_all,'UniformOutput',false);
        end
%         end

        feature_names = [feature_names{:}];
        if isempty(feature_names)
          feature_names = cell(1,0);
        end
        obj.windowdata(iCls) = WindowData.windowdataSetFeaturenames(...
          obj.windowdata(iCls),{feature_names});
      end
        
      obj.ClearStatus();

    end  % method

    
    %     % ---------------------------------------------------------------------    
%     function [success,msg] = setFeatureLexiconAndTargetSpeciesRaw(obj,featureLexicon,targetSpecies,varargin)
%       % This sets the feature lexicon to the given one, and also sets the target species.  If the
%       % featureLexicon is not one of the named ones, then either no
%       % featureLexiconName should be given, or the name should be 'custom'.
%             
%       obj.targettype=targetSpecies;  % what species the targets are
%       [success,msg] = ...
%         obj.setFeatureLexiconAndFLName(featureLexicon,varargin{:});
%     end  % method
    
    
%     % ---------------------------------------------------------------------    
%     function [success,msg] = setFeatureLexiconAndTargetSpeciesFromFLName(obj,featureLexiconName)
%       % This sets the feature lexicon to the one named by
%       % featureLexiconName
%       
% %       % Only do stuff if the new lexicon name is different than the current
% %       % one
% %       if isequal(featureLexiconName,obj.featureLexiconName)
% %         success=true;
% %         return
% %       end
%       
%       % Get the lexicon itself and the associated animal type
%       [featureLexicon,animalType]= ...
%         featureLexiconFromFeatureLexiconName(featureLexiconName);     
% 
%       % Store the lexicon-associated stuff in obj
%       %[success,msg] = obj.setFeatureLexiconAndTargetSpeciesRaw(featureLexicon,animalType,featureLexiconName);
%       [success,msg] = obj.setFeatureLexiconAndFLName(featureLexicon,featureLexiconName);
%       
%       % Set the animal type
%       obj.targettype=animalType;
%     end  % method
    
    
%     % ---------------------------------------------------------------------    
%     function [success,msg] = setFeatureLexiconAndFLName(obj,featureLexicon,featureLexiconName)
%       % This sets the feature lexicon to the given one.  If the
%       % featureLexicon is not one of the named ones, then either no
%       % featureLexiconName should be given, or the name should be 'custom'.
%       
%       % process args
%       if ~exist('featureLexiconName','var')
%         featureLexiconName='custom';
%       end
%       
%       % Setup the default return values
%       %success = false;
%       msg = '';
% 
%       % Store the lexicon-associated stuff in obj
%       obj.featureLexiconName=featureLexiconName;
%       obj.featureLexicon=featureLexicon;  % save to obj
%       %obj.targettype=targetSpecies;  % what species the targets are
%       
%       % Update obj.perframe_params based on the new feature lexicon
%       if isfield(featureLexicon,'perframe_params'),
%         obj.perframe_params=featureLexicon.perframe_params;
%         % pf_fields = fieldnames(featureLexicon.perframe_params);
%         % for ndx = 1:numel(pf_fields),
%         %   obj.perframe_params.(pf_fields{ndx}) = ...
%         %     featureLexicon.perframe_params.(pf_fields{ndx});
%         % end
%       end
% 
%       % Update obj.allperframefns based on the new feature lexicon
%       obj.allperframefns =  fieldnames(featureLexicon.perframe);
%       
%       % Re-load the perframe feature signals, since the PFFs may have changed
%       obj.loadPerframeData(obj.expi,obj.flies);
% 
%       % Clear the classifier, since the feature lexicon has changed
%       % This also clears the features currently in use by the classifier
%       % trainer
%       obj.clearClassifierProper();
%       
%       % If we got this far, all is good
%       obj.needsave=true;
%       success = true;
%     end  % method

  end
  
  
  %% Experiment/Filesystem

  methods 
        

    % ---------------------------------------------------------------------
    function [success,msg] = AddExpDir(obj, ...
                                       expDirName, ...
                                       interactivemode) %#ok
      % Add a new experiment to the GUI.  This does not change self.expi,
      % and does not do any preloading.
      
      %MERGESTOK
      
      % Set default return values
      success = false; msg = '';
      
      if isnumeric(expDirName)
        return; 
      end
      
      % interactive mode is not stored in the instance vars
      if exist('interactivemode','var')
        warning(['JLabelData.AddExpDir: Whether or not the JLabelData is interactive is now stored in the object.  ' ...
                 'Setting it for the call to AddExpDir() is ignored.']);
      end
      interactivemode = obj.isInteractive;

      % make sure directory exists
      obj.SetStatus('Checking that %s exists...',expDirName);
      if ~exist(expDirName,'file'),
        error('JLabelData:expDirDoesNotExist', '%s', ...
              expDirName);
      end
      
      if obj.fromAPT && ~strcmp(obj.aptInfo.apt_trx_type,'custom') % ~obj.aptInfo.has_trx && 
      % Maybe generate trx file if project is imported from APT and doesn't
      % have trx.
        trxFileNameAbs = fullfile(expDirName,obj.GetFileName('trx'));
        trkFileName =  fullfile(expDirName,obj.GetFileName('trk'));
        if ~exist(trxFileNameAbs,'file')          
          [trx,success1,msg1] = create_trx_apt(trkFileName,obj.aptInfo);  %#ok<PROPLC>
          if ~success1
            success = false;
            msg = msg;
            return;
          end
          A = struct('trx',trx); %#ok<PROPLC>
          save(trxFileNameAbs,'-struct','A');
        end        
      end
 
      [~,expName] = myfileparts(expDirName);
      
      % create clips dir
      clipsdir = obj.GetFileName('clipsdir');
      outclipsdir = fullfile(expDirName,clipsdir);  %#ok

      % okay, checks succeeded, start storing stuff
      obj.expdirs{end+1} = expDirName;

      % If we call remove, it's to roll-back a failed experiment add, so we
      % don't want to set needsave in this case.
      needSaveIfSuccessfulRemoval = false;
      
      % Update the status table        
      obj.SetStatus('Updating status table for %s...',expName);
      [success,msg,missingfiles] = obj.UpdateStatusTable('',obj.nexps);
      missingfiles = missingfiles{obj.nexps};
      if ~success,
        obj.SetStatus('Bad experiment directory %s...',expDirName);
        obj.RemoveExpDirs(obj.nexps,needSaveIfSuccessfulRemoval);
        return;
      end
      
      % check for existence of necessary files in this directory
      if ~obj.filesfixable,
        msg = sprintf(['Experiment %s is missing required files that cannot '...
          'be generated within this interface. Removing...'],expDirName);
        obj.SetStatus('Bad experiment directory %s...',expDirName);
        success = false;
        % undo
        obj.RemoveExpDirs(obj.nexps,needSaveIfSuccessfulRemoval);
        return;
      end
      
      % Always regenerate the score feature perframe files, to make sure
      % they're from the right classifier
      % Have to do this _before_ calling GenerateMissingFiles() so that
      % when UpdateStatusTable() gets called wihtin GenerateMissingFiles(), 
      % any score feature perframe files are already in place
      [success,msg] = obj.GenerateScoreFeaturePerframeFiles(obj.nexps);
      if ~success,
        %msg = sprintf(['Error generating missing required files %s '...
        %  'for experiment %s: %s. Removing...'],...
        %  sprintf(' %s',missingfiles{:}),expDirName,msg);
        obj.SetStatus('Error generating score feature perframe files for %s...',expName);
        obj.RemoveExpDirs(obj.nexps,needSaveIfSuccessfulRemoval);
        return
      end
      
      % If some files are missing, and we can generate them, do so
      if obj.filesfixable && ~obj.allfilesexist,
        obj.SetStatus('Some files missing for %s...',expName);
        if interactivemode && isdisplay(),
          if isempty(obj.GetGenerateMissingFiles) || ~obj.GetGenerateMissingFiles()
            if numel(missingfiles)>10,
              missingfiles = missingfiles(1:10);
              missingfiles{end+1} = ' and more ';
            end
            % Would be good to move UI stuff out of JLabelData, which is
            % essentially a model in MVC terms --ALT, Apr 30, 2013
            res = questdlg(sprintf(['Experiment %s is missing required files:%s. '...
              'Generate now?'],expDirName,sprintf(' %s',missingfiles{:})),...
              'Generate missing files?','Yes','Cancel','Yes');
            if strcmpi(res,'Yes')
              obj.SetGenerateMissingFiles();
            end
          else 
            res = 'Yes';
          end
        else
          res = 'Yes';
        end
        if strcmpi(res,'Yes'),
          obj.SetStatus('Generating missing files for %s...',expName);
          [success,msg] = obj.GenerateMissingFiles(obj.nexps);
          if ~success,
            msg = sprintf(['Error generating missing required files %s '...
              'for experiment %s: %s. Removing...'],...
              sprintf(' %s',missingfiles{:}),expDirName,msg);
            obj.SetStatus('Error generating missing files for %s...',expName);
            obj.RemoveExpDirs(obj.nexps,needSaveIfSuccessfulRemoval);
            return;
          end
        else
          obj.SetStatus('Not generating missing files for %s, not adding...',expName);
          obj.RemoveExpDirs(obj.nexps,needSaveIfSuccessfulRemoval);
          return;
        end
      end
      
%       % Convert the scores file into perframe files.      
%       for i = 1:numel(obj.scoreFeatures)
%         obj.SetStatus('Generating score-based per-frame feature file %s for %s...',obj.scoreFeatures(i).scorefilename,expName);
%         [success,msg] = obj.ScoresToPerframe(obj.nexps, ...
%                                              obj.scoreFeatures(i).scorefilename, ...
%                                              obj.scoreFeatures(i).ts);
%           if ~success,
%             obj.SetStatus('Error generating score-based per-frame file %s for %s...',obj.scoreFeatures(i).scorefilename,expName);
%             obj.RemoveExpDirs(obj.nexps,needSaveIfSuccessfulRemoval);
%             return;
%           end
%       end
      
      % Read from the trx file
      obj.SetStatus('Getting basic trx info for %s...',expName);
      trxFileNameAbs = fullfile(expDirName,obj.GetFileName('trx'));
      try
        [nFlies,firstFrames,endFrames,~,hasSex,fracSex,sex,hasperframesex] = ...
          JLabelData.readTrxInfoFromFile(trxFileNameAbs);
      catch err
         if (strcmp(err.identifier,'JAABA:JLabelData:readTrxInfoFromFile:errorReadingTrxFile'))
           msg = sprintf('Error getting basic trx info: %s',msg1);
           obj.SetStatus('Error getting basic trx info for %s, not adding...',expName);
           obj.RemoveExpDirs(obj.nexps,needSaveIfSuccessfulRemoval);
           return;
         else
           rethrow(err);
         end
      end
      obj.hassex = obj.hassex || hasSex;
      obj.hasperframesex = hasperframesex;
      obj.nflies_per_exp(end+1) = nFlies;
      obj.sex_per_exp{end+1} = sex;
      obj.frac_sex_per_exp{end+1} = fracSex;
      obj.firstframes_per_exp{end+1} = firstFrames;
      obj.endframes_per_exp{end+1} = endFrames;
      
      % Initialize the labels for the current labeling mode
      iExp = obj.nexps;
      obj.labels(iExp) = Labels.labels(1);
%       obj.labelstats(iExp).nflies_labeled = 0;
%       obj.labelstats(iExp).nbouts_labeled = 0;
      
      % Initialize the prediction data, if needed (?? --ALT, Mar 5 2013)
      if numel(obj.predictdata)<obj.nexps
        obj.SetStatus('Initializing prediction data for %s...',expName);
        obj.PredictDataInit(obj.nexps);
      end
            
      % % Set the default path to the experiment directory's parent
      obj.expdefaultpath = fileparts(expDirName);

      % Update the status
      obj.SetStatus('Successfully added experiment %s...',expDirName);
      
      % Declare victory
      obj.needsave = true;
      success = true;
    end
    
    
    function [success,msg] = AddExpDirAndLabelsFromJab(obj,jabfilename,importlabels)
      %MERGESTUPDATED
      
      assert(obj.nclassifiers==1,'Current project must be a single-classifier project.');
      % AL 20141219: In the multiclassifier world one can imagine a
      % different merge operation, where the incoming Jab contains new
      % classifiers/behaviors to be added to the current Jab.
      
      assert(~obj.isST,'ST: unsupported.'); % AL: for no good reason
      
      success = true; 
      msg = '';
      
      try
        Q = load(jabfilename,'-mat');
      catch ME
        success = false;
        msg = ME.message;
        return;
      end
      
      assert(isscalar(Q.x.classifierStuff),...
        'Jabfile ''%s'' is not a single-classifier project.',jabfilename);
      
      if obj.IsGTMode
        origExpDirNames = Q.x.gtExpDirNames;
        origLabels = Q.x.gtLabels;
      else
        origExpDirNames = Q.x.expDirNames;
        origLabels = Q.x.labels;
      end
      origLabels = Labels.modernizeLabels(origLabels,Q.x.behaviors.names);
      assert(numel(origExpDirNames)==numel(origLabels));
      
      for ndx = 1:numel(origExpDirNames)
        expdirname = origExpDirNames{ndx};
        labels = origLabels(ndx);
        if importlabels % Change the names
          assert(numel(Q.x.behaviors.names)==2);
          origBehName = Q.x.behaviors.names{1};
          origNoBehName = Q.x.behaviors.names{2};
          assert(strcmpi(origNoBehName,'none'));
          newBehName = obj.labelnames{1};
          newNoBehName = obj.labelnames{2};
          assert(strcmpi(newNoBehName,'none'));
          labels = Labels.renameBehavior(labels,origBehName,newBehName,origNoBehName,newNoBehName);
        end
        
        curexp = find(strcmp(expdirname,obj.expdirs));
        if isempty(curexp)
          [success,msg] = obj.AddExpDir(expdirname);
          if ~success
            return;
          end
          if importlabels
            obj.labels(end) = labels;
          end
        elseif importlabels
          tfLabelsUpdated = false;
          if obj.haslabels(curexp)
            dlgstr = sprintf('Experiment %s already has labels. Discard the existing (current) labels and load the new ones?',...
              expdirname);
            res = questdlg(dlgstr,'Load new labels?',...
              'Keep Existing','Load New','Cancel','Keep Existing');
            if strcmp(res,'Load New')
              obj.labels(curexp) = labels;
              tfLabelsUpdated = true;
            end
          else
            obj.labels(curexp) = labels;
            tfLabelsUpdated = true;
          end
          
          if tfLabelsUpdated && curexp==obj.expi
            % initialize .labelidx from .labels
            % AL 20141219: added tfLabelsUpdated condition; seems like we
            % only want to update .labelidx (especially from .labels) if
            % the new labels were accepted
            
            labelShort = Labels.labelsShort;
            [labelsShort,tffly] = Labels.labelsShortInit(labelShort,labels,obj.flies);            
            if tffly
              T0 = max(obj.GetTrxFirstFrame(obj.expi,obj.flies));
              T1 = min(obj.GetTrxEndFrame(obj.expi,obj.flies));
              labelidx = Labels.labelIdx(obj.labelnames,T0,T1);
              labelidx = Labels.labelIdxInit(labelidx,labelsShort);
              obj.labelidx = labelidx;
            else
              % none. obj.flies not present in new labels. .labelidx not 
              % updated, so may contain bouts that differ from .labels
            end
          end          
        end
      end
    end    


    % ---------------------------------------------------------------------
    function [success,msg] = RemoveExpDirs(obj,expi,needSaveIfSuccessful)
      % [success,msg] = RemoveExpDirs(obj,expi,[needSaveIfSuccessful])
      % Removes experiments in expi from the GUI. If the currently loaded
      % experiment is removed, then a different experiment may be preloaded. 
      % If needSaveIfSuccessful is defined and true, the JLabelData's
      % needsave property will be set to true if the experiment is
      % successfully removed.  If needSaveIfSuccessful is not defined, is
      % empty, or is false, needsave will not be changed.  It is useful to
      % set needSaveIfSuccessful to false if RemoveExpDirs() is being called
      % in order to roll-back the partially-completed addition of an
      % experiment.

      %MERGESTUPDATED
      
      success = false;
      msg = '';
      
      if ~exist('needSaveIfSuccessful','var') || isempty(needSaveIfSuccessful)
        needSaveIfSuccessful = true;
      end
      
      if any(obj.nexps < expi) || any(expi < 1),
        msg = sprintf('expi = %s must be in the range 1 < expi < nexps = %d',mat2str(expi),obj.nexps);
        return;
      end
      
      origNExp = obj.nexps;
      newExpNumbers = zeros(1,obj.nexps);
      for ndx = 1:obj.nexps
        if ismember(ndx,expi);
          newExpNumbers(1,ndx) = 0;
        else
          newExpNumbers(1,ndx) = ndx-nnz(expi<ndx);
        end
      end
      
      if ~(numel(obj.expdirs)<expi)
        obj.expdirs(expi) = [];   
      end
      
      assert(nnz(newExpNumbers>0)==numel(obj.expdirs));
      
      if ~(numel(obj.nflies_per_exp)<expi); obj.nflies_per_exp(expi) = []; end
      if ~(numel(obj.sex_per_exp)<expi); obj.sex_per_exp(expi) = []; end
      if ~(numel(obj.frac_sex_per_exp)<expi); obj.frac_sex_per_exp(expi) = []; end
      if ~(numel(obj.firstframes_per_exp)<expi); obj.firstframes_per_exp(expi) = []; end
      if ~(numel(obj.endframes_per_exp)<expi); obj.endframes_per_exp(expi) = []; end
      if ~(numel(obj.labels)<expi); obj.labels(expi) = []; end
%       if ~(numel(obj.labelstats)<expi); obj.labelstats(expi) = []; end

      % Clean window data
      for iCls = 1:numel(obj.windowdata)
        wdExp = obj.windowdata(iCls).exp;
        assert(iscolumn(wdExp) || isequal(wdExp,[]));
        newExp = newExpNumbers(wdExp);
        obj.windowdata(iCls).exp = newExp(:);
      end
      % removed exps now have obj.windowdata.exp==0
      obj.windowdata = WindowData.windowdataTrim(obj.windowdata,@(x)x.exp==0);
      % AL 20141121. There used to be some odd code here:
      % - Not all windowdata fields trimmed, eg labelidx_old
      % - Code implied that many fields (eg labelidx_cur) did not have to 
      % size consistent with windowdata.X, eg size of labelidx_cur implied
      % to be different than labelidx_new
      WindowData.windowdataVerify(obj.windowdata);      

      for curex = sort(expi(:)','descend') %#ok<UDIM>
        if numel(obj.predictdata)>=curex
          obj.predictdata(curex) = [];
        end
      end
      
      nCls = obj.nclassifiers;
      pbfnames = fieldnames(obj.predictblocks);
      for iCls = 1:nCls
        if ~isempty(obj.predictblocks(iCls).expi)
          idxcurr = ismember(obj.predictblocks(iCls).expi,expi);
          for fndx = 1:numel(pbfnames)
            obj.predictblocks(iCls).(pbfnames{fndx})(idxcurr) = [];
            % all fields of predictblocks are currently row vecs
          end
          tmp = newExpNumbers(obj.predictblocks(iCls).expi);
          obj.predictblocks(iCls).expi = tmp(:)';
        end
      end
      
      if obj.gtMode,
        if ~isempty(obj.randomGTSuggestions)
          dirstoremove = expi;
          dirstoremove(dirstoremove > numel(obj.randomGTSuggestions))=[];
          obj.randomGTSuggestions(dirstoremove) = [];
        end
        
        if ~isempty(obj.loadedGTSuggestions)
          dirstoremove = expi;
          dirstoremove(dirstoremove > numel(obj.loadedGTSuggestions))=[];
          obj.loadedGTSuggestions(dirstoremove) = [];
        end
        
        if ~isempty(obj.balancedGTSuggestions)
          expidx = [obj.balancedGTSuggestions.exp];
          toremove = ismember(expidx,expi);
          obj.balancedGTSuggestions(toremove) = [];
          for ndx = 1:numel(obj.balancedGTSuggestions)
            obj.balancedGTSuggestions(ndx).exp = newExpNumbers(obj.balancedGTSuggestions(ndx).exp);
          end
        end
      end
      % update current exp, flies
      if ~isempty(obj.expi) && obj.expi > 0 
        if ismember(obj.expi,expi), % The current experiment was removed.
          newexpi = find(newExpNumbers(obj.expi+1:end),1);
          if isempty(newexpi), % No next experiment.
            newexpi = find(newExpNumbers(1:obj.expi-1),1,'last');
            if isempty(newexpi), % No previous experiment either.
              newexpi = 0;
            else
              newexpi = newExpNumbers(newexpi);
            end
          else
            newexpi = newExpNumbers(obj.expi+newexpi);
          end
          obj.expi = 0;
          obj.flies = nan(size(obj.flies));
          if obj.nexps > 0,
            obj.setCurrentTarget(newexpi,1);
          end
        else
          obj.expi = obj.expi - nnz(ismember(1:obj.expi,expi));
        end
      end
      
      % Set needsave, if called for
      if needSaveIfSuccessful
        obj.needsave = true;
      end
      
      % Declare victory
      success = true;
    end  % method

    
    % ---------------------------------------------------------------------
    function [success,msg] = SetExpDirs(obj,expdirs)
    % Changes what experiments are currently being used for this
    % classifier. This function calls RemoveExpDirs to remove all current
    % experiments not in expdirs, then calls AddExpDirs to add the new
    % experiment directories. 

      % MERGESTUPDATED
      
      success = false;
      msg = '';
      
      if isnumeric(expdirs), % AL ?
        return;
      end
      
      if nargin < 2,
        error('Usage: obj.SetExpDirs(expdirs)');
      end
                  
      oldexpdirs = obj.expdirs;
      
      % remove all oldexpdirs that aren't in expdirs
      expi = find(~ismember(oldexpdirs,expdirs));
      if ~isempty(expi)
        [success1,msg] = obj.RemoveExpDirs(expi);
        if ~success1,
          return;
        end
      end

      % add all new expdirs that weren't in oldexpdirs
      idx = find(~ismember(expdirs,oldexpdirs));
      success = true;
      for i = idx(:)'
        [success1,msg1] = obj.AddExpDir(expdirs{i});
        success = success && success1;
        if isempty(msg),
          msg = msg1;
        else
          msg = sprintf('%s\n%s',msg,msg1);
        end        
      end
    end

        
    % ---------------------------------------------------------------------
    function res = GetFileName(obj,fileType)
    % res = GetFileName(obj,fileType)
    % Get base name of file of the input type fileType.
      switch fileType,
        case 'movie',
          res = obj.moviefilename;
        case 'movieindex'
          res = obj.movieindexfilename;
        case 'trx',
          res = obj.trxfilename;
        case {'perframedir','perframe'},
          res = obj.perframedir;
        case {'clipsdir','clips'},
          if ischar(obj.clipsdir)
            res = strtrim(obj.clipsdir);
          else
            res = 'clips';
          end            
        case 'scores',
          res = obj.scorefilename;
        case 'stfeatures',
          res = obj.stfeatures;
        case 'trk'
          res = obj.trkfilename;
        otherwise
          error('Unknown file type %s',fileType);
      end
    end

    
    % ---------------------------------------------------------------------    
    function [filename,timestamp,tffound] = GetFile(obj,fileType,expi)
    % filename: full path(s). Either char or cellstr
    % timestamp: timestamps for files that exist. Size corresponding to filename.
    % tffound: logical vector. Size corresponding to filename.
    %
    % If tffound is true, filename is the full filename and timestamp is the filesystem timestamp.
    % If tffound is false, filename is the file-that-was-tried and timestamp=-inf.
    %
    % Vector output occurs eg for fileType='score' when using multiple
    % classifiers.
  
    %ST OK
    
      % base name
      fileNameLocal = obj.GetFileName(fileType);
      
      % if this is an output file, only look in output experiment directory
%       if dowrite && JLabelData.IsOutputFile(file),
%         %expdirs_try = obj.outexpdirs(expi);
%         expdirs_try = obj.expdirs(expi);
%       else
%         % otherwise, first look in output directory, then look in input
%         % directory
%         %expdirs_try = {obj.outexpdirs{expi},obj.expdirs{expi}};
%         expdirs_try = obj.expdirs(expi);
%       end
      expdirs_try = obj.expdirs(expi);
      
      % non-char or empty filename signals its absence
      if ~ischar(fileNameLocal) && ~iscellstr(fileNameLocal) || isempty(fileNameLocal)
        filename = fileNameLocal;
        timestamp = -inf;
        tffound = false;
        return
      end
      
      if ischar(fileNameLocal)
        [tffound,filename,timestamp] = JLabelData.GetFileRaw(fileType,expdirs_try,fileNameLocal);
      elseif iscellstr(fileNameLocal)
        [tffound,filename,timestamp] = cellfun(@(x)JLabelData.GetFileRaw(fileType,expdirs_try,x),fileNameLocal,'uni',0);
        timestamp = cell2mat(timestamp);
        tffound = cell2mat(tffound);
      end
    end    
  
    
    % ---------------------------------------------------------------------
    function SetGenerateMissingFiles(obj)
      obj.perframeGenerate = true;
    end
    
    
    % ---------------------------------------------------------------------
    function perframeGenerate = GetGenerateMissingFiles(obj)
      perframeGenerate = obj.perframeGenerate;
    end

    
    % ---------------------------------------------------------------------
    function [success,msg] = GenerateMissingFiles(obj,expi)
    % [success,msg] = GenerateMissingFiles(obj,expi)
    % Generate required, missing files for experiments expi. 
    % TODO: implement this!
    
    % ST OK
      
      success = true;
      msg = '';
            
      for i = 1:numel(obj.filetypes),
        file = obj.filetypes{i};
        if obj.IsRequiredFile(file) && obj.CanGenerateFile(file) && ...
            ~obj.FileExists(file,expi),
          fprintf('Generating %s for %s...\n',file,obj.expnames{expi});
          switch file
            case 'perframedir'
              [success1,msg1] = obj.GeneratePerframeFilesExceptScoreFeatures(expi);
              success = success && success1;
              if ~success1,
                msg = [msg,'\n',msg1]; %#ok<AGROW>
              end
%               [success2,msg2] = obj.GenerateScoreFeaturePerframeFiles(expi);
%               success = success && success2;
%               if ~success2,
%                 msg = [msg,'\n',msg2]; %#ok<AGROW>
%               end
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

    
    % ---------------------------------------------------------------------
    function RemoveArenaPFs(obj)
      %Update .allperframefns, .curperframefns, .windowfeaturesparams,
      %.windowfeaturescellparams
      
      %MERGESTUPDATED
      
      % Determine which PFFs to remove
      settings = obj.featureLexicon;
      apff = obj.allperframefns;      
      Napff = numel(apff);
      tfRm = false(Napff,1);
      for i = 1:Napff
        curpf = apff{i};
        if any(strcmp(curpf,{obj.scoreFeatures(:).scorefilename}))
          continue;
        end
        curtypes = settings.perframe.(curpf).type;
        tfRm(i) = any(strcmpi(curtypes,'arena')) || any(strcmpi(curtypes,'position'));
      end
      
      % Remove
      apffRm = apff(tfRm);
      obj.allperframefns(tfRm,:) = [];
      nCls = obj.nclassifiers;
      for iCls = 1:nCls
        fldsRm = intersect(fieldnames(obj.windowfeaturesparams{iCls}),apffRm);
        obj.windowfeaturesparams{iCls} = rmfield(obj.windowfeaturesparams{iCls},fldsRm);
        % AL 20141215: next line is legacy, don't understand why it is 
        % necessary given we are just removing fields
        obj.windowfeaturesparams{iCls} = JLabelData.convertTransTypes2Cell(obj.windowfeaturesparams{iCls});
        obj.windowfeaturescellparams{iCls} = JLabelData.convertParams2CellParams(obj.windowfeaturesparams{iCls});
        
        tf = ismember(obj.curperframefns{iCls},apffRm);
        obj.curperframefns{iCls}(tf,:) = [];
      end      
    end
    
    
    % ---------------------------------------------------------------------
    function [success,msg] = GeneratePerframeFilesExceptScoreFeatures(obj,expi)
      success = false; %#ok<NASGU>
      msg = '';

      perframedir = obj.GetFile('perframedir',expi);
      
      isInteractive=obj.isInteractive;
      if ~isInteractive,
        dooverwrite = false;
      elseif ~isempty(obj.perframeOverwrite) 
        dooverwrite = obj.perframeOverwrite;
      elseif exist(perframedir,'dir'),
        res = questdlg('Do you want to leave existing files alone or regenerate them?',...
                       'Regenerate existing files?', ...
                       'Leave Alone','Regenerate', ...
                       'Leave Alone');
        dooverwrite = strcmpi(res,'Regenerate');
        obj.perframeOverwrite = dooverwrite;
      else
        dooverwrite = true;
      end
      
      expdir = obj.expdirs{expi};
      
      hwait=-1;  % make sure hwait var exists, and is not a graphics handle
                 % unless it gets overwritten
      if isInteractive
        hwait = mywaitbar(0,sprintf('Generating perframe files for %s',expdir),'interpreter','none');
      else
        fprintf('Generating perframe files for %s\n',expdir);
      end
      
      perframetrx = Trx('trxfilestr',obj.GetFileName('trx'),...
                        'moviefilestr',obj.GetFileName('movie'),...
                        'perframedir',obj.GetFileName('perframedir'),...
                        'default_landmark_params',obj.landmark_params,...
                        'perframe_params',obj.perframe_params,...
                        'trkfilestr',obj.GetFileName('trk'),...
                        'stFeatures',obj.stFeatures,...
                        'stInfo',obj.stInfo);
                        %'rootwritedir',expdir);
%                        'rootwritedir',obj.rootoutputdir);
      
      perframetrx.AddExpDir(expdir,'dooverwrite',dooverwrite,'openmovie',false,'aptInfo',obj.aptInfo,'fromAPT',obj.fromAPT);
      
      if isempty(fieldnames(obj.landmark_params)) && ~perframetrx.HasLandmarkParams && obj.arenawarn,
        if expi>1,
          success = false;
          msg = ['Landmark params were not defined in the configuration file or in the trx file for the current experiment. '...
              'Cannot compute arena perframe features for the experiment and removing it'];
          if isInteractive && ishandle(hwait),
              delete(hwait);
          end
          return;
        end

        if isInteractive,
          uiwait(warndlg(['Landmark params were not defined in the configuration file'...
            ' or in the trx file. Not computing arena features and removing them from the perframe list']));
        else
          fprintf('Landmark params were not defined in the configuration file. Not computing arena features and removing them from the perframe list.\n');
        end
        obj.RemoveArenaPFs();
        obj.arenawarn = false;
      end
      
      perframefiles = obj.GetPerframeFiles(expi);
      for i = 1:numel(obj.allperframefns),
        fn = obj.allperframefns{i};
        %ndx = find(strcmp(fn,obj.allperframefns));
        file = perframefiles{i};
        if ~dooverwrite && exist(file,'file'),
          continue;
        end        
        if isInteractive
          hwait = mywaitbar(i/numel(obj.allperframefns),hwait,...
            sprintf('Computing %s and saving to file %s',fn,file));
        else
          fprintf('Computing %s and saving to file %s\n',fn,file);
        end
        
        % Generate the per-frame files that are not score features
        % Don't generate the per-frame files from scores here anymore..
        
%          try
           if isempty(obj.scoreFeatures) || ~any(strcmp(fn,{obj.scoreFeatures(:).scorefilename}))
             fprintf('%s...\n',fn);
            perframetrx.(fn);
           end        
%          catch excp
%            if isInteractive && ishandle(hwait),
%              delete(hwait);
%            end
%            if isequal(excp.identifier,'MATLAB:UndefinedFunction')
%              msg=sprintf('Unable to calculate per-frame feature %s.',fn);
%              success=false;
%              return
%           else
%              rethrow(excp);
%            end
%          end
      end
      
      if isInteractive && ishandle(hwait),
        delete(hwait);
      end
      
      success = true;
      
    end  % method

    
    % ---------------------------------------------------------------------
    function [success,msg] = GenerateScoreFeaturePerframeFiles(obj,expi)
      % MERGEST OK
      
      success = true;
      msg = '';
      % Convert the scores file into perframe files.
      for i = 1:numel(obj.scoreFeatures)
        %obj.SetStatus('Generating score-based per-frame feature file %s for %s...',obj.scoreFeatures(i).scorefilename,expName);
        [success,msg] = obj.ScoresToPerframe(expi, ...
          obj.scoreFeatures(i).scorefilename, ...
          obj.scoreFeatures(i).ts);
        if ~success,
          %obj.SetStatus('Error generating score-based per-frame file %s for %s...',obj.scoreFeatures(i).scorefilename,expName);
          %obj.RemoveExpDirs(obj.nexps,needSaveIfSuccessfulRemoval);
          return
        end
      end
    end

    
    % ---------------------------------------------------------------------
    function [success,msg,ts,projectName] = ...
        ScoresToPerframe(obj,expi,fileName,ts,projectName)
      % Create perframe file from scorefile in toplevel-dir
      %
      % ts: timestamp
      
      % MERGESTUPDATED
      
      persistent perframescoresfile_didselectyes;
      
      success = true;
      msg = '';

      outdir = obj.expdirs{expi};
      scoresFileIn = [fullfile(outdir,fileName) '.mat'];
      scoresFileOut = [fullfile(outdir,obj.GetFileName('perframe'),fileName) '.mat'];
      if ~exist(scoresFileIn,'file')
        success = false; 
        msg = sprintf('Scores file %s does not exist to be used as perframe feature',scoresFileIn);
        return;
      end
      
      Q = load(scoresFileIn);
      if obj.isInteractive
        if Q.timestamp > ts
          res = questdlg(sprintf('The timestamp for scores file %s (%s) is newer than that expected for this project (%s). Do you want to update the current project, proceed without updating, or cancel?',...
              scoresFileIn,datestr(Q.timestamp),datestr(ts)),'Scores file timestamp mismatch','Update project','Proceed without updating','Cancel','Update project');
          if strcmpi(res,'Update project')
            
            while true
              [f,p] = uigetfile('*.jab',sprintf('Choose the project corresponding to %s',fileName),projectName);
              if ~ischar(f)
                break;
              end
              newProjectName = fullfile(p,f);
              if ~exist(newProjectName,'file'),
                uiwait(warndlg(sprintf('File %s does not exist',newProjectName)));
                continue;
              end
              projectData = loadAnonymous(newProjectName);
              projectData = Macguffin(projectData);
              projectData.modernize();
              assert(isscalar(projectData.classifierStuff),...
                'Project ''%s'' is a multi-classifier project; this is currently unsupported.',newProjectName);              
              if projectData.classifierStuff.timeStamp ~= Q.timestamp
                uiwait(warndlg(sprintf('Timestamp in project file %s = %s does not match time stamp in score file %s = %s',newProjectName,datestr(projectData.classifierStuff.timeStamp),fileName,datestr(ts))));
                continue;
              end
              if ~strcmp(projectData.file.scorefilename{1},[fileName,'.mat']),
                uiwait(warndlg(sprintf('Score file name in %s = %s does not match score file name %s',newProjectName,projectData.file.scorefilename,[fileName,'.mat'])));
                continue;
              end
              projectName = newProjectName;
              ts = projectData.classifierStuff.timeStamp;
              break;
            end
            
          elseif strcmpi(res,'Proceed without updating')
            Q.timestamp = ts;
          end
        elseif Q.timestamp < ts,
          res = questdlg(sprintf('The timestamp for scores file %s (%s) is older than that expected for this project (%s). Do you want to add this experiment anyways or cancel?',...
            scoresFileIn,datestr(Q.timestamp),datestr(ts)),'Scores file timestamp mismatch','Add experiment anyway','Cancel','Cancel');
          if strcmpi(res,'Add experiment anyway'),
            Q.timestamp = ts;
          end
        end
      end
      if Q.timestamp ~= ts % check the timestamps match the classifier's timestamp.
        success = false; 
        msg = sprintf(['The scores file %s was generated using a classifier' ...
          ' that was saved on %s while the classifier chosen was saved on %s'],...
          scoresFileIn,datestr(Q.timestamp),datestr(ts));
      end
      
      OUT = struct();
      OUT.units = struct();
      OUT.units.num = {'scores'};
      OUT.units.den = {''};
      % TODO: Allen, figure out how to handle this correctly in the
      % multi-class case
      if iscell(Q.allScores),
        if numel(Q.allScores) > 1,
          error('Not implemented: do not how to do this with multiple classifiers.');
        else
          Q.allScores = Q.allScores{1};
        end
      end
      for ndx = 1:numel(Q.allScores.scores)
        t0 = Q.allScores.tStart(ndx);
        t1 = Q.allScores.tEnd(ndx);
        OUT.data{ndx} = Q.allScores.scores{ndx}(t0:t1);
      end
      try
        save(scoresFileOut,'-struct','OUT');
      catch ME,
        questmsg = sprintf('Could not write perframe file from scores file: %s.',fileName);
        if obj.isInteractive && isempty(perframescoresfile_didselectyes)
          button = questdlg([questmsg,' Continue?'],'Continue','Yes');
          if ~strcmp(button,'Yes')
            success = false;
            msg = ME.message;
          else
            perframescoresfile_didselectyes = true;
          end
        else
          warning(questmsg); %#ok<SPWRN>
        end
      end
    end

    
    % ---------------------------------------------------------------------
    function [windowfeaturesparams,windowfeaturescellparams] = GetPerframeParams(obj)
      windowfeaturesparams = obj.windowfeaturesparams; 
      windowfeaturescellparams = obj.windowfeaturescellparams; 
    end  
    
    
    % ---------------------------------------------------------------------
    function [filenames,timestamps] = GetPerframeFiles(obj,expi)
    % [filenames,timestamps] = GetPerFrameFiles(obj,file,expi)
    % Get the full path to the per-frame mat files for experiment expi
      
%       if nargin < 3,
%         dowrite = false;
%       end
      
      fn = obj.GetFileName('perframedir');
      
      % if this is an output file, only look in output experiment directory
%       if dowrite && JLabelData.IsOutputFile('perframedir'),
%         %expdirs_try = obj.outexpdirs(expi);
%         expdirs_try = obj.expdirs(expi);
%       else
%         % otherwise, first look in output directory, then look in input
%         % directory
%         %expdirs_try = {obj.outexpdirs{expi},obj.expdirs{expi}};
%         expdirs_try = obj.expdirs(expi);
%       end
      expdirs_try = obj.expdirs(expi);
      
      filenames = cell(1,numel(obj.allperframefns));
      timestamps = -inf(1,numel(obj.allperframefns));
      
      is_pc = ispc; % repeatedly calling ispc is slow
      
      fsep = filesep;
      pffns = obj.allperframefns;
      nago = nargout;
      for i = 1:numel(obj.allperframefns),

        % loop through directories to look in
        expdir = expdirs_try{1};
        %           perframedir = fullfile(expdir,fn);  % BJA: slow by 10x
        perframedir = [expdir fsep fn];
        if is_pc && ~exist(perframedir,'dir'),
          [actualperframedir,didfind] = GetPCShortcutFileActualPath(perframedir);
          if didfind,
            perframedir = actualperframedir;
          end
        end
        %           filename = fullfile(perframedir,[obj.allperframefns{i},'.mat']);  % BJA: slow by 10x
        filename = [perframedir fsep pffns{i} '.mat'];
        if is_pc && ~exist(filename,'file'),
          [actualfilename,didfind] = GetPCShortcutFileActualPath(filename);
          if didfind,
            filename = actualfilename;
          end
        end
        
        filenames{i} = filename;
        if nago > 1 && exist(filename,'file');
          tmp = dir(filename);
          timestamps(i) = tmp.datenum;
        end
      end
    end
    
    
    % ---------------------------------------------------------------------
    function [fe,ft] = FileExists(obj,fileType,expi)
    % [fe,ft] = FileExists(obj,file,expi)
    % Returns whether the input file exists for the input experiment. 
      filei = find(strcmpi(fileType,obj.filetypes),1);
      if isempty(filei),
        error('file type %s does not match any known file type',fileType);
      end
      if nargin < 3,
        expi = 1:obj.nexps;
      end
      fe = obj.fileexists(expi,filei);
      ft = obj.filetimestamps(expi,filei);
    end  % method
    
    % ---------------------------------------------------------------------
    function varargout = get_readframe_fcn(obj,expi)
    % [readframe,nframes,fid,headerinfo] = get_readframe_fcn(obj,expi)
    % Calls get_readframe_fcn with the appropriate arguments for experiment
    % expi. This function does not do any checks for existence of video,
    % etc. 
    
      moviefilename = obj.GetFile('movie',expi);
      if ischar(moviefilename)
        [~,~,ext] = fileparts(moviefilename);
      else
        [~,~,ext] = fileparts(moviefilename{1});
      end
      ismaybeindexedfile = ismember(lower(ext),{'.mjpg','.mjpeg'});
      if ismaybeindexedfile,
        movieindexfilename = obj.GetFile('movieindex',expi);
      end
      
      varargout = cell(1,nargout);
      if ismaybeindexedfile && ischar(movieindexfilename),
        [varargout{:}] = get_readframe_fcn(moviefilename,'indexfilename',movieindexfilename);
      else
        [varargout{:}] = get_readframe_fcn(moviefilename);
      end
      

    end % method
    
    
    
    % ---------------------------------------------------------------------
    function [successes,msg] = CheckMovies(obj,expis)
    % [successes,msg] = CheckMovies(obj,expis)
    % check that the movie files exist and can be read for the input
    % experiments.
      
      successes = []; msg = '';
      
      if nargin < 2,
        expis = 1:obj.nexps;
      end
      
      if isempty(expis),
        return;
      end
      
      successes = true(1,numel(expis));
      
      if ~obj.ismovie,
        return;
      end
      
      for i = 1:numel(expis),
        all_moviefilenames = obj.GetFile('movie',expis(i));
        all_movieindexfilenames = obj.GetFile('movieindex',expis(i));
        obj.SetStatus('Checking movie %s...',moviefilename);
        
        if ischar(all_moviefilenames),
          all_moviefilenames = {all_moviefilenames};
          all_movieindexfilenames = {all_movieindexfilenames};
        end
        
        for ndx = 1:numel(all_moviefilenames)
          moviefilename = all_moviefilenames{ndx};
          movieindexfilename = all_movieindexfilenames{ndx};
          % check for file existence
          if ~exist(moviefilename,'file'),
            successes(i) = false;
            msg1 = sprintf('File %s missing',moviefilename);
            if isempty(msg),
              msg = msg1;
            else
              msg = sprintf('%s\n%s',msg,msg1);
            end
          else

            if ischar(movieindexfilename) && ~exist(movieindexfilename,'file'),
              successes(i) = false;
              msg1 = sprintf('File %s missing',movieindexfilename);
              if isempty(msg),
                msg = msg1;
              else
                msg = sprintf('%s\n%s',msg,msg1);
              end

            else

              % try reading a frame
              %           try
              [readframe,~,movie_fid] = ...
                obj.get_readframe_fcn(expis(i));
              if movie_fid <= 0,
                error('Could not open movie %s for reading',moviefilename);
              end
              readframe(1);
              fclose(movie_fid);
              %           catch ME,
              %             successes(i) = false;
              %             msg1 = sprintf('Could not parse movie %s: %s',moviefilename,getReport(ME));
              %             if isempty(msg),
              %               msg = msg1;
              %             else
              %               msg = sprintf('%s\n%s',msg,msg1);
              %             end
              %           end

            end
          end
        end
      end
      
      obj.ClearStatus();
      
    end  % method
    
    
    % ---------------------------------------------------------------------
    function [success,msg,missingfiles] = UpdateStatusTable(obj,filetypes,expis)
    % [success,msg] = UpdateStatusTable(obj,filetypes,expis)
    % Update the tables of what files exist for what experiments. This
    % returns false if all files were in existence or could be generated
    % and now they are/can not.  It mutates the instance variables
    % fileexists and filetimestamps, and sets filesfixable and
    % allfilesexist.
    %
    % missingfiles: Cell vector of length obj.nexps. Each element is a 
    % cellstr. QUIRKY API: Only missingfiles(expis) are meaningful, 
    % containing missing files for expis. Remaining elements are
    % empty/indeterminate.

      % Initialize return vars
      msg = '';
      success = false;
      missingfiles = cell(1,obj.nexps);
      for i = 1:obj.nexps,
        missingfiles{i} = {};
      end

      % Process arguments
      if ~exist('filetypes','var') || isempty(filetypes)
        filetypes = obj.filetypes;
      end
      if ~exist('expis','var') || isempty(expis)
        expis = 1:obj.nexps;
      end

      % Check that the file types requested are valid
      [fileTypeIsValid,fileTypeIndices] = ismember(filetypes,obj.filetypes);
      if any(~fileTypeIsValid),
        msg = 'Unknown filetypes';
        return
      end

      % Get the file existance and time-stamp tables
      [fileExists,fileTimeStamps,missingFileNames] = ...
        obj.fileOfGivenTypesExistForGivenExps(filetypes,expis);
      
      % Update the state vars accordingly
      obj.fileexists(expis,fileTypeIndices) = fileExists;
      obj.filetimestamps(expis,fileTypeIndices) = fileTimeStamps;
      missingfiles(expis) = missingFileNames;
      
      % store old values to see if latest change broke something
      old_filesfixable = obj.filesfixable;
      old_allfilesexist = obj.allfilesexist;

      % initialize summaries to true
      [allRequiredFilesExist, ...
       missingFilesCanBeGenerated, ...
       oneMissingFileTypeThatCantBeGenerated] = ...
        obj.allRequiredFilesExist(obj.fileexists);
      
      % Write the file status summaries into the object state, and
      % modify the error message if necessary to reflect missing files that
      % cannot be generated.
      obj.allfilesexist = allRequiredFilesExist;
      if allRequiredFilesExist
        obj.filesfixable = true;
      else
        obj.filesfixable = missingFilesCanBeGenerated;
        if ~missingFilesCanBeGenerated
          msg1 = sprintf('%s missing and cannot be generated.',oneMissingFileTypeThatCantBeGenerated);
          if isempty(msg),
            msg = msg1;
          else
            msg = sprintf('%s\n%s',msg,msg1);
          end
        end
      end

      % Modify the error message to reflect missing files for each exp
      nFileNamesMissingTotal = sum(cellfun(@(c)length(c),missingFileNames));
      if nFileNamesMissingTotal>0
        stringForEachExpWithMissingFiles={};        
        for i=1:length(missingFileNames)
          expi=expis(i);
          missingFileNamesThis=missingFileNames{i};
          if ~isempty(missingFileNamesThis)
            stringForEachExpWithMissingFiles{end+1}= ...
              [civilizedStringFromCellArrayOfStrings(missingFileNamesThis) ' from experiment ' obj.expnames{expi}];  %#ok
          end
        end
        msg= [msg ' Missing ' civilizedStringFromCellArrayOfStrings(stringForEachExpWithMissingFiles,';')];
      end
      
      % fail if was ok and now not ok
      success = ~(old_allfilesexist || old_filesfixable) || ...
                 (obj.allfilesexist || obj.filesfixable);      
    end  % method
    
    
    % ---------------------------------------------------------------------
    function [success,msg] = SetTrxFileName(obj,trxfilename)
      % [success,msg] = SetTrxFileName(obj,trxfilename)
      % set the name of the trx file within the experiment directory. this
      % does not currently check for missing/bad trx files, or replace
      % preloaded trx data, so you really shouldn't call it if expdirs are
      % loaded. (TODO)

      success = false;
      msg = '';
      if ischar(trxfilename),
        if ischar(obj.trxfilename) && strcmp(trxfilename,obj.trxfilename),
          success = true;
          return;
        end
        obj.trxfilename = trxfilename;
        [success,msg] = obj.UpdateStatusTable('trx');        
        % TODO: check that trx are parsable, remove bad experiments, update
        % preloaded trx
      end
    end  % method
    
    
    % ---------------------------------------------------------------------
    function [success,msg] = SetMovieFileName(obj,moviefilename,movieindexfilename)
    % change/set the name of the movie within the experiment directory
    % will fail if movie files don't exist for any of the current
    % experiment directories (checked by CheckMovies)

      if nargin < 3,
        movieindexfilename = 0;
      end
    
      success = false; msg = '';
      
      indexfilesmatch = (ischar(movieindexfilename) && ischar(obj.movieindexfilename) && ...
          strcmp(movieindexfilename,obj.movieindexfilename)) || ...
          (~ischar(movieindexfilename) && ~ischar(obj.movieindexfilename));

      if ischar(moviefilename),
        if indexfilesmatch && ischar(obj.moviefilename) && strcmp(moviefilename,obj.moviefilename),
          success = true;
          return;
        end
        oldmovieindexfilename = obj.movieindexfilename;
        oldmoviefilename = obj.moviefilename;
        obj.moviefilename = moviefilename;
        obj.movieindexfilename = movieindexfilename;
        %obj.ismovie = ~isempty(moviefilename) && obj.openmovie;
        [success1,msg] = obj.CheckMovies();
        if ~success1,
          obj.moviefilename = oldmoviefilename;
          obj.movieindexfilename = oldmovieindexfilename;
          return;
        end
        [success,msg] = obj.UpdateStatusTable('movie');
      elseif iscell(moviefilename)
        if indexfilesmatch && ischar(obj.moviefilename) && strcmp(moviefilename,obj.moviefilename),
          success = true;
          return;
        end
        oldmovieindexfilename = obj.movieindexfilename;
        oldmoviefilename = obj.moviefilename;
        obj.moviefilename = moviefilename;
        obj.movieindexfilename = movieindexfilename;
        %obj.ismovie = ~isempty(moviefilename) && obj.openmovie;
        [success1,msg] = obj.CheckMovies();
        if ~success1,
          obj.moviefilename = oldmoviefilename;
          obj.movieindexfilename = oldmovieindexfilename;
          return;
        end
        [success,msg] = obj.UpdateStatusTable('movie');

      end
      
    end
    

    % ---------------------------------------------------------------------
    function scoreFileName = getScoreFileName(obj)
      scoreFileName = obj.scorefilename;
    end
    
    
    % ---------------------------------------------------------------------
    function [success,msg] = setScoreFileName(obj,sfn)
      assert(iscellstr(sfn) && numel(sfn)==obj.nclassifiers,...
        'Number of score filenames must match number of classifiers.');
      assert(numel(sfn)==numel(unique(sfn)),...
        'Score filenames must be unique.');
      
      tfValid = cellfun(@JLabelData.isValidScoreFileName,sfn);
      if ~all(tfValid)
        error('JLabelData:scoreFileName',...
          'The following are not valid score filenames: %s',...
          civilizedStringFromCellArrayOfStrings(sfn(~tfValid)));
      end
      
      obj.scorefilename = sfn(:)';
      obj.needsave = true;
      [success,msg] = obj.UpdateStatusTable('scores');
    end

    
    % ---------------------------------------------------------------------
    function res = IsRequiredFile(obj,file)
      if obj.isST
        res = any(strcmp(file,{'trx','perframedir','stfeatures'}));
      else
        res = ismember(file,{'trx','perframedir'});
      end
    end
    
   
    % ---------------------------------------------------------------------
    function expdirs_removed = removeExperimentsWithNoLabels(self)
      self.StoreLabelsForCurrentAnimal();  % make sure labels are commited
      nExps=self.nexps;
      markedForRemoval=false(1,nExps);
      for i=1:nExps
        t0sThis=self.labels(i).t0s;
        % t0sThis is a cell array, with one element per labeled target.
        % Each element holds a double array, with the start frames of each
        % labelled bout for that target.
        nBoutsPerLabeledTarget=cellfun(@(a)length(a),t0sThis);
        nBouts=sum(nBoutsPerLabeledTarget);
        markedForRemoval(i)=(nBouts==0);
      end
      expIndicesToRemove=find(markedForRemoval);
      expdirs_removed = self.expdirs(expIndicesToRemove);
      self.RemoveExpDirs(expIndicesToRemove);  
    end  % method

        
    % ---------------------------------------------------------------------
    function [fileExists,fileTimeStamps,missingFileNamesComplete] = ...
        fileOfGivenTypesExistForGivenExps(obj,filetypes,expis)
      % This gets whether the experiment files of filetypes (a cell array
      % of strings) exist, for the experiments indicated by expis.  On
      % return, fileExists is a length(expis) x length(filetypes) logical
      % array, indicating whether the given file type exists for the given
      % experiment.  fileTimeStamps has the same shape, and given the time
      % stamp of each file if it exists, and -inf otherwise.
      % missingFileNamesComplete is a cell array vector of length
      % length(expis).  Each element is a cell array of strings, listing
      % the missing files for that experiment.  These lists are not 
      % truncated to any particular length.  (That's why 'Complete' is in 
      % the name.  Note that this function is a getter---it does not change
      % the object state at all.
      %
      % For filetypes that can represent multiple files:
      %  - fileExists is only true if *all* files are present
      %  - fileTimeStamps is the *oldest* of all timestamps
      %      AL20150109: For fileType='perframedir', fileTimeStamps appears
      %      to be *newest* of all timestamps

      % ST OK
      
      % Process arguments
      if ~exist('filetypes','var') || isempty(filetypes)
        filetypes = obj.filetypes;
      end
      if ~exist('expis','var') || isempty(expis)
        expis = 1:obj.nexps;
      end

      % Check that the file types requested are valid
      [fileTypeIsValid,fileTypeIndices] = ismember(filetypes,obj.filetypes);
      assert(all(fileTypeIsValid),'JLabelData:unknownFileType', ...
        'Internal error: Unknown file type.  Please report to the JAABA developers.');
      
      % Initialize the return vars
      nExpsToCheck = length(expis);
      nFileTypesToCheck = length(fileTypeIndices);
      fileExists = false(nExpsToCheck,nFileTypesToCheck);
      fileTimeStamps = nan(nExpsToCheck,nFileTypesToCheck);
      missingFileNamesComplete = cell(1,nExpsToCheck);
      for i = 1:nExpsToCheck
        missingFileNamesComplete{i} = {};
      end
      
      % Loop through all file types, checking the existence of files of
      % that type.
      for j=1:nFileTypesToCheck
        fileTypeIndex = fileTypeIndices(j);
        fileType = obj.filetypes{fileTypeIndex};
        % loop through experiments
        for i=1:nExpsToCheck
          expi = expis(i);
          if strcmpi(fileType,'perframedir'),
            % fileType of 'perframedir' is a special case, and has to be
            % handled separately
            [perframeFileNames,timestamps] = obj.GetPerframeFiles(expi);
            if isempty(perframeFileNames),
              fileExists(i,j) = false;
              fileTimeStamps(i,j) = -inf;
            else
              perframeFileExists = cellfun(@(s) exist(s,'file'),perframeFileNames);
              fileExists(i,j) = all(perframeFileExists);
              if ~fileExists(i,j) && obj.IsRequiredFile(fileType),
                for indicesOfMissingPerframeFiles = find(~perframeFileExists(:)'),
                  [~,fileNameThis] = myfileparts(perframeFileNames{indicesOfMissingPerframeFiles});
                  missingFileNamesComplete{i}{end+1} = ['perframe_' fileNameThis];
                end
              end
              fileTimeStamps(i,j) = max(timestamps); % AL 20150109: This is newest timestamp, comments above say oldest
            end
          else
            % if fileType is anything besides 'perframedir'
            % check for existence of current file(s)
            [fn,fTS,tffound] = obj.GetFile(fileType,expi);
            assert(iscell(fn)||ischar(fn),'Legacy check.');
            fileExists(i,j) = all(tffound);
            fileTimeStamps(i,j) = min(fTS);
            if ~fileExists(i,j) && obj.IsRequiredFile(fileType)
              missingFileNamesComplete{i}{end+1} = fileType;
            end
          end
        end
      end
    end  % method


    % ---------------------------------------------------------------------
    function [allRequiredFilesExist, ...
              missingFilesCanBeGenerated, ...
              oneMissingFileTypeThatCantBeGenerated] = allRequiredFilesExist(obj,fileExists)
      % Check whether all required files exist for all experiments.
      % fileExists is optional, and if provided should be a complete matrix
      % of file existance information, e.g. as returned by
      % fileExists=obj.fileOfGivenTypesExistForGivenExps().  If this
      % argument is not provided, this method calls
      % fileOfGivenTypesExistForGivenExps() to get current file existance
      % information.  (The idea is that if you just called
      % fileOfGivenTypesExistForGivenExps(), you can give its output to
      % this method, and then this method won't have to figure that stuff
      % out again.)
      % On return, allRequiredFilesExist is a logical scalar.  If
      % allRequiredFilesExist is false, missingFilesCanBeGenerated is a
      % logical scalar indicating whether the missing files can be
      % generated.  If both allRequiredFilesExist and
      % missingFilesCanBeGenerated are false,
      % oneMissingFileTypeThatCantBeGenerated gives the name of one file
      % type that cannot be generated.  Note that this function is a
      % getter---it does not change the object state at all.

      % ST OK
      
      % Get detailed information about file existance, if not given as
      % argument
      if ~exist('fileExists','var')
        fileExists = obj.fileOfGivenTypesExistForGivenExps();
          % with no args, computed for all file types, all exps
      end
      
      nExp = obj.nexps;
      nType = numel(obj.filetypes);
      assert(isequal(size(fileExists),[nExp nType]));
            
      % initialize outputs
      allRequiredFilesExist = true;
      missingFilesCanBeGenerated = [];
      oneMissingFileTypeThatCantBeGenerated = '';
      
      for iType = 1:nType
        type = obj.filetypes{iType};
        tfReqd = obj.IsRequiredFile(type);
        if tfReqd
          for iExp = 1:nExp
            if ~fileExists(iExp,iType)
              allRequiredFilesExist = false;
              % if furthermore file can't be generated, then not fixable
              missingFilesCanBeGenerated = JLabelData.CanGenerateFile(type);
              if ~missingFilesCanBeGenerated
                oneMissingFileTypeThatCantBeGenerated = type;
                return;
              end
            end
          end
        end
      end
    end 

    % ---------------------------------------------------------------------
    function [success,msg] = SetPerFrameDir(obj,perframedir)
      % [success,msg] = SetPerFrameDir(obj,perframedir)
      % Sets the per-frame directory name within the experiment directory.
      % Currently, this does not change the cached per-frame data or check
      % that all the per-frame files necessary are within the directory
      % (TODO).
      
      success = false; msg = '';
      
      if ischar(perframedir),
        if ischar(obj.perframedir) && strcmp(perframedir,obj.perframedir),
          success = true;
          return;
        end
        
        obj.perframedir = perframedir;
        
        % TODO: check per-frame directories are okay, remove bad
        % experiments
        
        [success,msg] = obj.UpdateStatusTable('perframedir');
      end
    end
    
    % ---------------------------------------------------------------------
    function [success,msg] = SetClipsDir(obj,clipsdir)
    % [success,msg] = SetClipsDir(obj,clipsdir)
    % Sets the clips directory name within the experiment directory.
      
      success = true;
      msg = '';

      if ischar(clipsdir),
        for i = 1:numel(obj.expdirs),
          %clipsdircurr = fullfile(obj.expdirs{i},clipsdir);
%           if exist(obj.expdirs{i},'dir') && ~exist(clipsdircurr,'dir'),
%             mkdir(clipsdircurr);
%           end
        end
        if ischar(obj.clipsdir) && strcmp(clipsdir,obj.clipsdir),
          success = true;
          return;
        end

        obj.clipsdir = clipsdir;        
        [success,msg] = obj.UpdateStatusTable('clipsdir');
      end
      
    end

    function [success,msg] = SetTrkFileName(obj,trkfilename)
    % [success,msg] = SetClipsDir(obj,clipsdir)
    % Sets the trk filename within the experiment directory.
      
      success = true;
      msg = '';

      if ischar(trkfilename),
        for i = 1:numel(obj.expdirs),
          trkcurr = fullfile(obj.expdirs{i},trkfilename);
          if ~exist(trkcurr,'file')
            success = false;
            msg = sprintf('Trk file %s does not exist',trkcurr);
          end
        end
        if ischar(obj.trkfilename) && strcmp(trkfilename,obj.trkfilename),
          success = true;
          return;
        end

        obj.trkfilename = trkfilename;        
        [success,msg] = obj.UpdateStatusTable('trk');
      elseif iscell(trkfilename),
        if iscell(obj.trkfilename) && all(strcmp(trkfilename,obj.trkfilename)),
          success = true;
          return;
        end
        obj.trkfilename = trkfilename;        
        [success,msg] = obj.UpdateStatusTable('trk');
      end
      
    end

    
    % ---------------------------------------------------------------------
    function [success,msg] = loadPerframeData(obj,expi,indicesOfTargets)
      % Loads the per-frame features for experiment expi and target
      % iTarget into obj.perframedata and obj.perframeunits      
      
      % Test for various degenerate cases 
      if isempty(indicesOfTargets)
        success=true;  msg='';
        return
      end
      iTarget=indicesOfTargets(1);
      if isempty(expi) || ~((1<=expi)&&(expi<=obj.nexps))
        success=true;  msg='';
        return
      end        
      % If we got here, expi and iTarget are non-degenerate
      perframeFileNameList = obj.GetPerframeFiles(expi);
      nPerframeFeatures=numel(obj.allperframefns);
      perframedata=cell(1,nPerframeFeatures);
      perframeunits=cell(1,nPerframeFeatures);      
      
      for j = 1:nPerframeFeatures,
        perframeFileName=perframeFileNameList{j};
        if ~exist(perframeFileName,'file'),
          success = false;
          msg = sprintf('Per-frame data file %s does not exist',perframeFileNameList{j});
          return;
        end
      end
      
      stInfo = obj.stInfo;
      all_update = false(1,nPerframeFeatures);
      for j = 1:nPerframeFeatures % TODO: change to parfor if it is for
        perframeFileName=perframeFileNameList{j};
        [ftmp,funits,update] = readPFData(perframeFileName,iTarget,stInfo);
%         tmp = load(perframeFileName);
%         assert(isequal(ftmp,tmp.data(iTarget)));
%         assert(isequal(funits,tmp.units));
        
        perframedata{j} = ftmp{1};
        perframeunits{j} = funits;
        all_update(j) = update;        
      end
      if any(all_update)
          warndlg('The perframe features were generated with different parameters. To regenerate them with new parameters, delete them. JAABA, for now, will continue to use existing ones.'); 
      end
      obj.perframedata = perframedata;
      obj.perframeunits = perframeunits;
      success=true;
      msg='';
    end  % method
    
    
    % ---------------------------------------------------------------------
    function [perframedata,T0,T1] = GetPerFrameData(obj,expi,flies,prop,T0,T1)
    % [perframedata,T0,T1] = GetPerFrameData(obj,expi,flies,prop,T0,T1)
    % Returns the per-frame data for the input experiment, flies, and
    % property. 

      if ischar(prop),
        prop = find(strcmp(prop,obj.allperframefns),1);
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
      pffile = fullfile(perframedir,[obj.allperframefns{prop},'.mat']);
      perframedata = readPFData(pffile,flies(1),obj.stInfo);
      perframedata = perframedata{1};
      if nargin < 5,
        T0 = max(obj.GetTrxFirstFrame(expi,flies));
        % TODO: generalize to multi-fly
        T1 = T0 + numel(perframedata) - 1;
        return;
      end
      off = 1 - obj.GetTrxFirstFrame(expi,flies);
      i0 = T0 + off;
      i1 = T1 + off;
      perframedata = perframedata(i0:i1);
    end

    
    % ---------------------------------------------------------------------
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
    
  end
    
  methods (Static)
    
    % ---------------------------------------------------------------------
    function result = isValidScoreFileName(scoreFileName)
      [path,baseName,ext]=fileparts_platind(scoreFileName);
      if ~isempty(path) ,
        result=false;
      elseif ~isequal(ext,'.mat')
        result=false;
      else
        result=~isempty(regexp(baseName,'^[a-zA-Z_0-9]+$','once'));
      end
    end
    
    % movie, trx, and perframedir are required for each experiment
    % perframedir can be generated
    function res = CanGenerateFile(file)
      res = ismember(file,{'perframedir'});
    end
    
    % which files should go in the output directory
    function res = IsOutputFile(file)
      res = ismember(file,{'label','clipsdir','scores','gt_label'});
    end
   
  end
  
  methods (Static,Access=private)
    
    function [tffound,filename,timestamp] = GetFileRaw(fileType,parentDir,fileNameLocal)
      % Look for a filename in a parent dir.
      % fileType: The fileType enum used here in JLabelData. This is used
      % only for a quirk regarding .lnk/.seq files
      % fileNameLocal: char, single short filename.
      %
      % tffound: scalar logical.
      %   - If true, filename is the full filename and timestamp is the filesystem timestamp.
      %   - If false, filename is the file-that-was-tried; and timestamp=-inf.
      
      assert(isscalar(parentDir));
      
      for j = 1:numel(parentDir),
        ddir = parentDir{j};
        
        filename = fullfile(ddir,fileNameLocal);
        if exist(filename,'file'),
          tmp = dir(filename);
          timestamp = tmp.datenum;
          tffound = true;
          return;
        end
        % check for lnk files
        if ispc && exist([filename,'.lnk'],'file'),
          isseq = ~isempty(regexp(fileNameLocal,'\.seq$','once'));
          % for seq file, just keep the soft link, get_readframe_fcn will
          % deal with it
          [actualfilename,didfind] = GetPCShortcutFileActualPath(filename);
          if didfind,
            tmp = dir(actualfilename);
            timestamp = tmp.datenum;
            if ~isseq || ~strcmpi(fileType,'movie'),
              filename = actualfilename;
            end
            tffound = true;
            return;
          end
        end
      end
      
      tffound = false;
      % filename is the file-that-was-tried
      timestamp = -inf;
    end
    
  end 
  
  %% Tracking info
  
  methods

    % ---------------------------------------------------------------------
    function nflies = GetNumFlies(obj,expi)
      nflies = obj.nflies_per_exp(expi);
    end
    
    
    % --------------------------------------------------------------------------
    function [success,msg] = GetTrxInfo(obj,expi,canusecache,trx)
    % [success,msg] = GetTrxInfo(obj,expi)
    % Fills in nflies_per_exp, firstframes_per_exp, and endframes_per_exp
    % for experiment expi. This may require loading in trajectories. 
      success = true;
      msg = '';
      if nargin < 3,
        canusecache = true;
      end
%       canusecache = false;
      istrxinput = nargin >= 4;
      
      obj.SetStatus('Reading trx info for experiment %s',obj.expdirs{expi});
      if numel(obj.nflies_per_exp) < expi || ...
          numel(obj.sex_per_exp) < expi || ...
          numel(obj.frac_sex_per_exp) < expi || ...
          numel(obj.firstframes_per_exp) < expi || ...
          numel(obj.endframes_per_exp) < expi || ...
          isnan(obj.nflies_per_exp(expi)),
        if ~istrxinput,

          trxfile = fullfile(obj.expdirs{expi},obj.GetFileName('trx'));
          if ~exist(trxfile,'file'),
            msg = sprintf('Trx file %s does not exist, cannot count flies',trxfile);
            success = false;
            return;
          else
          
            if isempty(obj.expi) || obj.expi == 0,
              % TODO: make this work for multiple flies
              obj.setCurrentTarget(expi,1);
              trx = obj.trx;
            elseif canusecache && expi == obj.expi,
              trx = obj.trx;
            else
              trx = load_tracks(trxfile);
%               try
                % REMOVE THIS
%                 global CACHED_TRX; %#ok<TLEV>
%                 global CACHED_TRX_EXPNAME; %#ok<TLEV>
%                 if isempty(CACHED_TRX) || isempty(CACHED_TRX_EXPNAME) || ...
%                     ~strcmp(obj.expnames{expi},CACHED_TRX_EXPNAME),
%                   hwait = mywaitbar(0,sprintf('Loading trx to determine number of flies for %s',...
%                     obj.expnames{expi}),'interpreter','none');
%                   trx = load_tracks(trxfile);
%                   if ishandle(hwait), delete(hwait); end
%                   CACHED_TRX = trx;
%                   CACHED_TRX_EXPNAME = obj.expnames{expi};
%                 else
% %                   fprintf('DEBUG: Using CACHED_TRX. REMOVE THIS\n');
%                   trx = CACHED_TRX;
%                 end
%               catch ME,
%                 msg = sprintf(['Could not load trx file for experiment %s '...
%                     'to count flies: %s'],obj.expdirs{expi},getReport(ME));
%               end
            end
          end
        end
        
        obj.nflies_per_exp(expi) = numel(trx);
        obj.firstframes_per_exp{expi} = [trx.firstframe];
        obj.endframes_per_exp{expi} = [trx.endframe];

        obj.hassex = obj.hassex || isfield(trx,'sex');
        
        % store sex info
        tmp = repmat({nan},[1,numel(trx)]);
        obj.frac_sex_per_exp{expi} = struct('M',tmp,'F',tmp);
        obj.sex_per_exp{expi} = repmat({'?'},[1,numel(trx)]);
        if isfield(trx,'sex'),
          obj.hasperframesex = iscell(trx(1).sex);
          if obj.hasperframesex,
            for fly = 1:numel(trx),
              n = numel(trx(fly).sex);
              nmale = nnz(strcmpi(trx(fly).sex,'M'));
              nfemale = nnz(strcmpi(trx(fly).sex,'F'));
              obj.frac_sex_per_exp{expi}(fly).M = nmale/n;
              obj.frac_sex_per_exp{expi}(fly).F = nfemale/n;
              if nmale > nfemale,
                obj.sex_per_exp{expi}{fly} = 'M';
              elseif nfemale > nmale,
                obj.sex_per_exp{expi}{fly} = 'F';
              else
                obj.sex_per_exp{expi}{fly} = '?';
              end
            end
          else
            for fly = 1:numel(trx),
              obj.sex_per_exp{expi}{fly} = trx(fly).sex;
              if strcmpi(trx(fly).sex,'M'),
                obj.frac_sex_per_exp{expi}(fly).M = 1;
                obj.frac_sex_per_exp{expi}(fly).F = 0;
              elseif strcmpi(trx(fly).sex,'F'),
                obj.frac_sex_per_exp{expi}(fly).M = 0;
                obj.frac_sex_per_exp{expi}(fly).F = 1;
              end
            end
          end
        end
      end
      if isfield(trx,'arena')
        obj.hasarenaparams(expi) = true;
      else
        obj.hasarenaparams(expi) = false;
      end
      
      obj.ClearStatus();
      
    end

    
    % --------------------------------------------------------------------------
    function out = GetTrxValues(obj,infoType,expi,flies,ts)
    % A generic function that return track info.

      if numel(expi) ~= 1,
        error('expi must be a scalar');
      end
      
      if expi ~= obj.expi,
        % TODO: generalize to multiple flies
        [success,msg] = obj.setCurrentTarget(expi,1);
        if ~success,
          error('Error loading trx for experiment %d: %s',expi,msg);
        end
      end

      if nargin < 4,     % No flies given
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

      
      elseif nargin < 5, % No ts given
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
    
    
    % ---------------------------------------------------------------------
    function pos = GetTrxPos1(varargin)
    % [x,y,theta,a,b] = GetTrxPos1(obj,expi,fly,ts)
    % Returns the position for the input experiment, SINGLE fly, and
    % frames. If ts is not input, then all frames are returned. 

      % moved to separate file so that function could be easily modified
      pos = JLabelData_GetTrxPos(varargin{:});

    end

    
    % ---------------------------------------------------------------------
    function sex = GetSex(obj,expi,fly,ts,fast)
    % x = GetSex(obj,expi,fly,ts)
    % Returns the sex for the input experiment, SINGLE fly, and
    % frames. If ts is not input, then all frames are returned. 

      if ~obj.hassex,
        sex = '?';
        return;
      end
      
      if nargin < 5,
        fast = false;
      end
      
      if ~obj.hasperframesex || fast,
        sex = obj.sex_per_exp{expi}(fly);
        return;
      end
      
      if expi ~= obj.expi,
        % TODO: generalize to multiple flies
        [success,msg] = obj.setCurrentTarget(expi,fly);
        if ~success,
          error('Error loading trx for experiment %d: %s',expi,msg);
        end
      end
      
      if nargin < 4,
        sex = obj.trx(fly).sex;
        return;
      end
      
      sex = obj.trx(fly).sex(ts + obj.trx(fly).off);

    end

    
    % ---------------------------------------------------------------------
    function sex = GetSex1(obj,expi,fly,t)
    % x = GetSex1(obj,expi,fly,t)
    % Returns the sex for the input experiment, SINGLE fly, and
    % SINGLE frame. 

      if ~obj.hassex,
        sex = '?';
        return;
      end
            
      if ~obj.hasperframesex,
        sex = obj.sex_per_exp{expi}(fly);
        if iscell(sex),
          sex = sex{1};
        end
        return;
      end
      
      if expi ~= obj.expi,
        % TODO: generalize to multiple flies
        [success,msg] = obj.setCurrentTarget(expi,fly);
        if ~success,
          error('Error loading trx for experiment %d: %s',expi,msg);
        end
      end
            
      sex = obj.trx(fly).sex{t + obj.trx(fly).off};

    end
    
    
    % ---------------------------------------------------------------------
    function sexfrac = GetSexFrac(obj,expi,fly)
    % x = GetSexFrac(obj,expi,fly)
    % Returns a struct indicating the fraction of frames for which the sex
    % of the fly is M, F

      sexfrac = obj.frac_sex_per_exp{expi}(fly);

    end
    
    
    % ---------------------------------------------------------------------
    function t0 = GetTrxFirstFrame(obj,expi,flies)
    % t0 = GetTrxFirstFrame(obj,expi,flies)
    % Returns the firstframes for the input experiment and flies. If flies
    % is not input, then all flies are returned. 
      if numel(expi) ~= 1,
        error('expi must be a scalar');
      end
      
      if nargin < 3,
        t0 = double(obj.firstframes_per_exp{expi});
        return;
      end

      t0 = double(obj.firstframes_per_exp{expi}(flies));
      
    end

    
    % ---------------------------------------------------------------------
    function t1 = GetTrxEndFrame(obj,expi,flies)
    % t1 = GetTrxEndFrame(obj,expi,flies)
    % Returns the endframes for the input experiment and flies. If flies
    % is not input, then all flies are returned. 

      if numel(expi) ~= 1,
        error('expi must be a scalar');
      end
      
      if nargin < 3,
        t1 = double(obj.endframes_per_exp{expi});
        return;
      end

      t1 = double(obj.endframes_per_exp{expi}(flies));
      
    end
    
    
    % ---------------------------------------------------------------------
    function firstframes = GetFirstFrames(obj,expi,flies)
      
      if nargin < 2,
        firstframes = obj.firstframes_per_exp;
      elseif nargin < 3,
        firstframes = obj.firstframes_per_exp(expi);
      else
        firstframes = obj.firstframes_per_exp{expi}(flies);
      end
      
    end

    
    % ---------------------------------------------------------------------
    function endframes = GetEndFrames(obj,expi,flies)
      
      if nargin < 2,
        endframes = obj.endframes_per_exp;
      elseif nargin < 3,
        endframes = obj.endframes_per_exp(expi);
      else
        endframes = obj.endframes_per_exp{expi}(flies);
      end
      
    end

    
    % ---------------------------------------------------------------------
    function minFirstframes = GetMinFirstFrame(obj)
      minFirstframes = min(obj.firstframes_per_exp{obj.expi});
    end

    
    % ---------------------------------------------------------------------
    function maxEndframe = GetMaxEndFrame(obj)
      maxEndframe = max(obj.endframes_per_exp{obj.expi});
    end
    
  end
  
  methods (Static)
    
     % --------------------------------------------------------------------------
    function [nFlies,firstFrames,endFrames,hasArenaParams,hasSex,fracSex,sex,hasPerFrameSex] = ...
        readTrxInfoFromFile(trxFileName,varargin)
      % Read the trx file
      [trx,~,success] = load_tracks(trxFileName);
      if ~success,
        error('JAABA:JLabelData:readTrxInfoFromFile:errorReadingTrxFile', ...
              'Unable to read .trx file');
      end
      
      % Set the returned things which can be read directly from the .trx
      % file
      nFlies = numel(trx);
      firstFrames = [trx.firstframe];
      endFrames = [trx.endframe];
      hasArenaParams=isfield(trx,'arena');
      hasSex = isfield(trx,'sex');
      hasPerFrameSex = false;
        
      % Compute fracSex and sex from the trx file contents
      fracSex = struct('M',repmat({nan},[1 nFlies]), ...
                       'F',repmat({nan},[1 nFlies]));
      sex = repmat({'?'},[1 nFlies]);
      if hasSex,
        % this track file has sex
        if nFlies==0,
          hasPerFrameSex=true;  % vacuous truth
        else
          hasPerFrameSex = iscell(trx(1).sex);
        end
        if hasPerFrameSex,
          % this trx file has per-frame sex
          for iFly = 1:nFlies,
            n = numel(trx(iFly).sex);
            nMale = nnz(max(strcmpi(trx(iFly).sex,'M'),strcmpi(trx(iFly).sex,'male')));
            nFemale = nnz(max(strcmpi(trx(iFly).sex,'F'),strcmpi(trx(iFly).sex,'female')));
            fracSex(iFly).M = nMale/n;
            fracSex(iFly).F = nFemale/n;
            if nMale > nFemale,
              sex{iFly} = 'M';
            elseif nFemale > nMale,
              sex{iFly} = 'F';
            else
              sex{iFly} = '?';
            end
          end
        else
          % this trx files does not have per-frame sex
          for iFly = 1:nFlies,
            sex{iFly} = trx(iFly).sex;
            if strcmpi(trx(iFly).sex,'M'),
              fracSex(iFly).M = 1;
              fracSex(iFly).F = 0;
            elseif strcmpi(trx(iFly).sex,'F'),
              fracSex(iFly).M = 0;
              fracSex(iFly).F = 1;
            end
          end
        end
      end
    end  % method   
    
  end
  
  %% WindowData
  
  methods
    
    
    function res = HasWindowdata(self,iCls)
      %MERGESTUPDATED
      names = self.iCls2LblNames{iCls};
      res = self.getAtLeastOneNormalLabelOfEachClassExists(names) && ...
            ~isempty(self.windowdata(iCls).X);
    end

    
    % ---------------------------------------------------------------------
    function idx = FlyNdx(obj,expi,flies,iCls)
      %MERGESTUPDATED      
      
      if isempty(obj.windowdata(iCls).exp)
        idx = [];
        return;
      end
      idx = obj.windowdata(iCls).exp==expi & all(bsxfun(@eq,obj.windowdata(iCls).flies,flies),2);
    end
    
    
    % ---------------------------------------------------------------------
    function object = createPreLoadWindowDataObj(obj,expi,flies,iCls)
      % create a simple struct for PreLoadWindowData, which is parfor-ed
      
      obj.CheckExp(expi); 
      obj.CheckFlies(flies);

      object = struct();
      object.isempty_fieldnames_windowfeaturesparams = isempty(fieldnames(obj.windowfeaturesparams{iCls})); % ALTODO can almost definitely be removed
      object.GetTrxFirstFrame = obj.GetTrxFirstFrame(expi,flies);
      object.GetTrxEndFrame = obj.GetTrxEndFrame(expi,flies);
      object.windowdatachunk_radius = obj.windowdatachunk_radius;
      object.not_isempty_windowdata_exp = ~isempty(obj.windowdata(iCls).exp);
      object.windowdata_t_flyndx = obj.windowdata(iCls).t(obj.FlyNdx(expi,flies,iCls));
      %object.gettrxfirstframe = obj.GetTrxFirstFrame(expi,flies);
      object.curperframefns = obj.curperframefns{iCls};
      object.allperframefns = obj.allperframefns;
      object.perframefile = obj.GetPerframeFiles(expi);
      object.windowfeaturescellparams = obj.windowfeaturescellparams{iCls};
      
      % ALTODO factor this out (see also ComputeWindowDataChunk) into eg ensurePerFrameFiles.
      % Call the result before calling createPreLoadWindowDataObj. This
      % should also eliminate perframefn consistency issues in caller
      % (Preloadperi)
      for j = 1:numel(object.curperframefns)
        fn = object.curperframefns{j};

        % get per-frame data
        ndx = find(strcmp(fn,object.allperframefns));
        assert(~isempty(ndx),...
          'Internal error: There is at least one per-frame feature in the vocabulary (%s) that is not in the subdialect.',fn);

        if ~exist(object.perframefile{ndx},'file'),
          if ~isdeployed 
            if isempty(obj.GetGenerateMissingFiles())
              res = questdlg(sprintf(['Experiment %s is missing some perframe files '...
                '(%s, possibly more). Generate now?'],obj.expnames{expi},object.perframefile{ndx}),...
                'Generate missing files?','Yes','Cancel','Yes');
              if strcmpi(res,'Yes');
                obj.SetGenerateMissingFiles();
              end
            else 
              res = fif(obj.GetGenerateMissingFiles(),'Yes','No');
            end
          else
            res = 'Yes';
          end

          if strcmpi(res,'Yes'),
            for ndx = 1:obj.nexps  
              [success1,msg1] = obj.GenerateMissingFiles(ndx);
              if ~success1
                error(msg1);
              end
            end
          else
            error('Cannot compute window data for %s.',obj.expnames{expi});
          end
        end
      end
    end
      
      
    function [success,msg,t0,t1,X,feature_names] = ComputeWindowDataChunk(obj,expi,flies,t,mode,forceCalc,iCls)
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

      % MERGESTUPDATED
    
      success = false; msg = '';  %#ok
      
      if ~exist('mode','var') || isempty(mode), mode = 'center'; end
      if ~exist('forceCalc','var') || isempty(forceCalc), forceCalc = false; end
      
      % Check if the features have been configured.
      % I really don't like this.  The JLabelData is pretty close to being a
      % model in the MVC sense.  As such, it shouldn't be creating a view.
      % Not clear to me what the best way to fix this is, though.
      % --ALT, Feb 4, 2013
      if isempty(fieldnames(obj.windowfeaturesparams{iCls}))
        figureJLabel = findall(0,'tag','figure_JLabel');
        figureSelectFeatures = SelectFeatures(figureJLabel);
        uiwait(figureSelectFeatures);
        if isempty(fieldnames(obj.windowfeaturesparams{iCls}))
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
      if ~isempty(obj.windowdata(iCls).exp) && ~forceCalc,
        tscomputed = obj.windowdata(iCls).t(obj.FlyNdx(expi,flies,iCls));
        tscomputed = tscomputed(tscomputed >= t0 & tscomputed <= t1);
        docompute(tscomputed+off) = false;
      end
      
      X = single([]);
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
      
      curperframefns = obj.curperframefns{iCls};
      allperframefns = obj.allperframefns;
      perframeInMemory = ~isempty(obj.flies) && obj.IsCurFly(expi,flies);
      perframedata_all = obj.perframedata;
      perframefile = obj.GetPerframeFiles(expi);
      x_curr_all = cell(1,numel(curperframefns));
      feature_names_all = cell(1,numel(curperframefns));
      windowfeaturescellparams = obj.windowfeaturescellparams{iCls};
      
      % loop through per-frame fields to check that they exist.
      % ALTODO: factor this out, see createPreLoadWindowDataObj
      for j = 1:numel(curperframefns),
        fn = curperframefns{j};
        
        % get per-frame data
        ndx = find(strcmp(fn,allperframefns));
        if isempty(ndx),
          success = false;
          msg = sprintf('Internal error: There is at least one per-frame feature in the vocabulary (%s) that is not in the subdialect.',fn);
          return;
        end
        
        if ~exist(perframefile{ndx},'file'),
          if ~isdeployed 
            if isempty(obj.GetGenerateMissingFiles())
              res = questdlg(sprintf(['Experiment %s is missing some perframe files '...
                '(%s, possibly more). Generate now?'],obj.expnames{expi},perframefile{ndx}),...
                'Generate missing files?','Yes','Cancel','Yes');
              if strcmpi(res,'Yes');
                obj.SetGenerateMissingFiles();
              end
            else 
              res = fif(obj.GetGenerateMissingFiles(),'Yes','No');
            end
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
            msg = sprintf('Cannot compute window data for %s ',obj.expnames{expi});
            return;
          end
        end
        
      end
      
      parfor j = 1:numel(curperframefns),
      %for j = 1:numel(curperframefns),
        fn = curperframefns{j};
%        fprintf('Computing window data for per-frame feature %d: %s\n',j,fn);
        
        % get per-frame data
        ndx = find(strcmp(fn,allperframefns));
        if perframeInMemory,
          perframedata = perframedata_all{ndx};  %#ok
        else
          perframedata = load(perframefile{ndx});  %#ok
          perframedata = perframedata.data{flies(1)};  %#ok
        end
        
        i11 = min(i1,numel(perframedata));
        [x_curr,feature_names_curr] = ...
          ComputeWindowFeatures(perframedata,windowfeaturescellparams.(fn){:},'t0',i0,'t1',i11);  %#ok
        if any(imag(x_curr(:)))
          fprintf('Feature values are complex, check input\n');
        end
        
        if i11 < i1,
          x_curr(:,end+1:end+i1-i11) = nan;
        end
        
        x_curr_all{j} = single(x_curr);
        feature_names_all{j} = feature_names_curr;        
      end
      
      feature_names = cell(1,numel(curperframefns));
      for j = 1:numel(curperframefns),
        fn = curperframefns{j};
        x_curr = x_curr_all{j};
        feature_names_curr = feature_names_all{j};
        % add the window data for this per-frame feature to X
        nold = size(X,1);
        nnew = size(x_curr,2);
        if nold > nnew,
          warning(['Number of examples for per-frame feature %s does not match number '...
            'of examples for previous features'],fn);
          x_curr(:,end+1:end+nold-nnew) = nan;
        elseif nnew > nold && ~isempty(X),
          warning(['Number of examples for per-frame feature %s does not match number '...
            'of examples for previous features'],fn);
          X(end+1:end+nnew-nold,:) = nan;
        end
        X(:,end+1:end+size(x_curr,1)) = x_curr';
        % add the feature names
        if nargout>5
          feature_names{j} = cellfun(@(s) [{fn},s],feature_names_curr,'UniformOutput',false); 
        end
      end
      feature_names = [feature_names{:}];
      %       catch ME,
      %         msg = getReport(ME);
      %         return;
      %       end
      %X = single(X);
      success = true;
     
    end
   
    
    % ---------------------------------------------------------------------
    function initWindowData(obj)
      % Inits windowdata and predictblocks. Clears predictions.
      
      %MERGESTUPDATED

      obj.windowdata = WindowData.windowdata(obj.nclassifiers);
      useselfeatures = obj.selFeatures(1).use;
      for ndx = 1:obj.nclassifiers
        obj.selFeatures(ndx) = SelFeatures.createEmpty();
      end
      [obj.selFeatures(:).use] = deal(useselfeatures);
      if useselfeatures,
        [obj.selFeatures(:).do] = deal(true);
      end
      obj.predictblocks = Predict.predictblocks(obj.nclassifiers);      
      obj.UpdatePredictedIdx();
    end
    
    
    % ---------------------------------------------------------------------
    function [success,msg] = PreLoadPeriLabelWindowData(obj)
      % [success,msg] = PreLoadPeriLabelWindowData(obj)
      % This function precomputes any missing window data for all labeled
      % training examples by calling PreLoadWindowData on all labeled frames.
      %
      % Effect: Update obj.predictblocks, obj.windowdata
      % SideEffect: In obscure situations, obj.curperframefns can get
      % updated
      
      %MERGESTUPDATED
      
      for expi = 1:obj.nexps
        for iCls = 1:obj.nclassifiers
          
          obj.SetStatus('Computing windowdata for %s: %s',...
            obj.expnames{expi},obj.labelnames{iCls});
          obj.CheckExp(expi);
          
          flies_curr = obj.labels(expi).flies; % labeled flies for this exp (ANY behavior)
          Nfliescurr = size(flies_curr,1);
          
          % load up state for parfor-windowdata computation
          labelIdxVals = cell(1,Nfliescurr); % vectors, labelidx.vals for this exp/fly/classifier
          labelIdxImp = cell(1,Nfliescurr); % vectors, labelidx.imp for this exp/fly/classifier
          labelIdxT0 = cell(1,Nfliescurr); % scalars, labelidx T0 offsets for this exp/fly
          missingts = cell(1,Nfliescurr); % vectors, labeled frames for this classiifer which are not in this classifier's windowdata
          object = cell(1,Nfliescurr); % structs, dummy object for this exp/fly/classifier
          for flyi = 1:Nfliescurr
            flies = flies_curr(flyi,:);  % BJA: is this ever 2-D ?
            obj.CheckFlies(flies);
            
            [labelIdx,labelIdxT0{flyi}] = obj.GetLabelIdx(expi,flies);
            assert(isequal(size(labelIdx.vals,1),size(labelIdx.imp,1),obj.nclassifiers));
            labelIdxVals{flyi} = labelIdx.vals(iCls,:);
            labelIdxImp{flyi} = labelIdx.imp(iCls,:);
%             object{flyi} = obj.createPreLoadWindowDataObj(expi,flies,iCls);
%           MK: moved this below so that we compute the object only if
%           there are some ts that are missing. The object creating tends
%           to be slow.

            % Find all labeled frames for this exp/classifier/fly
            if obj.nclassifiers==1
              % Legacy codepath: build up ts from labels_curr. This should
              % end up the same as labelidxVals{flyi}; assert this to
              % verify.
              ts = zeros(1,0);
              labels_curr = obj.GetLabels(expi,flies);
              for j = 1:numel(labels_curr.t0s),
                ts = [ts,labels_curr.t0s(j):labels_curr.t1s(j)-1]; %#ok<AGROW>
              end
              assert(isequal(sort(ts+double(labelIdx.off)),sort(find(labelIdxVals{flyi}))));
            else
              % Just use labelIdxVals{flyi} and trust that it is the right
              % thing.
              ts = find(labelIdxVals{flyi})-labelIdx.off;
            end
              
            % Determine which frames are missing windowdata
            if isempty(obj.windowdata(iCls).exp)
              missingts{flyi} = ts;
            else
              idxcurr = obj.FlyNdx(expi,flies,iCls);
              tscurr = obj.windowdata(iCls).t(idxcurr); % object{flyi}.windowdata_t_flyndx
              t0_labelidx = labelIdxT0{flyi};
              obj.windowdata(iCls).labelidx_new(idxcurr) = labelIdxVals{flyi}(tscurr-double(t0_labelidx)+1);
              obj.windowdata(iCls).labelidx_imp(idxcurr) = labelIdxImp{flyi}(tscurr-double(t0_labelidx)+1);
              missingts{flyi} = setdiff(ts,tscurr);
            end

          end
          
          % Trim unnecessary flies.
          toremove = cellfun(@(x) isempty(x),missingts);
          flies_curr = flies_curr(~toremove);
          Nfliescurr = size(flies_curr,1);
          labelIdxVals = labelIdxVals(~toremove) ; % vectors, labelidx.vals for this exp/fly/classifier
          labelIdxImp = labelIdxImp(~toremove); % vectors, labelidx.imp for this exp/fly/classifier
          labelIdxT0 = labelIdxT0(~toremove); % scalars, labelidx T0 offsets for this exp/fly
          missingts = missingts(~toremove); % vectors, labeled frames for this classiifer which are not in this classifier's windowdata
          object = object(~toremove); % structs, dummy object for this exp/fly/classifier
         
          if Nfliescurr==0, continue; end
          
          for flyi = 1:Nfliescurr
            flies = flies_curr(flyi,:); 
            object{flyi} = obj.createPreLoadWindowDataObj(expi,flies,iCls);
          end
          
          curperframefns = obj.curperframefns{iCls};
          allperframefns = obj.allperframefns;
          Ncpff = numel(curperframefns);
          parfor_predictblocks = cell(1,Ncpff); % parfor_predictblocks{perframei}{flyi}
          parfor_windowdata = cell(1,Ncpff); % parfor_windowdata{perframei}{flyi}
          obj_getperframefiles = obj.GetPerframeFiles(expi);
          assert(numel(allperframefns)==numel(obj_getperframefiles));
          
          % AL 20140313 Ensure consistency of curperframefns/allperframefns
          % lists. At the moment there appears to be an obscure potential
          % codepath from CreateObjectForComputeWindowDataChunk to
          % removeArenaPFs which could lead to mutation of the perframefn
          % lists across object{1}, object{2}, etc.
          % ALTODO this will be cleaned up createPreLoadWindowDataObj is
          % factored.
          for flyi = 1:Nfliescurr
            object{flyi}.curperframefns = curperframefns;
            object{flyi}.allperframefns = allperframefns;
          end
          
          usecacheperframe = false;
          if Nfliescurr ==1
            usecacheperframe = obj.IsCurFly(expi,flies_curr);
          end
          cacheperframedata = obj.perframedata;
          
          [~,pfidx] = ismember(curperframefns,allperframefns);
          tmp_cacheperframedata = cacheperframedata(pfidx);
          
          stInfo = obj.stInfo;
          parfor perframei = 1:Ncpff % TODO: Switch back to parfor
            perframefn = curperframefns{perframei};
            ndx = find(strcmp(perframefn,allperframefns));
            
            % KB: usecacheperframe was false when Nflies_curr was 0, so was
            % loading in per-frame data            
            if Nfliescurr > 0,
              if usecacheperframe
                parperframedata = tmp_cacheperframedata(perframei);
              else
                parperframedata = readPFData(obj_getperframefiles{ndx},flies_curr,stInfo);
              end
            else
              parperframedata = [];
            end

            parfor_predictblocks{perframei} = cell(1,Nfliescurr);
            parfor_windowdata{perframei} = cell(1,Nfliescurr);
            
            for flyi = 1:Nfliescurr
%               assert(isequal(perframedata{flies},parperframedata{flyi}));
              [tmpsuccess,~,predictblocks,windowdata] = ...
                PreLoadWindowData(object{flyi},perframefn,parperframedata{flyi},...
                missingts{flyi},labelIdxVals{flyi},labelIdxImp{flyi},labelIdxT0{flyi}); %#ok<PFBNS>
              assert(tmpsuccess,'Loading windowdata failed for perframe function %s.',perframefn);
              
              if perframei==1
                parfor_predictblocks{perframei}{flyi}.t0 = predictblocks.t0;
                parfor_predictblocks{perframei}{flyi}.t1 = predictblocks.t1;
              end
              parfor_windowdata{perframei}{flyi}.X = windowdata.X;
              parfor_windowdata{perframei}{flyi}.t = windowdata.t;
              parfor_windowdata{perframei}{flyi}.labelidx_new = windowdata.labelidx_new;
              parfor_windowdata{perframei}{flyi}.labelidx_imp = windowdata.labelidx_imp;
            end
          end % parfor perframei
          
          for flyi = 1:Nfliescurr
            flies = flies_curr(flyi,:);
            
            tmp2 = length(parfor_predictblocks{1}{flyi}.t0);
            obj.predictblocks(iCls).expi = [obj.predictblocks(iCls).expi repmat(expi,1,tmp2)];
            obj.predictblocks(iCls).flies = [obj.predictblocks(iCls).flies repmat(flies,1,tmp2)];
            obj.predictblocks(iCls).t0 = [obj.predictblocks(iCls).t0 parfor_predictblocks{1}{flyi}.t0];
            obj.predictblocks(iCls).t1 = [obj.predictblocks(iCls).t1 parfor_predictblocks{1}{flyi}.t1];
            
            tmp = cellfun(@(x) x{flyi}, parfor_windowdata); % should be 1xNcpff struct array
            nframes = size(tmp(1).X,1);
            %ALTODO: How are we assured that the columns of [tmp.X]
            %correspond to the columns of windowdata.X? In fact
            %curperframefns may have mutated during this method...
            obj.windowdata(iCls).X = [obj.windowdata(iCls).X; [tmp.X]];
            obj.windowdata(iCls).exp = [obj.windowdata(iCls).exp; repmat(expi,nframes,1)];
            obj.windowdata(iCls).flies = [obj.windowdata(iCls).flies; repmat(flies,nframes,1)];
            obj.windowdata(iCls).t = [obj.windowdata(iCls).t; tmp(1).t];
            obj.windowdata(iCls).labelidx_cur = [obj.windowdata(iCls).labelidx_cur; zeros(nframes,1)];
            obj.windowdata(iCls).labelidx_new = [obj.windowdata(iCls).labelidx_new; tmp(1).labelidx_new];
            obj.windowdata(iCls).labelidx_imp = [obj.windowdata(iCls).labelidx_imp; tmp(1).labelidx_imp];
            obj.windowdata(iCls).labelidx_old = [obj.windowdata(iCls).labelidx_old; zeros(nframes,1)];
            obj.windowdata(iCls).scores_validated = [obj.windowdata(iCls).scores_validated; zeros(nframes,1)];
            if ~isempty(obj.windowdata(iCls).binVals),
              tmpbins = findThresholdBins([tmp.X], obj.windowdata(iCls).binVals);
            else
              tmpbins = [];
            end
            obj.windowdata(iCls).bins = [obj.windowdata(iCls).bins tmpbins];
%             obj.windowdata(iCls).predicted = [obj.windowdata(iCls).predicted; zeros(nframes,1)];
%             obj.windowdata(iCls).scores = [obj.windowdata(iCls).scores; zeros(nframes,1)];
%             obj.windowdata(iCls).scores_old = [obj.windowdata(iCls).scores_old; zeros(nframes,1)];
%             obj.windowdata(iCls).postprocessed = [obj.windowdata(iCls).postprocessed; zeros(nframes,1)];
%             obj.windowdata(iCls).isvalidprediction = [obj.windowdata(iCls).isvalidprediction; false(nframes,1)];
            %obj.windowdata(iCls).featurenames = [obj.windowdata(iCls).featurenames tmp(1).featurenames];
          end          
        end % iCls

        % AL: Move outside exp loop?
%         obj.windowdata = WindowData.windowdataTrim(obj.windowdata,...
%           @(x)x.labelidx_new==0);        
      end % expi
      
      % MK: Moved this out of exp loop as per AL's suggestion
      obj.windowdata = WindowData.windowdataTrim(obj.windowdata,...
        @(x)x.labelidx_new==0);        
      
      success = true;
      msg = '';
      
      obj.ClearStatus();
    end  % function/method
    
    
    % ---------------------------------------------------------------------
    function UpdateBoostingBins(obj,iCls)
      %MERGESTUPDATED
      islabeled = obj.windowdata(iCls).labelidx_new~=0;
      obj.windowdata(iCls).binVals = findThresholds(obj.windowdata(iCls).X(islabeled,:),...
        obj.classifier_params{iCls},'deterministic',obj.deterministic);
      obj.windowdata(iCls).bins = findThresholdBins(obj.windowdata(iCls).X(islabeled,:),obj.windowdata(iCls).binVals );
      
    end
    
    
    % Deprecated

    % ---------------------------------------------------------------------
    function MaybeStoreLabelsAndPreLoadWindowDataNow(self)
      % This is a hint to JLabelData that right now might be a good time to
      % write-back the currrent labels to the main store, and to pre-load the
      % window data.
      % This method is named as if it's a hint, but in the one place it's
      % currently called, it may well be required for proper behavior.
      % This method is deprecated b/c this the sort of book-keeping that
      % should be internal to JLabelData.  Callers shouldn't need to tell
      % JLabelData to get its house in order.
      self.StoreLabelsAndPreLoadWindowData();
    end  % method
    
  end
  
  methods (Static)
    
    function [X,feature_names] = ...
      ComputeWindowDataChunkStatic(curperframefns,allperframefns,perframefile,...
        flies,windowfeaturescellparams,t0,t1)
      % Computes a chunk of windowdata between frames t0 and t1 for flies
      % flies. Loop through all the per-frame features, and call
      % ComputeWindowFeatures to compute all the window data for that
      % per-frame feature.
      %
      % X: nframes x nfeatures window data
      % feature_names: 1 x nfeatures cell array labeling columns of X
      
      %MERGESTUPDATED
      
      assert(numel(allperframefns)==numel(perframefile));
      assert(isequal(fieldnames(windowfeaturescellparams),curperframefns(:)));
      
      X = [];
      feature_names = cell(1,numel(curperframefns));
      for j = 1:numel(curperframefns)
        fn = curperframefns{j};
        ndx = find(strcmp(fn,allperframefns));
        
        perframedata = load(perframefile{ndx}); %#ok
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
        if nargout>1
          feature_names{j} = cellfun(@(s) [{fn},s],feature_names_curr,'UniformOutput',false); 
        end
      end
      X = single(X);
      feature_names = [feature_names{:}];
    end
    
  end
  
  %% Target
  
  methods
    
     % ---------------------------------------------------------------------
    function val = IsCurFly(obj,expi,flies)
      val = all(flies == obj.flies) && (expi==obj.expi);
    end

    
    % ---------------------------------------------------------------------
    function expi = GetExp(obj)
      expi = obj.expi;
    end

    
    % ---------------------------------------------------------------------
    function flies = GetFlies(obj)
      flies = obj.flies;
    end
    
    
    % ---------------------------------------------------------------------
    function [success,msg] = setCurrentTarget(obj,expi,flies,headerinfo,force)
      % This is the method formerly known as PreLoad(). Sets the current
      % target to experiment expi, animal flies.  This implies preloading
      % data associated with the input experiment and flies. If neither the
      % experiment nor flies are changing, then we do nothing. If there is
      % currently a preloaded experiment, then we store the labels in
      % labelidx into labels using StoreLabels. We then load from labels into
      % labelidx for the new experiment and flies. We load the per-frame data
      % for this experiment and flies. If this is a different experiment,
      % then we load in the trajectories for this experiment.
      
      success = false;
      msg = '';
      
      if ~exist('force','var')
        force=false;
      end
      
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
      if ~diffflies && ~force,
        success = true;
        return;
      end

      if ~isempty(obj.expi) && obj.expi > 0,
        % store labels currently in labelidx to labels
        obj.StoreLabelsAndPreLoadWindowData();
      end
      
      if diffexpi || force,
        
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
          % Please DO NOT USE caching. Creates a lot of bugs if the user
          % doesn't know about it or forgets about it.
%           global CACHED_TRX; %#ok<TLEV>
%           global CACHED_TRX_EXPNAME; %#ok<TLEV>
%           if isempty(CACHED_TRX) || isempty(CACHED_TRX_EXPNAME) || ...
%               ~strcmp(obj.expnames{expi},CACHED_TRX_EXPNAME),
            [trx,~,~,~,trxrest] = load_tracks(trxfilename);
            ff = fieldnames(trx);
            for fnum = 1:numel(ff)
              if numel(trx(1).(ff{fnum})) == trx(1).nframes && ~strcmpi(ff{fnum},'sex');
                for fly = 1:numel(trx)
                  trx(fly).(ff{fnum}) = trx(fly).(ff{fnum})(:)';
                end
              end
            end
            
            if isfield(trxrest,'skeleton'),
              ntotal = numel(trxrest.skeleton) + sum(cellfun(@numel,trxrest.skeleton)) - 1;
              skeleton_edges = zeros(ntotal,1);
              off = 0;
              for i = 1:numel(trxrest.skeleton),
                skeleton_edges(off+1:off+numel(trxrest.skeleton{i})) = trxrest.skeleton{i};
                off = off + numel(trxrest.skeleton{i}) + 1;
              end
              for i = 1:numel(trx),
                trx(i).skeleton_edges = skeleton_edges;
              end
            end

            if obj.fromAPT
              trkfilename = obj.GetFile('trk',expi);
%               if any(cellfun(@(x) ~exist(x,'file'),trkfilename))
%                 msg = sprintf('APT Trk file %s does not exist',trkfilename);
%                 success = false;
%                 return;
%               end
              prev_width = 0;
              for ndx = 1:numel(trkfilename)
                %MK 20220629
                % not required to load trk info anymore. 
                % data is now stored in trx file in kpts

                % MK 20230303 - It seems we might have to add kpts for apt
                % projects with trx in case the trx file doens't have them
                if obj.aptInfo.has_trx && ~any(strcmp(fieldnames(trx),'kpts'))
                  trx = addAPTTrk2Trx(trx,trkfilename{ndx},'view',ndx,'aptInfo',...
                    obj.aptInfo,'prev_width',prev_width);
                end

                if exist('headerinfo','var') && iscell(headerinfo)
                  prev_width = prev_width + headerinfo{ndx}.nc;
                end
              end
              
            end

            obj.trx = trx;
%             CACHED_TRX = obj.trx;
%             CACHED_TRX_EXPNAME = obj.expnames{expi};
%           else
% %            fprintf('DEBUG: Using CACHED_TRX. REMOVE THIS\n');
%             obj.trx = CACHED_TRX;
%           end
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
 
      end  % if diffexpi

      % set labelidx from labels
      obj.SetStatus('Caching labels for experiment %s, flies%s',obj.expnames{expi},sprintf(' %d',flies));
      [obj.labelidx,obj.t0_curr,obj.t1_curr] = obj.GetLabelIdx(expi,flies);
      obj.labelidx_off = 1 - obj.t0_curr;
      
      % load perframedata
      obj.SetStatus('Loading per-frame data for %s, flies %s',obj.expdirs{expi},mat2str(flies));
      [success,msg]=obj.loadPerframeData(expi,flies);
      if ~success, return;  end
      
      obj.expi = expi;
      obj.flies = flies;

      if numel(obj.predictdata)<obj.expi
        obj.PredictDataInit(obj.expi);
      end
      obj.UpdatePredictedIdx();
      obj.ClearStatus();
           
      success = true;
    end
    
    
    % ---------------------------------------------------------------------
    function unsetCurrentTarget(obj)
      % Sets the object to a state where no target is currently selected.
      % This also clears the cached data for the currently loaded
      % experiment.
      obj.StoreLabelsForCurrentAnimal();
      obj.trx = {};
      obj.expi = 0;
      obj.flies = [];
      obj.perframedata = {};
      obj.labelidx = struct('vals',[],'imp',[],'timestamp',[]);
      obj.labelidx_off = 0;
      obj.t0_curr = 0;
      obj.t1_curr = 0;
      obj.predictedidx = [];
      obj.scoresidx = [];
      obj.scoresidx_old = [];
      %obj.erroridx = [];
    end    
             
  end
  
  methods (Access=public,Static=true)
    
    function valid = CheckExp(expi)
      if numel(expi) ~= 1,
        error('Usage: expi must be a scalar');
        %valid = false;
      else
        valid = true;
      end
    end
    
    function valid = CheckFlies(flies)
      if size(flies,1) ~= 1,
        error('Usage: one set of flies must be selected');
        %valid = false;
      else
        valid = true;
      end
    end
    
  end
  
  
  %% Labels
  
  methods
    
    
    % ---------------------------------------------------------------------
    function timestamp = GetLabelTimestamps(obj,expis,flies,ts)
      timestamp = nan(size(ts));      
      for expi = 1:obj.nexps,
        expidx = expis == expi;
        if ~any(expidx),
          continue;
        end
        % TODO: extend to multiple flies
        for fly = 1:obj.nflies_per_exp(expi),
          flyidx = expidx & flies == fly;
          if ~any(flyidx),
            continue;
          end
          [labelidx,T0] = obj.GetLabelIdx(expi,fly);
          timestamp(flyidx) = labelidx.timestamp(ts(flyidx)-T0+1);
        end
      end
    end
    
    
    % -----------------------------------------------------------
    function [labelidx,T0,T1] = GetLabelIdx(obj,expi,flies,T0,T1)
    % Returns the labelidx for the input experiment and flies read from
    % labels. 
    %
    % labelidx only guaranteed to have fields .vals, .imp, .timestamp.
    % labelidx.vals takes values in {0 (unlabeled),1,2,...,obj.nbehaviors}.
    % See Labels.labelIdx for more info.
    
      % Access labelidx cache if appropriate
      if ~isempty(obj.expi) && numel(flies)==numel(obj.flies) && obj.IsCurFly(expi,flies)
        if nargin < 4
          labelidx = obj.labelidx;
          T0 = obj.t0_curr;
          T1 = obj.t1_curr;
        else
          idx = T0+obj.labelidx_off:T1+obj.labelidx_off;
          labelidx = struct(); % cropped labelIdx, doesn't have all the usual labelidx fields
          labelidx.nbeh = obj.nbehaviors;
          labelidx.vals = obj.labelidx.vals(:,idx);
          labelidx.imp = obj.labelidx.imp(:,idx);
          labelidx.timestamp = obj.labelidx.timestamp(:,idx);
        end
        return;
      end
      
      if nargin < 4
        T0 = max(obj.GetTrxFirstFrame(expi,flies));
        T1 = min(obj.GetTrxEndFrame(expi,flies));
      end
      
      labelsShort = obj.GetLabels(expi,flies);
      %assert(obj.nbehaviors==numel(obj.labelnames));
      labelidx = Labels.labelIdx(obj.labelnames,T0,T1);
      labelidx = Labels.labelIdxInit(labelidx,labelsShort);
    end
    
    
    %----------------------------------------------------------------------
    function [labels,gtLabels] = getLabelsAndGTLabels(self)
      % Returns a single structure containing all the labels, suitable for
      % saving.
    
      % % short-circuit if no labels
      %if isempty(self.labels) && isempty(self.gt_labels), 
      %  labels=self.labels;
      %  gtLabels=self.gt_labels;
      %  return; 
      %end
          
      % Take the labels currently only in labelidx, and commit them to
      % self.labels
      self.StoreLabelsForCurrentAnimal();
      
      % Put the right labels in the right place
      if self.gtMode ,
        gtLabels=self.labels;
        labels=self.otherModeLabelsEtc.labels;
      else
        labels=self.labels;
        gtLabels=self.otherModeLabelsEtc.labels;
      end
        
    end
   
        
    % ---------------------------------------------------------------------
    function has = haslabels(obj,expnum)
      obj.StoreLabelsForCurrentAnimal();
      has = false;
      for fly = 1:numel(obj.labels(expnum).flies)
        if ~isempty(obj.labels(expnum).t0s{fly})
          has = true;
        end
      end
    end
    
    
    % ---------------------------------------------------------------------
    function tf = classifierHasLabels(obj,iClsVec)
      % iClsVec: vector of classifier indices
      % tf: logical vec, same size as iCls. tf(i) is true if there are any
      % labels (beh or no-beh) for obj.classifiernames{i}, across all exps.

      obj.StoreLabelsForCurrentAnimal();
      labels = obj.labels;
      iCls2LblNames = obj.iCls2LblNames;
      
      tf = false(size(iClsVec));
      for i = 1:numel(iClsVec)
        iCls = iClsVec(i);
        lbls = iCls2LblNames{iCls};
        tfTmp = Labels.labelsSeen(labels,lbls);
        tf(i) = any(tfTmp);
      end
    end
    
    
    % ---------------------------------------------------------------------
    function setAllLabels(self, ...
                          everythingParams)
      % MERGEST OK
      
      % AL 20141125: Unnecessary, SetExpDirs() call handles this
%       if self.nexps>0
%         self.RemoveExpDirs(1:self.nexps);
%       end

      if self.gtMode
        dirNames = everythingParams.gtExpDirNames;
        labels = everythingParams.gtLabels;
        self.otherModeLabelsEtc = struct('expDirNames',{everythingParams.expDirNames}, ...
                                       'labels',{everythingParams.labels});
      else
        % Normal labeling mode (not GT)
        dirNames = everythingParams.expDirNames;
        labels = everythingParams.labels;
        self.otherModeLabelsEtc = struct('expDirNames',{everythingParams.gtExpDirNames}, ...
                                       'labels',{everythingParams.gtLabels});
      end
      assert(numel(dirNames)==numel(labels));

      [success,msg] = self.SetExpDirs(dirNames);
      if ~success, error(msg); end
      
      [success,msg] = self.UpdateStatusTable();
      if ~success, error(msg); end

      self.setLabelsFromStructForAllExps(labels);
      
      % % set the GT labels
      % self.setGTLabelsFromStructForAllExps(everythingParams.gtLabels);

      % Preload the first track of the first video, which sets the current
      % experiment and track to experiment 1, track 1
      if (self.nexps>0)
        self.SetStatus('Pre-loading experiment %s...',self.expnames{1});
        [success1,msg1] = self.setCurrentTarget(1,1);
        if ~success1
          msg = sprintf('Error getting basic trx info: %s',msg1);
          self.SetStatus('Error getting basic trx info for %s.',self.expnames{1});
          uiwait(warndlg(msg));
          self.RemoveExpDirs(1:obj.nexps);
          self.ClearStatus();
          return;
        end
      end
                  
      % % update cached dataset
      % [success,msg] = self.PreLoadPeriLabelWindowData();  % need to move!!
      % if ~success,error(msg);end   
      
      % clear the cached per-frame, trx data
      self.ClearCachedPerExpData();
    end        
    
    
    % ---------------------------------------------------------------------
    function labels_curr = GetLabels(obj,expi,flies,tfSkipUpdateFromCache)
      % Returns the labels for the given experiment index, fly(s) index.
      % This takes into account the current GT mode (normal/GT), and
      % returns the appropriate labels.
      %
      % tfSkipUpdateFromCache: optional scalar logical, defaults to false.
      % If true, do not update labels from .labelIdx for current target 
      % (when relevant).
      %
      % labels_curr: a 'labelsShort', see Labels.m
      

      % MERGESTUPDATED
      
      labels_curr = Labels.labelsShort();
              
      if nargin < 2 || isempty(expi)
        expi = obj.expi;
      end
      if expi < 1
        % AL: questionable legacy codepath
        return;
      end
      if nargin < 3 || isempty(flies)
        flies = obj.flies;
      end
      if nargin < 4
        tfSkipUpdateFromCache = false;
      else
        assert(islogical(tfSkipUpdateFromCache) && isscalar(tfSkipUpdateFromCache));
      end

      % cache these labels if current experiment and flies selected
      if expi==obj.expi && all(flies==obj.flies) && ~tfSkipUpdateFromCache
        obj.StoreLabelsForCurrentAnimal();
      end
      
      [labels_curr,tffly] = Labels.labelsShortInit(labels_curr,obj.labels(expi),flies);
      if ~tffly
        t0_curr = max(obj.GetTrxFirstFrame(expi,flies));
        labels_curr.off = 1-t0_curr;
      end
    end
    
      
    % ---------------------------------------------------------------------
    function SetLabel(obj,expi,flies,ts,behaviori,important)
    % SetLabel(obj,expi,flies,ts,behaviori,important)
    % Set behavior label for experiment expi, flies, and frames ts. If
    % expi, flies match current expi, flies, then we only set labelidx.
    % Otherwise, we set labels.
    %
    % expi: scalar
    % flies: scalar
    % ts: vector 
    % behaviori: scalar, same as JLabel/SetLabelPlot:
    % - If greater than 0, behavior index to label. Valid vals: [1..numbehaviors]
    % - If == 0, clear labels for all timelines/behaviors.
    % - If < 0, timeline index to clear. Valid vals: [1..numtimelines]
    % important: scalar
    
      assert(isscalar(behaviori));
      
      if behaviori<=0 % clear labels
        assert(important==0);        
        if behaviori==0
          iTL = 1:obj.ntimelines;
        else
          iTL = -behaviori;
        end  
        lblVal = 0;
      else % label
        iTL = obj.labelidx.idxBeh2idxTL(behaviori);
        lblVal = behaviori;
      end
     
      obj.setLabelsRaw(expi,flies,ts,iTL,lblVal,important);      
      %obj.UpdateErrorIdx();
      obj.needsave = true;
    end
    
    
    function SetLabelsToMatchPrediction(obj,expi,flies,t0,t1,important)
      % Set the labels (currently, for all classifiers/timelines) to match
      % the current predictions.
      %
      % important: scalar
      
      assert(isscalar(important));
      
      [prediction,t0new,t1new] = obj.GetPredictedIdx(expi,flies,t0,t1);
      if t0new~=t0 || t1new~=t1
        warningNoTrace('JLabelData:modifiedTime','Using modified time interval.');
        t0 = t0new;
        t1 = t1new;
      end
      predictedidx = obj.PredictedIdxExpandValues(prediction.predictedidx);
            
      obj.setLabelsRaw(expi,flies,t0:t1,1:obj.ntimelines,predictedidx,important);
      obj.needsave = true;
    end
    
    
    % ---------------------------------------------------------------------
    function tf = getAtLeastOneNormalLabelOfEachClassExists(self,labelNames)
      % Returns true iff at least one normal (non-GT) label exists for each
      % member of labelNames. 
      %
      % labelNames: cellstr, subset of self.labelnames. Defaults to
      % self.labelnames
      %
      % tf: scalar logical

      %MERGESTUPDATED
      
      if ~exist('labelNames','var') || isempty(labelNames)
        labelNames = self.labelnames;
      else
        assert(iscellstr(labelNames) && all(ismember(labelNames,self.labelnames)));
      end
      
      tf = Labels.labelsSeen(self.labels,labelNames);      
      if ~all(tf)
        % Haven't seen all labels yet, check labelidx
        iLbl = unique(self.labelidx.vals(:)); % indices into self.labelnames (or 0)
        iLbl = iLbl(iLbl>0);
        names = self.labelnames(iLbl); % all unique names seen in .labelidx      
        tf = tf | ismember(labelNames,names);
      end
      
      % MAYANK_JAN16_2015: For multiclassifier, 
      % this should return true if labels for behavior and not behavior 
      % exist for a single behavior. Should require labels for all
      % behaviors. Below assumes that behaviors are in the first half and 
      % not behaviors are in the second half.
      tf = tf(1:end/2) & tf(end/2+1:end);
      
      tf = any(tf);      
    end
        
    
%     % ---------------------------------------------------------------------
%     function atLeastOneNormalLabelExists=getAtLeastOneNormalLabelExists(self)
%       % Returns true iff at least one normal (non-GT) label exists.
%       atLeastOneNormalLabelExists=false;
%       for i=1:self.nexps
%         if ~isempty(self.labels(i).t0s)
%           atLeastOneNormalLabelExists=true;
%           break
%         end
%       end
%       % If we haven't found any labels yet, check labelidx
%       if ~atLeastOneNormalLabelExists,
%         if any(self.labelidx.vals),
%           atLeastOneNormalLabelExists=true;
%         end
%       end
%     end
    % ---------------------------------------------------------------------
%     function SetLabel_KB(obj,expi,flies,ts,behaviori,important)
% 
%       % SetLabel_KB(obj,expi,flies,ts,behaviori)
%       % Set label for experiment expi, flies, and frames ts to behaviori.
%       % Store everywhere.
%       
%       % store in labelidx
%       if obj.IsCurFly(expi,flies),
%         obj.labelidx.vals(ts+obj.labelidx_off) = behaviori;
%         obj.labelidx.imp(ts+obj.labelidx_off) = important;
%         obj.labelidx.timestamp(ts+obj.labelidx_off) = now;
%       end      
% 
%       % store in labels
%       [labelidx,T0] = obj.GetLabelIdx(expi,flies);
%       labelidx.vals(ts+1-T0) = behaviori;
%       labelidx.imp(ts+1-T0) = important;
%       labelidx.timestamp(ts+1-T0) = now;
%       obj.StoreLabels1(expi,flies,labelidx,1-T0);
%       
%       % store in windowdata
%       idxcurr = obj.windowdata.exp == expi & obj.windowdata.flies == flies & ismember(obj.windowdata.t,ts);
%       obj.windowdata.labelidx_new(idxcurr) = behaviori;
%       obj.windowdata.labelidx_imp(idxcurr) = important;
%       
%     end

        % ---------------------------------------------------------------------
%     function isstart = IsLabelStart(obj,expi,flies,ts)
%       % AL: appears to be useless
%       
%       if obj.expi == expi && all(flies == obj.flies),
%         isstart = obj.labelidx.vals(ts+obj.labelidx_off) ~= 0 & ...
%           obj.labelidx.vals(ts+obj.labelidx_off-1) ~= obj.labelidx.vals(ts+obj.labelidx_off);
%       else
%         
% %         if obj.IsGTMode(),
% %           labelsToUse = 'gt_labels';
% %         else
% %           labelsToUse = 'labels';
% %         end
% 
%         if isempty(obj.labels(expi).flies)
%           ism = false;
%         else
%           [ism,fliesi] = ismember(flies,obj.labels(expi).flies,'rows');
%         end
% 
%         if ism,
%           isstart = ismember(ts,obj.labels(expi).t0s{fliesi});
%         else
%           isstart = false(size(ts));
%         end
%       end
%       
%         end

    
%     % ---------------------------------------------------------------------
%     function ClearLabels(obj,expi,flies)
%       
%       
%       if obj.nexps == 0,
%         return;
%       end
%       
%       
%       %timestamp = now;
%       
%       % use all experiments by default
%       if nargin < 2,
%         expi = 1:obj.nexps;
%       end
%       
%       % delete all flies by default
%       if nargin < 3,
%         for i = expi(:)',
%           obj.labels(expi).t0s = {};
%           obj.labels(expi).t1s = {};
%           obj.labels(expi).names = {};
%           obj.labels(expi).flies = [];
%           obj.labels(expi).off = [];
%           obj.labels(expi).timestamp = {};
%           obj.labelstats(expi).nflies_labeled = 0;
%           obj.labelstats(expi).nbouts_labeled = 0;
%           obj.labels(expi).imp_t0s = {};
%           obj.labels(expi).imp_t1s = {};
%         end
%       else
%         if numel(expi) > 1,
%           error('If flies input to ClearLabels, expi must be a single experiment');
%         end
%         % no labels
%         if numel(obj.labels) < expi,
%           return;
%         end
%         % which index of labels
%         [~,flyis] = ismember(obj.labels(expi).flies,flies,'rows');
%         for flyi = flyis(:)',
%           % keep track of number of bouts so that we can update stats
%           ncurr = numel(obj.labels(expi).t0s{flyi});
%           obj.labels(expi).t0s{flyi} = [];
%           obj.labels(expi).t1s{flyi} = [];
%           obj.labels(expi).names{flyi} = {};
%           obj.labels(expi).timestamp{flyi} = [];
%           obj.labels(expi).imp_t0s{flyi} = [];
%           obj.labels(expi).imp_t1s{flyi} = [];
%           % update stats
%           obj.labelstats(expi).nflies_labeled = obj.labelstats(expi).nflies_labeled - 1;
%           obj.labelstats(expi).nbouts_labeled = obj.labelstats(expi).nbouts_labeled - ncurr;
%         end
%       end
%       
%       % clear labelidx if nec
%       if ismember(obj.expi,expi) && ((nargin < 3) || ismember(obj.flies,flies,'rows')),
%         obj.labelidx.vals(:) = 0;
%         obj.labelidx.imp(:) = 0;
%         obj.labelidx.timestamp(:) = 0;
%       end
%       
%       % clear windowdata labelidx_new
%       for i = expi(:)',
%         if nargin < 3,
%           idx = obj.windowdata.exp == i;
%         else
%           idx = obj.windowdata.exp == i & ismember(obj.windowdata.flies,flies,'rows');
%         end
%         obj.windowdata.labelidx_new(idx) = 0;
%         obj.windowdata.labelidx_imp(idx) = 0;
%         obj.UpdateErrorIdx();
%       end
%       
%     end

%     ---------------------------------------------------------------------
%     function [idx,T0,T1] = IsBehavior(obj,behaviori,expi,flies,T0,T1)
%     [idx,T0,T1] = IsBehavior(obj,behaviori,expi,flies,T0,T1)
%     Returns whether the behavior is labeled as behaviori for experiment
%     expi, flies from frames T0 to T1. If T0 and T1 are not input, then
%     firstframe to endframe are used. 
% 
%           
%       if ~isempty(obj.expi) && expi == obj.expi && numel(flies) == numel(obj.flies) && all(flies == obj.flies),
%         if nargin < 4,
%           idx = obj.labelidx.vals == behaviori;
%           T0 = obj.t0_curr;
%           T1 = obj.t1_curr;
%         else
%           idx = obj.labelidx.vals(T0+obj.labelidx_off:T1+obj.labelidx_off) == behaviori;
%         end
%         return;
%       end
%       
%       if nargin < 4,
%         T0 = max(obj.GetTrxFirstFrame(expi,flies));
%         T1 = min(obj.GetTrxEndFrame(expi,flies));
%       end
%       n = T1-T0+1;
%       off = 1 - T0;
%       labels_curr = obj.GetLabels(expi,flies);
%       idx = false(1,n);
%       for j = find(strcmp(labels_curr.names,obj.labelnames{behaviori})),
%         t0 = labels_curr.t0s(j);
%         t1 = labels_curr.t1s(j);
%         idx(t0+off:t1-1+off) = true;
%       end
%       
%     end

  end
  
  methods % more private
    
    
    function setLabelsRaw(obj,expi,flies,ts,iTL,lblVal,imp)
      % ts: vector of time indices
      % iTL: vector of timelines
      % lblVal: either a scalar, or array of size numel(iTL)-by-numel(ts)
      % imp: either a scalar, or array of size numel(iTL)-by-numel(ts)
      
      if obj.IsCurFly(expi,flies)
        obj.labelidx.vals(iTL,ts+obj.labelidx_off) = lblVal;
        obj.labelidx.imp(iTL,ts+obj.labelidx_off) = imp;
        obj.labelidx.timestamp(iTL,ts+obj.labelidx_off) = now;
      else
        [labelidx,T0] = obj.GetLabelIdx(expi,flies); %#ok<*PROP>
        labelidx.vals(iTL,ts+1-T0) = lblVal;
        labelidx.imp(iTL,ts+1-T0) = imp;
        labelidx.timestamp(iTL,ts+1-T0) = now;
        obj.StoreLabelsForGivenAnimal(expi,flies,labelidx,1-T0);
      end
    end
    
    
    % ---------------------------------------------------------------------
    function setLabelsFromStructForAllExps(self,labelsForAll)
      statusTableString = fif(self.gtMode,'gt_label','label');
      nExps = length(self.expdirs);
      for expi = 1:nExps
        self.loadLabelsFromStructForOneExp(expi,labelsForAll(expi));
        self.UpdateStatusTable(statusTableString);   
      end
    end
    
    
    % ---------------------------------------------------------------------
    function loadLabelsFromStructForOneExp(self,expi,labels)
      % Load the labels for a single experiment into self.
            
      self.SetStatus('Loading labels for %s',self.expdirs{expi});

      self.labels(expi).t0s = labels.t0s;
      self.labels(expi).t1s = labels.t1s;
      self.labels(expi).names = labels.names;
      self.labels(expi).flies = labels.flies;
      self.labels(expi).off = labels.off;
%       self.labelstats(expi).nflies_labeled = size(labels.flies,1);
%       self.labelstats(expi).nbouts_labeled = numel([labels.t0s{:}]);
      Nfly = numel(labels.flies);
      if iscell(labels.timestamp)
        self.labels(expi).timestamp = labels.timestamp;
      else
        for ndx = 1:Nfly
          nBouts = numel(labels.t0s{ndx});
          if isempty(labels.timestamp)
            self.labels(expi).timestamp{ndx}(1:nBouts) = now;
          else
            self.labels(expi).timestamp{ndx}(1:nBouts) = labels.timestamp;
          end
        end
      end
      if isfield(labels,'timelinetimestamp')
        timelineTS = labels.timelinetimestamp;
        assert(iscell(timelineTS) && isequal(size(timelineTS),[1 Nfly]));
        self.labels(expi).timelinetimestamp = timelineTS;
      else
        self.labels(expi).timelinetimestamp = cell(1,Nfly);
        for ndx = 1:Nfly
          self.labels(expi).timelinetimestamp{ndx} = struct();
        end
      end
      if isfield(labels,'imp_t0s');
        self.labels(expi).imp_t0s = labels.imp_t0s;
        self.labels(expi).imp_t1s = labels.imp_t1s;
      else
        self.labels(expi).imp_t0s = cell(1,numel(labels.flies));
        self.labels(expi).imp_t1s = cell(1,numel(labels.flies));
      end

      self.ClearStatus();
    end  % method
   
    
    % ---------------------------------------------------------------------
    function StoreLabelsAndPreLoadWindowData(obj)
      % Store labels cached in labelidx for the current experiment and flies
      % to labels structure. This is when the timestamp on labels gets
      % updated.  Also preloads the window data if not in GT mode.
      
      % flies not yet initialized
      if isempty(obj.flies) || all(isnan(obj.flies)) || isempty(obj.labelidx.vals),
        return;
      end
      
      obj.StoreLabelsForCurrentAnimal();
            
%       % preload labeled window data while we have the per-frame data loaded
%       ts = find(obj.labelidx.vals~=0) - obj.labelidx_off;
%       if ~obj.IsGTMode(),
%         [success,msg] = obj.PreLoadWindowData(obj.expi,obj.flies,ts);
%         if ~success,
%           warning(msg);
%         end
%       end
%       
%       % update windowdata's labelidx_new
%       if ~isempty(obj.windowdata.exp),
%         idxcurr = obj.windowdata.exp == obj.expi & ...
%           all(bsxfun(@eq,obj.windowdata.flies,obj.flies),2);
%         obj.windowdata.labelidx_new(idxcurr) = obj.labelidx.vals(obj.windowdata.t(idxcurr)+obj.labelidx_off);
%         obj.windowdata.labelidx_imp(idxcurr) = obj.labelidx.imp(obj.windowdata.t(idxcurr)+obj.labelidx_off);
%       end
      
      %obj.UpdateWindowDataLabeled(obj.expi,obj.flies);
      
    end  % method
  
    
    % ---------------------------------------------------------------------
    function StoreLabelsForCurrentAnimal(obj)
      % Store labels cached in labelidx for the current experiment and flies
      % to labels structure. This is when the timestamp on labels gets
      % updated. 
      if isempty(obj.flies) || obj.expi==0 || all(isnan(obj.flies)) || isempty(obj.labelidx.vals),
        return
      end      
      obj.StoreLabelsForGivenAnimal(obj.expi,obj.flies,obj.labelidx,obj.labelidx_off);
    end  % method

    
    % ---------------------------------------------------------------------
    function StoreLabelsForGivenAnimal(obj,expi,flies,labelidx,labelidx_off)
      
      % MERGESTUPDATED
      
      assert(labelidx.off==labelidx_off);
      assert(isequal(labelidx.labelnames,obj.labelnames));
      assert(labelidx.nbeh==obj.nbehaviors);      
      
      labelsShort = Labels.labelsShortFromLabelIdx(labelidx);
      obj.labels(expi) = Labels.assignFlyLabelsRaw(obj.labels(expi),labelsShort,flies);

      % Update labels.timelinetimestamp
      maxtimestamps = max(labelidx.timestamp,[],2); % most recent timestamp each timeline was edited
      [tf,ifly] = ismember(flies,obj.labels(expi).flies,'rows');
      assert(tf);
      for iTL = 1:labelidx.nTL
        classifiername = obj.labelnames{iTL};
        if ~isfield(obj.labels(expi).timelinetimestamp{ifly},classifiername)
          obj.labels(expi).timelinetimestamp{ifly}.(classifiername) = 0;
        end
        obj.labels(expi).timelinetimestamp{ifly}.(classifiername) = ...
          max(obj.labels(expi).timelinetimestamp{ifly}.(classifiername),maxtimestamps(iTL));
      end
    end

    
    % ---------------------------------------------------------------------
    function bouts = getLabeledBouts(obj,iCls)
    % Find window data for labeled bouts.
    %
    % bouts
    % .ndx: nBout-by-nsamp logical. bouts.ndx(iBout,:) indexes obj.windowdata(iCls).t etc
    % .label: 1-by-nBout vector of 1/2 for positive/negative lbls for this classifier
    % .timestamp: 1-by-nBout vector
        
    %MERGESTUPDATED
    
      bouts = struct('ndx',[],'label',[],'timestamp',[]);
      
      wd = obj.windowdata(iCls);
      clsLblNames = obj.iCls2LblNames{iCls};
      % clsLblNames = {posLbl negLbl}
      
      for expNdx = 1:obj.nexps
        for flyNdx = 1:obj.nflies_per_exp(expNdx)
          
          labelsShort = obj.GetLabels(expNdx,flyNdx);
          nBout = numel(labelsShort.t0s);
          tfFlyNdx = obj.FlyNdx(expNdx,flyNdx,iCls);
          for iBout = 1:nBout
            tfClsLbl = strcmp(labelsShort.names{iBout},clsLblNames);
            if any(tfClsLbl)
              idx = tfFlyNdx & ...
                wd.t >= labelsShort.t0s(iBout) & ...
                wd.t < labelsShort.t1s(iBout) & ...
                wd.labelidx_imp;
              if ~all(wd.labelidx_new(idx))
                continue;
              end
              bouts.ndx(end+1,:) = idx;
              bouts.label(1,end+1) = find(tfClsLbl); % 1/2 for posLbl/negLbl resp
              bouts.timestamp(1,end+1) = labelsShort.timestamp(iBout);
            end
          end
        end
      end      
    end

    
    function reinitLabelIdx(obj)
      % reinitialize obj.labelIdx from obj.labels for the current fly/exp.
      % obj.labels must be up-to-date.
      
      expi = obj.expi;
      flies = obj.flies;
      T0 = max(obj.GetTrxFirstFrame(expi,flies));
      T1 = min(obj.GetTrxEndFrame(expi,flies));
      
      labelsShort = obj.GetLabels(expi,flies,true);
      lblIdx = Labels.labelIdx(obj.labelnames,T0,T1);
      obj.labelidx = Labels.labelIdxInit(lblIdx,labelsShort);
    end
    
    
    
%     % ---------------------------------------------------------------------
%     function bouts = GetLabeledBouts_KB(obj)
%       
%       bouts = struct('t0s',[],'t1s',[],'flies',[],'expis',[],'timestamps',[],'names',{{}});
%       for expi = 1:numel(obj.labels),
%         for flyi = 1:size(obj.labels(expi).flies,1),
%           flies = obj.labels(expi).flies(flyi,:);
%           t0s = obj.labels(expi).t0s{flyi};
%           t1s = obj.labels(expi).t1s{flyi};
%           if isempty(t0s),
%             continue;
%           end
%           n = numel(t0s);
%           bouts.t0s(end+1:end+n) = t0s;
%           bouts.t1s(end+1:end+n) = t1s;
%           bouts.flies(end+1:end+n,:) = flies;
%           bouts.expis(end+1:end+n) = expi;
%           bouts.timestamps(end+1:end+n) = obj.labels(expi).timestamp{flyi};
%           bouts.names(end+1:end+n) = obj.labels(expi).names{flyi};
%         end
%       end
%       
%     end  % method
  end  
  
  %% Train/Predict
  
  methods % (more public)

    
    % ---------------------------------------------------------------------
    function Train(obj)
    % Train(obj)
    % Updates the classifier to reflect the current labels. This involves
    % first loading/precomputing the training features. Then, the clasifier
    % is trained/updated. Finally, predictions for the currently loaded
    % window data are updated. Currently, the only implemented classifier is 
    % boosting. If the classifier exists, then it is updated instead of
    % retrained from scratch. This involves three steps -- replacing labels
    % for frames which have changed label, removing examples for frames
    % which have been removed the training set, and adding new examples for
    % newly labeled frames. If the classifier has not yet been trained, it
    % is trained from scratch. 
      
    % MERGEST UPDATED
    
      assert(~obj.isST);
    
      obj.StoreLabelsAndPreLoadWindowData();
      
      % load all labeled data
      [success,msg] = obj.PreLoadPeriLabelWindowData();
      if ~success
        error('JLabelData:unableToLoadPerLabelWindowData',msg);
      end
      
      cls2IdxBeh = obj.iCls2iLbl;
      assert(numel(cls2IdxBeh)==obj.nclassifiers);
      assert(iscell(obj.trainstats) && numel(obj.trainstats)==obj.nclassifiers);
      didwarnopt = false;
      for iCls = 1:obj.nclassifiers
        islabeled = obj.windowdata(iCls).labelidx_new~=0 & obj.windowdata(iCls).labelidx_imp;
        if ~any(islabeled) && ~obj.selFeatures(iCls).do
          continue;
        end

        switch obj.classifiertype{iCls}
          case {'boosting','fancyboosting'}
            %fprintf('!!REMOVE THIS: resetting the random number generator for repeatability!!\n');
            %stream = RandStream.getGlobalStream;
            %reset(stream);

            if true % obj.DoFullTraining
              if obj.selFeatures(iCls).use 
                if obj.selFeatures(iCls).do,
                  pstr = sprintf('Optimizing window features used for %s classifier from %d examples ...',obj.labelnames{iCls},nnz(islabeled)); 
                else
                  pstr = sprintf('Training optimized %s classifier from %d examples ...',obj.labelnames{iCls},nnz(islabeled)); 
                end
              else
                pstr = sprintf('Training %s classifier from %d examples...',obj.labelnames{iCls},nnz(islabeled)); 
              end
              obj.SetStatus(pstr);

              % form label vec that has 1 for 'behavior present'
              labelidxnew = obj.windowdata(iCls).labelidx_new(islabeled);
              assert(numel(cls2IdxBeh{iCls})==2);
              valBeh = cls2IdxBeh{iCls}(1);
              valNoBeh = cls2IdxBeh{iCls}(2);
              assert(all(labelidxnew==valBeh | labelidxnew==valNoBeh));
              labels12 = 2*ones(size(labelidxnew));
              labels12(labelidxnew==valBeh) = 1;
              % check for presence of both positive and negative labels
              npos = nnz(labels12==1);
              nneg = nnz(labels12~=1);
              if npos < 1 || nneg < 1
                warnstr = sprintf('Classifier %s: Only behavior or nones have been labeled. Not training classifier.',...
                  obj.labelnames{iCls});
                uiwait(warndlg(warnstr));
                continue;
              end
              
              obj.classifier_old{iCls} = obj.classifier{iCls};
              if checkThresholds(obj.windowdata(iCls).X(islabeled,:),...
                  obj.classifier_params{iCls},obj.windowdata(iCls).binVals),
                [obj.windowdata(iCls).binVals] = findThresholds(...
                  obj.windowdata(iCls).X(islabeled,:),...
                  obj.classifier_params{iCls},'deterministic',obj.deterministic);
                obj.windowdata(iCls).bins = findThresholdBins(obj.windowdata(iCls).X,...
                  obj.windowdata(iCls).binVals);
              end
              
              %fprintf('!!REMOVE THIS: resetting the random number generator for repeatability!!\n');
              %stream = RandStream.getGlobalStream;
              %reset(stream);
              assert(~isempty(obj.windowdata(iCls).bins));
              
              % Do feature selection.
             if obj.selFeatures(iCls).use && obj.selFeatures(iCls).do,
               bins = obj.windowdata(iCls).bins(:,islabeled);
                obj.selFeatures(iCls) =...
                  SelFeatures.select(obj.selFeatures(iCls),obj.windowdata(iCls).X(islabeled,:), ...
                                  labels12,obj,...
                                  obj.windowdata(iCls).binVals,...
                                  bins, ...
                                  obj.classifier_params{iCls},pstr);
                 
             end
             
             WindowData.windowdataVerify(obj.windowdata(iCls));

             if strcmp(obj.classifiertype,'boosting') 
                 if obj.selFeatures(iCls).use,
                   if ~didwarnopt % Warn about the optimization only once.
                    [obj.selFeatures(iCls),didwarnopt] = SelFeatures.checkForOpt(obj.selFeatures(iCls),labels12);
                   end

                   curD = obj.windowdata(iCls).X(islabeled,obj.selFeatures(iCls).f);
                   binVals = obj.windowdata(iCls).binVals(:,obj.selFeatures(iCls).f);
                   bins = obj.windowdata(iCls).bins(obj.selFeatures(iCls).f,islabeled);
                   [curcls,outscores,trainstats] =...
                     boostingWrapper(curD, ...
                     labels12,obj,...
                     binVals,...
                     bins, ...
                     obj.classifier_params{iCls},pstr);
                   obj.classifier{iCls} = SelFeatures.convertClassifier(curcls,obj.selFeatures(iCls));
                 else
                   curD = obj.windowdata(iCls).X(islabeled,:);
                  bins = obj.windowdata(iCls).bins(:,islabeled);
                  [obj.classifier{iCls},outscores,trainstats] =...
                    boostingWrapper(curD, ...
                                    labels12,obj,...
                                    obj.windowdata(iCls).binVals,...
                                    bins, ...
                                    obj.classifier_params{iCls},pstr);
                 end
                 perr = nnz( outscores<0 & labels12==1);
                 nerr = nnz( outscores>0 & labels12~=1);
                 if ((perr+nerr)>0) 
                   obj.trainWarnCount = obj.trainWarnCount + 1;
                   if obj.trainWarn && obj.trainWarnCount>5,
                     wstr = {sprintf('The current trained classifier made %d errors on the %d labeled frames',perr+nerr,nnz(labels12)),...
                             'Consider increasing the number of training iteratoins,',...
                             '(Using Iterations field from Menu -> Classifier -> Training Parameters).',...
                             'You might also want to review labels for inconsistencies.'
                             };
                     warndlg(wstr,'Increase training iterations');
                     obj.trainWarnCount = 0;
                   end
                 else % reset if there is no training error..
                   obj.trainWarnCount = 0;
                 end
                 
              else
                [obj.classifier{iCls}] = ...
                  fastBag(obj.windowdata(iCls).X(islabeled,:),...
                          labels12,...
                          obj.windowdata(iCls).binVals,...
                          obj.windowdata.bins(:,islabeled), ...
                          obj.classifier_params{iCls});
                trainstats = struct;
              end
              obj.lastFullClassifierTrainingSize = nnz(islabeled);

            else
              assert(false,'MERGE unreachable codepath not updated.');
%               oldNumPts = nnz(obj.windowdata.labelidx_cur ~= 0 & obj.windowdata.labelidx_imp );
%               newNumPts = nnz(obj.windowdata.labelidx_new ~= 0 & obj.windowdata.labelidx_imp );
%               newData = newNumPts - oldNumPts;
%               obj.SetStatus('Updating boosting classifier with %d examples...',newData);
% 
%               bins = findThresholdBins(obj.windowdata.X(islabeled,:),obj.windowdata.binVals);
% 
%               obj.classifier_old = obj.classifier;
%               if strcmp(obj.classifiertype,'boosting'),
%                 [obj.classifier, ~, trainstats] = ...
%                   boostingUpdate(obj.windowdata.X(islabeled,:),...
%                                  obj.windowdata.labelidx_new(islabeled),...
%                                  obj.classifier, ...
%                                  obj.windowdata.binVals,...
%                                  bins, ...
%                                  obj.classifier_params);
%               else
%                 error('Fast updates not defined for fancy boosting.');
%               end
            end
            
            obj.classifierTS(iCls) = now();

            % store training statistics
            if isempty(obj.trainstats{iCls})
              obj.trainstats{iCls} = trainstats;
              obj.trainstats{iCls}.timestamps = obj.classifierTS(iCls);
            else
              assert(isstruct(obj.trainstats{iCls}));
              trainstatfns = fieldnames(trainstats);
              for fn = trainstatfns(:)',fn=fn{1}; %#ok<FXSET>
                val = trainstats.(fn);
                if isempty(val)
                  val = nan;
                end
                obj.trainstats{iCls}.(fn)(end+1) = val;
              end
              obj.trainstats{iCls}.timestamps(end+1) = obj.classifierTS(iCls);
            end

            obj.windowdata(iCls).labelidx_old = obj.windowdata(iCls).labelidx_cur;
            obj.windowdata(iCls).labelidx_cur = obj.windowdata(iCls).labelidx_new;
            obj.windowdata(iCls).scoreNorm = nan;
            
            % To later find out where each example came from.
            
  %           obj.windowdata.isvalidprediction = false(numel(islabeled),1);
        end
      end
      
      obj.PredictDataMoveCurToOld();
      obj.FindFastPredictParams();
      obj.PredictLoaded(obj.predictOnlyCurrentFly);
 
      obj.needsave = true;
      obj.ClearStatus();      
    end 
    
    
    function TrainST(obj)
      assert(obj.isST);
      
      % Right now working off the jab, presumably better to just use obj
      Q = loadAnonymous(obj.everythingFileNameAbs);
      
      % Figure out which classifiers are out of date
      timelineTSes = Labels.mostRecentTimelineTimestamps(Q.labels);
      classifierTSes = getstructarrayfield(Q.classifierStuff,'timeStamp',...
        'numericscalar',true);
      
      [realbehnames,nobehnames] = Labels.verifyBehaviorNames(Q.behaviors.names);
      Nrealbeh = numel(realbehnames);
      assert(all(ismember(fieldnames(timelineTSes),realbehnames)));
      assert(Nrealbeh==numel(classifierTSes));
      
      deltaTS = nan(Nrealbeh,1);
      for iBeh = 1:Nrealbeh
        bname = realbehnames{iBeh};
        if isfield(timelineTSes,bname)
          lblTS = timelineTSes.(bname);
        else
          lblTS = inf;
        end
        clsTS = classifierTSes(iBeh);
        deltaTS(iBeh) = lblTS-clsTS;
      end
      trainTS = deltaTS>0;
      
      % Train
      if isfield(Q.extra,'usePastOnly'),
        usePastOnly = Q.extra.usePastOnly;
      else
        usePastOnly = false;
      end
      trainIdx = find(trainTS);
      trainIdx = trainIdx(:)';
      if ~isempty(trainIdx)
        for tIdx = trainIdx
          fprintf('Training classifier ''%s'', %s out of date.\n',...
            realbehnames{tIdx},deltadatenumstr(deltaTS(tIdx)));
        end
      else
        fprintf('All classifiers up-to-date, no training necessary.\n');
      end
      pause(2.0);
      for tIdx = trainIdx
        obj.classifier{tIdx} = trainDetectorSTCore(Q.expDirNames,Q.labels,realbehnames{tIdx},nobehnames{tIdx},usePastOnly);
        obj.classifierTS(tIdx) = now;
      end
      
      obj.needsave = true;
    end
      
    
    % ---------------------------------------------------------------------
    function Predict(obj,expi,flies,t0,t1)
    % Predict(obj,expi,flies,ts)
    % Runs the behavior classifier on the input experiment, flies, and
    % frames. This involves first precomputing the window data for these
    % frames, then applying the classifier. 
     
      %MERGESTUPDATED
      
      assert(~obj.isST,'Time-restricted prediction unsupported for multi-classifier projects.');
      
      % TODO: don't store window data just because predicting. 

      if isempty(t0) || t0>t1
        return;
      end
            
      % Get the relevant window data in main memory so that we can predict on it
      % [success,msg] = obj.PreLoadWindowData(expi,flies,t0:t1);
      % if ~success ,
      % error('JLabelData:unableToLoadWindowDataForPrediction',msg);
      % end
            
      didpredict = false(obj.nclassifiers,1);
      for iCls = 1:obj.nclassifiers
        if isempty(obj.classifier{iCls})
          continue;
        end
        
        switch obj.classifiertype{iCls}
          case {'boosting','fancyboosting'}
            if obj.nclassifiers==1
              obj.SetStatus('Updating Predictions ...');
            else
              obj.SetStatus('Updating Predictions for %s...',obj.labelnames{iCls});
            end
            didpredict(iCls) = obj.PredictFast(expi,flies,t0,t1,iCls);
            if didpredict(iCls)
              obj.SetStatus('Predictions updated, applying postprocessing...\n');
            end
        end        
      end
      
      if any(didpredict)
        obj.ApplyPostprocessing(expi,flies);
        obj.UpdatePredictedIdx();
      end
      obj.ClearStatus();
    end
      
    
    % ---------------------------------------------------------------------
    function predictForCurrentTargetAndTimeSpan(obj)
      obj.Predict(obj.expi,obj.flies,obj.t0_curr,obj.t1_curr)
    end
    
    
    % ---------------------------------------------------------------------
    function allScoresCell = PredictWholeMovie(obj,expi)
      % Predict entire movie: all classifiers, all files
      %
      % allScoresCell: cell array of length nclassifiers
      
      % MERGESTUPDATED
      
      assert(isscalar(expi));

      if isempty(obj.classifier)
        return;
      end
           
      numFlies = obj.GetNumFlies(expi);
      tStartAll = obj.GetTrxFirstFrame(expi);
      tEndAll = obj.GetTrxEndFrame(expi);
      
      [~,firstFly] = max(tEndAll-tStartAll+1);
      perframefile = obj.GetPerframeFiles(expi);
      allperframefns = obj.allperframefns;
      clsNames = obj.classifiernames;
      allScoresCell = cell(obj.nclassifiers,1);

      for iCls = 1:obj.nclassifiers
      
        obj.SetStatus(sprintf('Classifying movie %d:%s, %s...',expi,...
          obj.expnames{expi},clsNames{iCls}));
        
        if obj.isST
          expdir = obj.expdirs{expi};
          classifier = obj.classifier{iCls};
          ppParams = obj.postprocessparams{iCls};
          scoreNorm = obj.windowdata(iCls).scoreNorm;
          usePastOnly = obj.usePastOnly;
          
          allScores = classifySTCore(expdir,classifier,ppParams,...
            scoreNorm,'usePastOnly',usePastOnly,'numTargets',numFlies);
          assert(all(tStartAll==allScores.tStart)); % allScores.tStart currently hardcoded to 1 in classifySTCore
          assert(all(tEndAll==allScores.tEnd));
          allScoresCell{iCls} = allScores;
        else
          windowfeaturescellparams = obj.fastPredict(iCls).windowfeaturescellparams;
          curperframefns = obj.fastPredict(iCls).pffs;
          classifier = obj.fastPredict(iCls).classifier;
          if ~obj.fastPredict(iCls).wfidx_valid
            [~,feature_names] = JLabelData.ComputeWindowDataChunkStatic(curperframefns,...
              allperframefns,perframefile,firstFly,windowfeaturescellparams,1,1);
            obj.fastPredict(iCls) = Predict.fastPredictFindWfidx(...
              obj.fastPredict(iCls),feature_names);
          end
          wfidx = obj.fastPredict(iCls).wfidx;
          
          %       tfile = tempname();
          %       wbar = waitbar(0,'Predicting..');
          %       pause(0.01);
          scoresA = cell(1,numFlies);
          postprocessedscoresA = cell(1,numFlies);
          parfor flies = 1:numFlies
            blockSize = 5000*2;
            tStart = tStartAll(flies);
            tEnd = tEndAll(flies);
            
            scores = nan(1,tEnd);
            
            if tEnd-tStart < 3,
              scores(tStart:tEnd) = -1;
              fprintf('Not predicting for %d fly, Trajectory is too short\n',flies);
            else
              for curt0 = tStart:blockSize:tEnd
                curt1 = min(curt0+blockSize-1,tEnd);
                X = JLabelData.ComputeWindowDataChunkStatic(curperframefns,...
                  allperframefns,perframefile,flies,windowfeaturescellparams,curt0-tStart+1,curt1-tStart+1);
                
                scores(curt0:curt1) = myBoostClassify(X(:,wfidx),classifier);
              end
              %         fid = fopen(tfile,'a');
              %         fprintf(fid,'.');
              %         fclose(fid);
              %         fid = fopen(tfile,'r');
              %         xx = fgetl(fid);
              %         fclose(fid);
              %         numdone = numel(xx);
              %         waitbar(numdone/numFlies,wbar);
              fprintf('Prediction done for %s: %d fly, total number of flies:%d\n',clsNames{iCls},flies,numFlies);
            end
            scoresA{flies} = scores;
          end
          for flies = 1:numFlies
            % AL: not sure why this isn't in the parfor
            postprocessedscoresA{flies} = nan(1,tEndAll(flies));
            postprocessedscoresA{flies}(tStartAll(flies):tEndAll(flies)) = ...
              PostProcessor.PostProcess(scoresA{flies}(tStartAll(flies):tEndAll(flies)),...
              obj.postprocessparams{iCls},obj.windowdata(iCls).scoreNorm);
          end
          
          allScores = struct; % See ScoreFile.allScrs
          allScores.scores = scoresA;
          allScores.tStart = tStartAll;
          allScores.tEnd = tEndAll;
          allScores.postprocessed = postprocessedscoresA;
          allScores.postprocessparams = obj.postprocessparams{iCls};
          for flies = 1:numFlies
            [i0s,i1s] = get_interval_ends(allScores.postprocessed{flies}>0);
            allScores.t0s{flies} = i0s;
            allScores.t1s{flies} = i1s;
          end
          allScores.scoreNorm = obj.windowdata(iCls).scoreNorm;
          
          allScoresCell{iCls} = allScores;
        end
        
      end
      
      obj.ClearStatus();
    end
   
    
    % ---------------------------------------------------------------------
    function PredictWholeMovieNoSave(obj,expi)
      allScores = obj.PredictWholeMovie(expi);
      obj.AddScores(expi,allScores,true);
    end
    
    
    % ---------------------------------------------------------------------
    function allScoresCell = PredictSaveMovie(self,expi,sfn)
    % Predicts for the whole movie and saves the scores.
    % 
    % expi: scalar, index of movie to predict
    % sfn: Optional. Either cellstr of length self.nclassifiers containing 
    %   scorefilenames, or 0, in which case scores are not saved (but are 
    %   returned). If not supplied, current/default scorefilenames for expi 
    %   are used.
    %
    % allScoresCell: cell array of length self.nclassifiers containing
    % scores for each classifier
    
    % MERGESTUPDATED
          
      if nargin < 3
        sfn = self.GetFile('scores',expi);
      end
      assert(isequal(sfn,0) || iscellstr(sfn) && numel(sfn)==self.nclassifiers);
      tfSaveScores = iscellstr(sfn);

      pdExp = self.predictdata{expi};
      nFly = self.nflies_per_exp(expi);
      firstFrms = self.firstframes_per_exp{expi};
      endFrms = self.endframes_per_exp{expi};
      assert(iscell(pdExp) && numel(pdExp)==nFly);
      assert(isequal(nFly,numel(firstFrms),numel(endFrms)));
      
      tf = JLabelData.AllPredictedScoresValid(pdExp,self.nclassifiers);
      if ~all(tf)
        allScoresCell = self.PredictWholeMovie(expi);
      else        
        self.SetStatus(sprintf('Exporting existing scores for movie %d:%s...',expi,self.expnames{expi}));
                
        % compile allScores from prediction data
        allScoresCell = cell(self.nclassifiers,1);
        for iCls = 1:self.nclassifiers
          allScores = ScoreFile.allScrs(nFly);
          allScores = ScoreFile.AllScrsInitFromPredData(allScores,pdExp,iCls,firstFrms,endFrms);
          allScores.postprocessparams = self.postprocessparams{iCls};
          allScores.scoreNorm = self.windowdata(iCls).scoreNorm;
          allScoresCell{iCls} = allScores;
        end
      end
      
      if tfSaveScores
        try
          self.SaveScores(allScoresCell,sfn);
        catch ME
          if nargout > 0
            warning('Could not save scores for experiment %s: %s',self.expnames{expi},getReport(ME));
          else
            error(getReport(ME));
          end
        end
      end
      self.AddScores(expi,allScoresCell,true);
      
      if self.predictdata{expi}{1}(1).loaded_valid(1) 
        % AL: not sure what intent is here; guess condition is whether 
        % scores have been loaded before. But we have already called 
        % .AddScores, so why load on top of the previously loaded scores,
        % which may be different?
        self.LoadScores(expi,sfn);
      end
    end

    
    % ---------------------------------------------------------------------
    function [prediction,T0,T1] = GetPredictedIdx(obj,expi,flies,T0,T1)
      % Get the prediction for experiment expi, target flies, over the time
      % span from T0 to T1.  The returned prediction variable is a scalar
      % structure with two fields: predictedidx and scoresidx, both with
      % nclassifiers rows.
      % - scoresidx holds the raw score output from the classifier.  
      % - predictedidx holds the index of the predicted behavior:
      %    * 1 represents the behavior
      %    * 2 represents no-behavior
      %    * 1.5 represents undecided between beh and no-beh
      %    * 0 indicates no/invalid/out-of-date prediction
      
      if ~isempty(obj.expi) && numel(flies) == numel(obj.flies) && obj.IsCurFly(expi,flies),
        assert(size(obj.predictedidx,1)==obj.ntimelines);
        assert(size(obj.scoresidx,1)==obj.ntimelines);
        if nargin < 4,
          prediction = struct('predictedidx',obj.predictedidx,...
                              'scoresidx', obj.scoresidx);
          T0 = double(obj.t0_curr);
          T1 = double(obj.t1_curr);
        else
          if T0<obj.t0_curr || T1>obj.t1_curr
            T0 = double(max(obj.t0_curr,T0));
            T1 = double(min(obj.t1_curr,T1));
          end
          prediction = struct(...
            'predictedidx', obj.predictedidx(:,T0+obj.labelidx_off:T1+obj.labelidx_off),...
            'scoresidx',  obj.scoresidx(:,T0+obj.labelidx_off:T1+obj.labelidx_off));
        end
        return;
      end
      
      if nargin < 4,
        T0 = max(obj.GetTrxFirstFrame(expi,flies));
        T1 = min(obj.GetTrxEndFrame(expi,flies));
      end
      
      n = T1-T0+1;
      off = 1 - T0;
      prediction = struct('predictedidx',zeros(obj.ntimelines,n),...
                         'scoresidx',zeros(obj.ntimelines,n));      
      for iTL = 1:obj.ntimelines
        idxcurr = obj.predictdata{expi}{flies}(iTL).cur_valid & ...
                  obj.predictdata{expi}{flies}(iTL).t>=T0 & ...
                  obj.predictdata{expi}{flies}(iTL).t<=T1;
      
        prediction.predictedidx(iTL,obj.predictdata{expi}{flies}(iTL).t(idxcurr)+off) = ...
          -sign(obj.predictdata{expi}{flies}(iTL).cur(idxcurr))*0.5+1.5;
        prediction.scoresidx(iTL,obj.predictdata{expi}{flies}(iTL).t(idxcurr)+off) = ...
          obj.predictdata{expi}{flies}(iTL).cur(idxcurr);
      end
    end
    
    
    function predictedidx = PredictedIdxExpandValues(obj,predictedidx)
      % Helper function for GetPredictedIdx.
      % The prediction.predictedidx returned by GetPredictedIdx takes
      % values in the set {0,1,1.5,2}. This function changes
      % representations so that the iTL'th row of predictedidx takes the
      % values:
      % * 0 for no-prediction or undecided (old value=1.5)
      % * obj.iCls2iLbls{iTL}(1) for behavior-predicted
      % * obj.iCls2iLbls{iTL}(2) for nobehavior-predicted.
      
      nTL = obj.ntimelines;
      iCls2iLbls = obj.iCls2iLbl;
      assert(isequal(nTL,size(predictedidx,1),numel(iCls2iLbls)));
      for iTL = 1:nTL
        p = predictedidx(iTL,:);
        predictedidx(iTL,:) = 0;
        iLbls = iCls2iLbls{iTL};
        for i12 = [1 2]
          tf = p==i12;
          predictedidx(iTL,tf) = iLbls(i12);
        end
      end      
    end
    
    
    function predictions = GetPredictionsAllFlies(obj,expi,curt,fliesinrange,iCls)
      % curt: scalar, time desired
      % fliesinrange: vector of flies of interest
      %
      % predictions: row vector, same size as fliesinrange. 1/2/0 for
      % beh/no-beh/no-prediction, resp.
      
      %MERGESTUPDATED
            
      predictions = zeros(size(fliesinrange));
      
      count = 1;
      for curfly = fliesinrange(:)'
        pd = obj.predictdata{expi}{curfly}(iCls);
        idxcurr = pd.cur_valid & pd.t==curt;
        if any(idxcurr)
          predictions(count) = 2 - pd.cur_pp(idxcurr);
        end
        count = count+1;
      end
    end
   
    
    % ---------------------------------------------------------------------
    function SetConfidenceThreshold(obj,thresholds,ndx)
      assert(obj.nclassifiers==1,'Unsupported for multi-classifier projects.');
      assert(isscalar(ndx) && any(ndx==[1 2]));
      obj.confThresholds(1,ndx) = thresholds;
    end
    
    
    % ---------------------------------------------------------------------
    function t = GetConfidenceThreshold(obj,ndx)
      % ndx: indices into obj.labelnames      
      t = obj.confThresholds(ndx); % single index into 2D array 
    end
    
    
  end
  
  methods % (more private)
        
    
    function pd = PredictDataCtor(obj,iExp,iTarget)
      % Construct a new/blank predictdata element.
      
      firstframe = obj.firstframes_per_exp{iExp}(iTarget);
      endframe = obj.endframes_per_exp{iExp}(iTarget);
      nframes = endframe-firstframe+1;
      nanArray = nan(1,nframes);
      falseArray = false(1,nframes);
      
      pd = struct();
      pd.t = (firstframe:endframe);
      pd.cur = nanArray;
      pd.cur_pp = nanArray;
      pd.cur_valid = falseArray;
      pd.old = nanArray;
      pd.old_pp = nanArray;
      pd.old_valid = falseArray;
      pd.loaded = nanArray;
      pd.loaded_pp = nanArray;
      pd.loaded_valid = falseArray;
      pd.timestamp = nanArray;          
    end
      
    % ---------------------------------------------------------------------
    function PredictDataInit(obj,iExp)
      % Dimensions the predictions for experiment iExp.
      
      % MERGESTOK
            
      % Create the cell array for this exeriment
      nTargets = obj.nflies_per_exp(iExp);
      obj.predictdata{iExp} = cell(1,nTargets);
      nTL = obj.ntimelines;
            
      for iTarget = 1:nTargets
        pd = obj.PredictDataCtor(iExp,iTarget);       
        for i = 1:nTL
          obj.predictdata{iExp}{iTarget}(1,i) = pd;
        end
        assert(numel(obj.predictdata{iExp}{iTarget})==nTL);
      end
    end
    
    
    function PredictDataAddClassifier(self)
      % Adds "blank" predictdata to end of each
      % obj.predictdata{iExp}{iTarget}
      
      %nCls = self.nclassifiers;
      nExps = self.nexps;
      for iExp = 1:nExps
        nTargets = self.nflies_per_exp(iExp);
        for iTarget = 1:nTargets
          %assert(numel(self.predictdata{iExp}{iTarget})==nCls);
          self.predictdata{iExp}{iTarget}(1,end+1) = self.PredictDataCtor(iExp,iTarget);
        end
      end
    end

    
    % ---------------------------------------------------------------------
    function PredictDataInvalidate(self)
      % Mark the current and old classifier predictions as invalid.
      nExps=self.nexps;
      for iExp=1:nExps
        nTargets=self.nflies_per_exp(iExp);
        for iTarget=1:nTargets
          nTL = numel(self.predictdata{iExp}{iTarget});
          for iPD = 1:nTL
            self.predictdata{iExp}{iTarget}(iPD).cur(:)=0;
            self.predictdata{iExp}{iTarget}(iPD).cur_pp(:)=0;
            self.predictdata{iExp}{iTarget}(iPD).cur_valid(:)=false;
            self.predictdata{iExp}{iTarget}(iPD).old(:)=0;
            self.predictdata{iExp}{iTarget}(iPD).old_pp(:)=0;          
            self.predictdata{iExp}{iTarget}(iPD).old_valid(:)=false;
          end
        end
      end
      % Do I need to set self.windowdata.scoreNorm to something?
    end  % method
      

    % ---------------------------------------------------------------------
    function PredictDataMoveCurToOld(obj)
      for expi = 1:numel(obj.predictdata)
        for flies = 1:numel(obj.predictdata{expi})
          for iTL = 1:numel(obj.predictdata{expi}{flies})
            obj.predictdata{expi}{flies}(iTL).old = obj.predictdata{expi}{flies}(iTL).cur;
            obj.predictdata{expi}{flies}(iTL).old_valid = obj.predictdata{expi}{flies}(iTL).cur_valid;
            obj.predictdata{expi}{flies}(iTL).old_pp = obj.predictdata{expi}{flies}(iTL).cur_pp;
            obj.predictdata{expi}{flies}(iTL).cur_valid(:) = false;
          end
        end
      end
    end  % method
   
    
    % ---------------------------------------------------------------------
    function PredictLoaded(obj,predictOnlyCurrentFly)
    % PredictLoaded(obj)
    % Runs the classifier on all preloaded window data.
            
    %MERGESTUPDATED

      if isempty(obj.classifier),
        return;
      end
      
      nCls = obj.nclassifiers;

      if nargin < 2,
        predictOnlyCurrentFly = false(1,nCls);
      end
      assert(numel(predictOnlyCurrentFly)==nCls);
      
      % Always predict on labeled data. Mayank June 17 2016
      for iCls = 1:nCls
        for expi = 1:obj.nexps
          for flyNum = 1:obj.nflies_per_exp(expi)
            obj.WindowDataPredictFast(expi,flyNum,iCls,...
              obj.firstframes_per_exp{expi}(flyNum),...
              obj.endframes_per_exp{expi}(flyNum));
          end
        end
      end
      
      for iCls = 1:nCls
        pbs = obj.predictblocks(iCls);
        switch obj.classifiertype{iCls},
          case {'boosting','fancyboosting'}
            for ndx = 1:numel(pbs.t0)
              curex = pbs.expi(ndx);
              flies = pbs.flies(ndx);
              if predictOnlyCurrentFly(iCls) && ~IsCurFly(obj,curex,flies),
                continue;
              end
              
              numcurex = nnz(pbs.expi(:)==curex & pbs.flies(:)==flies);
              numcurexdone = nnz(pbs.expi(1:ndx)==curex & pbs.flies(1:ndx)==flies);
              obj.SetStatus('Predicting for %s: exp %s fly %d ... %d%% done',...
                obj.labelnames{iCls},obj.expnames{curex},flies,round(100*numcurexdone/numcurex));
              obj.PredictFast(curex,flies,pbs.t0(ndx),pbs.t1(ndx),iCls);
            end
        end
      end
      
      obj.NormalizeScores(zeros(nCls,0));
      obj.ApplyPostprocessing();
      obj.ClearStatus();
            
      % transfer to predictidx for current fly
      if ~isempty(obj.expi) && obj.expi > 0 && ~isempty(obj.flies) && all(obj.flies > 0),
        obj.UpdatePredictedIdx();
      end
      
    end    
   
        
    % ---------------------------------------------------------------------
    function FindFastPredictParams(obj)
      obj.fastPredict = Predict.fastPredictInit(obj.fastPredict,...
        obj.classifier,obj.classifierTS,{obj.windowdata.featurenames}); 
    end
    
    
    % ---------------------------------------------------------------------
    function finished = WindowDataPredictFast(obj,expi,flies,clsIdx,t0,t1)
      % Try to predict all classifiers using data cached in windowdata.
      %
      % clsIdx: index of classifier to predict
      % 
      % Effect: update obj.predictdata{expi}{flies}(:) for t0:t1, as 
      % possible. This method will do a partial update.   
      %
      % finished: logical scalar. finished is true if classifier clsIdx
      % is fully predicted over t0:t1. finished will be false if there is a 
      % single timepoint in t0:t1 which is not
      % predicted.
      
      % MERGESTUPDATED 
 
      assert(isscalar(clsIdx) && any(clsIdx==1:obj.nclassifiers));

      finished = false;
      windata = obj.windowdata(clsIdx);
      
      if isempty(windata.exp)
        return;
      end
      idxfly = find(windata.exp == expi & all(bsxfun(@eq,windata.flies,flies),2));
      if isempty(idxfly)
        return;
      end
      
      % ism(i) is whether there is windowdata for t0+i-1
      % idxfly(idx1(i)) is which index of windowdata corresponds to t0+i-1
      [ism,idx1] = ismember(t0:t1,windata.t(idxfly));
      idxism = find(ism);
      if isempty(idxism)
        return;
      end
      % idx(i) is the index into windowdata for t0+idxism(i)-1
      idx = idxfly(idx1(idxism));
      
      % run the classifier: scores(i) is the score for t0+idxism(i)-1
      scores = myBoostClassify(windata.X(idx,:),obj.classifier{clsIdx});
      
      % store in predictdata
      % predictism(i) is whether t0+idxism(i)-1 is in predictdata
      % predictidx(i) is the location in predictdata{expi}{flies} of
      % t0+idxism-1
      [predictism,predictidx] = ismember(double(t0)+idxism-1,obj.predictdata{expi}{flies}(clsIdx).t);
      assert(all(predictism),'There are ts requested that are not in predictdata');
      obj.predictdata{expi}{flies}(clsIdx).cur(predictidx) = scores;
      obj.predictdata{expi}{flies}(clsIdx).timestamp(predictidx) = obj.classifierTS(clsIdx);
      obj.predictdata{expi}{flies}(clsIdx).cur_valid(predictidx) = true;
      
      finished = all(ism);
      % Even if finished==false, predictdata may be updated
    end
      
    
    % ---------------------------------------------------------------------
    function didpredict = PredictFast(obj,expi,flies,t0,t1,iCls)
    % Predict fast by computing only the required window features.
    %
    % iCls: scalar index for classifier to predict.
    %
    % Effect: update obj.predictdata{expi}{flies}(:) for t0:t1
    % Side Effects: update predictblocks, fastPredict
    %
    % didpredict: true if any (new) predictions were generated

    % MERGESTUPDATED
          
      nCls = obj.nclassifiers;
      assert(isscalar(iCls) && any(iCls==1:nCls));
      
      didpredict = false; 
      
      if isempty(obj.classifier) || t0>t1
        return;
      end
      
      % Determine which behaviors/classifiers need predicting;
      pd = obj.predictdata{expi}{flies}(iCls);
      idxcurr_t = pd.t>=t0 & pd.t<=t1;
      if all(pd.cur_valid(idxcurr_t))
        return;
      end
      
      didpredict = true;
      
      % Try WindowDatapredictFast
      finished = obj.WindowDataPredictFast(expi,flies,iCls,t0,t1);
      % ALTODO: optimization, WindowDataPredictFast might have
      % successfully predicted a lot of t0:t1
      if finished
        return;
      end
      
      %%% Full Prediction      
      if any(arrayfun(@(x)isempty(x.classifier),obj.fastPredict)) % isempty(obj.fastPredict(1).classifier)
        obj.FindFastPredictParams();
      end
            
      perframeInMemory = ~isempty(obj.flies) && obj.IsCurFly(expi,flies);
      perframefile = obj.GetPerframeFiles(expi);
      % There is at least one frame in t0:t1 that does not have a current 
      % prediction. We are going to predict on all frames in [t0,t1] even 
      % though some may already have current predictions.
      %
      % ALTODO: Optimization, be smarter about not repredicting everything.
      
      % Initialize missingts, ts which need new prediction
      
      obj.SetStatus('Predicting for classifier %s...',obj.labelnames{iCls});
      
      pd = obj.predictdata{expi}{flies}(iCls);
      missingts = double(pd.t(pd.t>=t0 & pd.t<=t1));
      nmissingts = inf;
      
      DEBUG = false;
      if DEBUG
        fprintf(1,'Exp/fly/cls %d/%d/%d. [t0 t1] = [%.3f %.3f]\n',expi,flies,iCls,t0,t1);
      end
      
      while true
        % - missingts is running state for this loop; each iteration 
        % predicts a subset of missingts and thereby reduces it
        % - Choose a frame missing window data
        % - Try to use an existing predict block if available
        
        t = double(missingts(1));
        curblockndx = obj.predictblocks(iCls).expi==expi & obj.predictblocks(iCls).flies==flies;
        curbs_t0 = obj.predictblocks(iCls).t0(curblockndx);
        curbs_t1 = obj.predictblocks(iCls).t1(curblockndx);
        
        if DEBUG
          fprintf(1,'New iter. t=%d. curbs_t0/curbs_t1: %s/%s\n',t,mat2str(curbs_t0(:)'),mat2str(curbs_t1(:)'));
          disp(obj.predictblocks(iCls));
        end
        
        tfExistingBlock = (t-curbs_t0)>=0 & (t-curbs_t1)<=0 & ...
          (curbs_t1-curbs_t0)>2*obj.predictwindowdatachunk_radius-2;
        if ~isempty(curbs_t0) && any(tfExistingBlock)
          tempndx = find(tfExistingBlock);
          t0 = curbs_t0(tempndx(1));
          t1 = curbs_t1(tempndx(1));
          
          if DEBUG
            fprintf(1,'Using existing block, t0/t1: %.3f/%.3f\n',t0,t1);
          end
        else
          t = median(missingts);
          if ~ismember(t,missingts),
            t = missingts(argmin(abs(t-missingts)));
          end
          
          % bound at start and end frame of these flies
          T0 = max(obj.GetTrxFirstFrame(expi,flies));
          T1 = min(obj.GetTrxEndFrame(expi,flies));
          
          t1 = min(t+obj.predictwindowdatachunk_radius,T1);
          % go backward 2*r to find the start of the chunk
          t0 = max(t1-2*obj.predictwindowdatachunk_radius,T0);
          % go forward 2*r again to find the end of the chunk
          t1 = min(t0+2*obj.predictwindowdatachunk_radius,T1);
          
          
          % Find blocks that overlap with the current interval and merge
          % them into one block.
          overlapping_blocks1 = find(curbs_t0-t0 >= 0 & curbs_t0-t1 <= 0);
          overlapping_blocks2 = find(curbs_t1-t0 >= 0 & curbs_t1-t1 <= 0);
          overlapping_blocks = unique([overlapping_blocks1(:);overlapping_blocks2(:)]);
          if ~isempty(overlapping_blocks),
            t0 = min(t0,min(curbs_t0(overlapping_blocks)));
            t1 = max(t1,max(curbs_t1(overlapping_blocks)));
            todelete = find(curblockndx);
            todelete = todelete(overlapping_blocks);
            obj.predictblocks(iCls).t0(todelete) = [];
            obj.predictblocks(iCls).t1(todelete) = [];
            obj.predictblocks(iCls).flies(todelete) = [];
            obj.predictblocks(iCls).expi(todelete) = [];
            
            % AL 20150529: The expansion of t0/t1 above can lead to
            % additional overlapping blocks, ie we could look for
            % overlapping blocks in a loop. Doesn't look important though.
            if DEBUG
              fprintf(1,'Deleting/merging %d blocks\n',numel(todelete));
            end
          end
          %             overlap_start = find( (t0-curbs_t0)>=0 & (t0-curbs_t1)<=0);
          %             if ~isempty(overlap_start),
          %               t0 = max(curbs_t1(overlap_start))+1;
          %             end
          %
          %             overlap_end = find( (t1-curbs_t0)>=0 & (t1-curbs_t1)<=0);
          %             if ~isempty(overlap_end),
          %               t1 = max(curbs_t0(overlap_end))-1;
          %             end
          
          if t0 <= t1,
            obj.predictblocks(iCls).t0(end+1) = t0;
            obj.predictblocks(iCls).t1(end+1) = t1;
            obj.predictblocks(iCls).expi(end+1) = expi;
            obj.predictblocks(iCls).flies(end+1) = flies;
            
            if DEBUG
              fprintf(1,'Adding block: [%d %d]\n',t0,t1);
            end
          else
            warning('Trying to add interval to predict with t0 = %d > t1 = %d, not doing this. MAYANK, IS THIS RIGHT??',t0,t1);
          end
        end
        
        % - Range of prediction [t0,t1] has been updated
        % - obj.predictblocks(iCls) has been updated as necessary to
        % include a block spanning precisely [t0,t1] (possibly merging blocks etc)
        
        i0 = t0 - obj.GetTrxFirstFrame(expi,flies) + 1;
        i1 = t1 - obj.GetTrxFirstFrame(expi,flies) + 1;
        
        X = [];
        
        %%% Compute window features
        perframedata_cur = obj.perframedata;
        allperframefns = obj.allperframefns;
        windowfeaturescellparams = obj.fastPredict(iCls).windowfeaturescellparams;
        pffs = obj.fastPredict(iCls).pffs;
        feature_names_list = cell(1,numel(pffs));
        x_curr_all = cell(1,numel(pffs));
        
        [~,pfidx] = ismember(pffs,obj.allperframefns);
        tmp_perframedata_cur = perframedata_cur(pfidx);
        stInfo = obj.stInfo;
        try
          parfor j = 1:numel(pffs),
            
            fn = pffs{j};
            
            ndx = find(strcmp(fn,allperframefns)); %#ok<PROPLC>
            if perframeInMemory,
%               perframedata = perframedata_cur{ndx};  %#ok
              perframedata = tmp_perframedata_cur{j}; %#ok<PROPLC>
            else
              perframedata = readPFData(perframefile{ndx},flies(1),stInfo);  %#ok
              perframedata = perframedata{1};  %#ok
            end
            
            i11 = min(i1,numel(perframedata)); %#ok<PROPLC>
            [x_curr,cur_f] = ...
              ComputeWindowFeatures(perframedata,...
              windowfeaturescellparams.(fn){:},'t0',i0,'t1',i11);  %#ok
            
            if i11 < i1,
              x_curr(:,end+1:end+i1-i11) = nan;
            end
            
            x_curr_all{j} = single(x_curr);
            feature_names_list{j} = cur_f;
          end  % parfor
          
        catch ME,
          
          uiwait(warndlg(sprintf('Could not predict:%s',ME.message)));
          
        end
        
        % Passing on to next stage: x_curr_all, feature_names_list
        
        % Form feature matrix X
        for j = 1:numel(pffs),
          fn = pffs{j};
          x_curr = x_curr_all{j}; % will be empty if err occurred
          % add the window data for this per-frame feature to X
          nold = size(X,1);
          nnew = size(x_curr,2);
          if nold > nnew,
            warning(['Number of examples for per-frame feature %s does not '...
              'match number of examples for previous features'],fn);
            x_curr(:,end+1:end+nold-nnew) = nan;
          elseif nnew > nold && ~isempty(X),
            warning(['Number of examples for per-frame feature %s does not '...
              'match number of examples for previous features'],fn);
            X(end+1:end+nnew-nold,:) = nan;
          end
          X = [X,x_curr']; %#ok<AGROW>
        end
        
        if ~obj.fastPredict(iCls).wfidx_valid,
          feature_names = {};
          for ndx = 1:numel(feature_names_list)
            fn = pffs{ndx};
            feature_names = [feature_names,cellfun(@(s) [{fn},s],...
              feature_names_list{ndx},'UniformOutput',false)]; %#ok<AGROW>
          end
          obj.fastPredict(iCls) = Predict.fastPredictFindWfidx(obj.fastPredict(iCls),feature_names);
          assert(obj.fastPredict(iCls).wfidx_valid);
        end
        
        scores = myBoostClassify(X(:,obj.fastPredict(iCls).wfidx),obj.fastPredict(iCls).classifier);
        
        curndx = obj.predictdata{expi}{flies}(iCls).t>=t0 & ...
          obj.predictdata{expi}{flies}(iCls).t<=t1;
        obj.predictdata{expi}{flies}(iCls).cur(curndx) = scores;
        obj.predictdata{expi}{flies}(iCls).timestamp(curndx) = obj.classifierTS(iCls);
        obj.predictdata{expi}{flies}(iCls).cur_valid(curndx) = true;
        
        missingts(missingts >= t0 & missingts <= t1) = [];
        
        if isempty(missingts),
          break;
        end
        
        nmissingtsnew = numel(missingts);
        if nmissingtsnew >= nmissingts,
          errordlg('Sanity check: Number of frames missing window features did not decrease. Breaking out of loop.');
          break;
        end
        nmissingts = nmissingtsnew;
      end
      
      
      obj.ClearStatus();
    end    
    
    
     % ---------------------------------------------------------------------
    function UpdatePredictedIdx(obj)
      % Updates obj.predictedidx, obj.scoresidx, obj.scoreTS 
      % to match what's in obj.predictdata, obj.expi, obj.flies, 
      % obj.t0_curr, and obj.t1_curr.  For instance, this might be used to 
      % update those variables after obj.expi changes.
      
      %MERGESTREVIEWED
      
      if obj.expi == 0,
        return;
      end
      
      n = obj.t1_curr - obj.t0_curr + 1;
      nTL = obj.labelidx.nTL;      
     
      obj.predictedidx = zeros(nTL,n);
      obj.scoresidx = zeros(nTL,n);
      obj.scoresidx_old = zeros(nTL,n);
      obj.scoreTS = zeros(nTL,n);
      
      if ~isempty(obj.predictdata) && ~isempty(obj.predictdata{obj.expi}) 
        % Overwrite by scores from predictdata.
        
        pdArr = obj.predictdata{obj.expi}{obj.flies};
        assert(numel(pdArr)==nTL);
        for iTL = 1:nTL 
          pd = pdArr(iTL);
          
          idxcurr = pd.cur_valid;
          idx = pd.t(idxcurr)-obj.t0_curr+1;
          obj.predictedidx(iTL,idx) = -sign(pd.cur(idxcurr))*0.5 + 1.5;
          obj.scoresidx(iTL,idx) = pd.cur(idxcurr);
          %obj.scoreTS(iTL,idx) = obj.classifierTS(iTL); 
          % ALTODO: .scoreTS does not appear to be used by anything.
          % This assignment causing problems due to order-of-initialization
          % (classifierTS is empty in certain call chains), just comment
          % out for now.
        end
      end
      %obj.UpdateErrorIdx();
    end

  end
  
  methods (Static,Access=private)
    
    function tf = AllPredictedScoresValid(predDataExp,nclassifier)
      % predDataExp: JLD.predictdata for one experiment, eg
      % obj.predictdata{expi}
      % tf: nclassifier x 1 logical. tf(i) is true if entire timeline for
      % all flies is predicted/valid for classifier i
      
      %MERGESTUPDATED
      
      nfly = numel(predDataExp);
      tf = true(nclassifier,1);
      for iCls = 1:nclassifier
        for fly = 1:nfly
          pd = predDataExp{fly};
          assert(numel(pd)==nclassifier);
          if ~all(pd(iCls).cur_valid)
            tf(iCls) = false;
            break;
          end
        end
      end
    end
    
  end  
  
  %% Scores/Postprocessing/Stats
  
  methods
    
    
    % ---------------------------------------------------------------------
    function AddScores(obj,expi,allScoresCell,updateCurrent) 
      % Set .predictdata from allScoresCell
      %
      % allScoresCell: cell array of length nclassifiers
      % updateCurrent: logical scalar. If true, then .predictdata{}{}().cur
      % and .cur_valid are updated; otherwise, .loaded and .loaded_valid
      % are updated.
      % 
      % This sure seems like it should be a private method, but it's called
      % by JLabelGUIData.  -- ALT, Apr 18, 2013
      
      % MERGESTUPDATED
      
      if isscalar(allScoresCell) && isstruct(allScoresCell)
        allScoresCell = {allScoresCell};
      end
      assert(iscell(allScoresCell) && numel(allScoresCell)==obj.ntimelines);
      
      if updateCurrent
        fldScore = 'cur';
        fldValid = 'cur_valid';
      else
        fldScore = 'loaded';
        fldValid = 'loaded_valid';
      end
        
      obj.SetStatus('Updating Predictions ...');
        
      for ibeh = 1:obj.ntimelines
        nFly = numel(allScoresCell{ibeh}.scores);
        assert(isequal(nFly,obj.GetNumFlies(expi),numel(obj.predictdata{expi})));
        for ndx = 1:nFly
          tStart = allScoresCell{ibeh}.tStart(ndx);
          tEnd = allScoresCell{ibeh}.tEnd(ndx);
          curScores = allScoresCell{ibeh}.scores{ndx}(tStart:tEnd);
          Nscore = tEnd-tStart+1;
          Npredict = numel(obj.predictdata{expi}{ndx}(ibeh).(fldValid));
          if Npredict < Nscore
            warningNoTrace('JLabelData:scoreSizeMismatch',...
              'Scores for experiment %d, behavior %d, fly %d have more elements (%d) than expected (%d). Truncating scores.',...
              expi,ibeh,ndx,Nscore,Npredict);
            curScores = curScores(1:Npredict);            
          elseif Npredict > Nscore
            warningNoTrace('JLabelData:scoreSizeMismatch',...
              'Scores for experiment %d, behavior %d, fly %d have fewer elements (%d) than expected (%d). Padding with nans.',...
              expi,ibeh,ndx,Nscore,Npredict);
            curScores(end+1:Npredict) = nan;
          end
          obj.predictdata{expi}{ndx}(ibeh).(fldScore)(:) = curScores;
          obj.predictdata{expi}{ndx}(ibeh).(fldValid)(:) = true;
        end
      end

      if ~isempty(obj.postprocessparams)
        [success,msg] = obj.ApplyPostprocessing();
        if ~success
          uiwait(warndlg(['Couldn''t apply postprocessing to the scores: ' msg]));
        end
      elseif ~updateCurrent,
        assert(false,'Deprecated codepath');
        for ibeh = 1:obj.ntimelines 
          for ndx = 1:numel(allScoresCell{ibeh}.loaded) % Out of date fieldname
            tStart = allScoresCell{ibeh}.tStart(ndx);
            tEnd = allScoresCell{ibeh}.tEnd(ndx);
            if isfield(allScoresCell{ibeh},'postprocessedscores'); % Out of date fieldname
              obj.postprocessparams = allScoresCell{ibeh}.postprocessparams;
              curpostprocessedscores = allScoresCell{ibeh}.postprocessedscores{ndx}(tStart:tEnd);
              obj.predictdata{expi}{ndx}(ibeh).loaded_pp(:) = curpostprocessedscores;
          else
              obj.predictdata{expi}{ndx}(ibeh).loaded_pp(:) = 0;
            end
          end
        end
      end
      
      obj.UpdatePredictedIdx();
      obj.ClearStatus();
    end

    
    % ---------------------------------------------------------------------
    function LoadScores(obj,expi,scorefns)
      
      %MERGESTUPDATED 
      
      assert(iscellstr(scorefns));
      Nscore = numel(scorefns);
      assert(Nscore==obj.nclassifiers); % must be one scorefile per classifier
      
      for i = 1:Nscore
        sfn = scorefns{i};
        if ~exist(sfn,'file')
          warndlg('Score file %s does not exist. Not loading scores',sfn);
          return;
        end
      end
      
      scorefnstr = civilizedStringFromCellArrayOfStrings(scorefns);
      obj.SetStatus('Loading scores for experiment %s from %s',obj.expnames{expi},scorefnstr);
      
      behNames = Labels.verifyBehaviorNames(obj.labelnames);
      assert(numel(behNames)==Nscore);
      allScoresCell = cell(1,Nscore);

%       winDataScoreNorm = obj.windowdata.scoreNorm;
%       % Windowdata.scoreNorm will be updated from scorefiles as appropriate
%       if numel(winDataScoreNorm)>Nscore
%         warningNoTrace('JLabelData:windataScoreNormTooBig','windowdata.scoreNorm has too many elements. Truncating, then loading from scorefile(s).');
%         obj.windowdata.scoreNorm = winDataScoreNorm(1:Nscore);
%       elseif numel(winDataScoreNorm<Nscore)
%         %warningNoTrace('JLabelData:windataScoreNormTooBig','windowdata.scoreNorm has too few elements. Padding with nans, then loading from scorefile(s).');
%         obj.windowdata.scoreNorm(end+1:Nscore) = nan;
%       end        
      for i = 1:Nscore
        sfn = scorefns{i};
        S = load(sfn); % Could use ScoreFile.load here
        
        % check the behavior name
        if isfield(S,'behaviorName')
          if ~strcmp(S.behaviorName,behNames{i})
            warningNoTrace('LoadScores:possibleBehaviorMismatch',...
              'Possible behavior mismatch. Behavior name in score file: %s. Expected: %s.',S.behaviorName,behNames{i});
          end
        end
        
        if ~isempty(obj.classifierTS) && obj.classifierTS(i)>0
          if S.timestamp~=obj.classifierTS(i)
            uiwait(warndlg(sprintf(['Scores were computed using a classifier trained on %s'...
              ' while the current classifier was trained on %s'],datestr(S.timestamp),...
              datestr(obj.classifierTS(i))))),
          end
        end
        % AL 20141211: classifierfilename appears to be deprecated field
        if isfield(S,'classifierfilename') %~isempty(whos('-file',sfn,'classifierfilename'))
          classifierfilename = S.classifierfilename;
        else
          classifierfilename = '';
        end
        
        % Set windowdata.scoreNorm if necessary
        %if isempty(winDataScoreNorm) || all(isnan(winDataScoreNorm)) || all(winDataScoreNorm==0),
        if isnan(obj.windowdata(i).scoreNorm) || obj.windowdata(i).scoreNorm==0
          if isfield(S.allScores,'scoreNorm')
            scoreNorm = S.allScores.scoreNorm;
          elseif isa(classifierfilename,'Macguffin')
            % AL 20141211: should be deprecated codepath
            scoreNorm = classifierfilename.classifierStuff.windowdata.scoreNorm;
          elseif exist(classifierfilename,'file') && ~isempty(whos('-file',classifierfilename,'scoreNorm'))
            scoreNorm = load(classifierfilename,'scoreNorm');
          else
            scoreNorm = 1;
            uiwait(warndlg(sprintf('Score file %s did not have the score normalization. Setting it to 1',sfn)));
          end
          
          assert(isscalar(scoreNorm),'Expected scalar scoreNorm.');
          obj.windowdata(i).scoreNorm = scoreNorm;
        end
        
        allScoresCell{i} = S.allScores;
      end
      
      obj.AddScores(expi,allScoresCell,false);
      
      obj.ClearStatus();
    end
    
    
    % ---------------------------------------------------------------------
    function LoadScoresDefault(obj,expi)
      % MERGESTUPDATED
      
      [scorefns,tffound] = obj.GetFile('scores',expi);
      assert(iscellstr(scorefns));
      if any(~tffound)
        warndlg(sprintf('Missing scorefile(s) for experiment %d:%s',expi,obj.expdirs{expi}));
      else
        obj.LoadScores(expi,scorefns);
      end
    end
    
    
    % ---------------------------------------------------------------------
    function tf = HasLoadedScores(obj,iCls)
      % tf: true if there is at least one exp, and there are loaded scores for all exps/flies
      
      %MERGESTUPDATED
      
      assert(isscalar(iCls) && any(iCls==1:obj.nclassifiers));
      
      if obj.nexps==0
        tf = false;
        return;
      end
        
      tf = true;
      for expi = 1:obj.nexps
        for flies = 1:obj.nflies_per_exp(expi)
          pdArr = obj.predictdata{expi}{flies};
          if ~pdArr(iCls).loaded_valid(1)
            tf = false;
            return;
          end
        end
      end      
    end    

    
    % ---------------------------------------------------------------------
    function [success,msg] = ApplyPostprocessing(obj,expis,allflies)
    % Applies postprocessing to current, loaded scores.
    %
    % Effect: update obj.predictdata{expis}{...}(:).cur_pp and
    % obj.predictdata{expis}{...}(:).loaded_pp
    
    %MERGESTUPDATED
            
      if nargin < 2,
        expis = 1:obj.nexps;
      end
      if nargin >= 3,
        if ~iscell(allflies),
          allflies = {allflies};
        end
      else
        allflies = {};
      end
      for i = 1:numel(expis),
        if numel(allflies) < i,
          endx = expis(i);
          allflies{i} = 1:obj.nflies_per_exp(endx);
        end
      end
      
      assert(numel(expis)==numel(allflies) && iscell(allflies));
      % allflies{i} contains list of flies for exp expis{i}
          
      %fprintf('Calling ApplyPostprocessing for %d experiments and %d flies...\n',numel(expis),sum(cellfun(@numel,allflies)));  
      
      nCls = obj.nclassifiers;
      assert(numel(obj.windowdata)==nCls);
      for ibeh = 1:nCls
        scoreNorm = obj.windowdata(ibeh).scoreNorm;
        ppparams = obj.postprocessparams{ibeh};
        for expii = 1:numel(expis),
          endx = expis(expii);
          for flies = allflies{expii},
            idx = find(obj.predictdata{endx}{flies}(ibeh).cur_valid);
            ts = obj.predictdata{endx}{flies}(ibeh).t(idx);
            [sortedts, idxorder] = sort(ts);
            gaps = find((sortedts(2:end) - sortedts(1:end-1))>1)+1; 
            gaps = [1;gaps';numel(ts)+1];
            for ndx = 1:numel(gaps)-1
              % Loop over "nongap" or "consecutive-t" subseqs
              % - idxorder(gaps(ndx):gaps(ndx+1)-1) appears to be indices into ts
              % for a consecutive subseq
              % - so curidx represent indices into fields of predictdata{}{}(ibeh)
              curidx = idx(idxorder(gaps(ndx):gaps(ndx+1)-1));
              curs = obj.predictdata{endx}{flies}(ibeh).cur(curidx);
              obj.predictdata{endx}{flies}(ibeh).cur_pp(curidx) = ...
                PostProcessor.PostProcess(curs,ppparams,scoreNorm);
            end
          end
        end
      end
      
      for ibeh = 1:nCls
        scoreNorm = obj.windowdata(ibeh).scoreNorm;
        ppparams = obj.postprocessparams{ibeh};
        for expii = 1:numel(expis),
          endx = expis(expii);
          for flies = allflies{expii},
            curidx = obj.predictdata{endx}{flies}(ibeh).loaded_valid;
            curt = obj.predictdata{endx}{flies}(ibeh).t(curidx);
            if any(curt(2:end)-curt(1:end-1) ~= 1)
              msg = 'Scores are not in order';
              success = false;
              return;
              % ALTODO early return has left obj.predictdata partially
              % updated. Is this an assert?
            end
            curs = obj.predictdata{endx}{flies}(ibeh).loaded(curidx);
            obj.predictdata{endx}{flies}(ibeh).loaded_pp(curidx) = ...
              PostProcessor.PostProcess(curs,ppparams,scoreNorm);
          end %flies
        end
      end
      
      msg = ''; 
      success = true;     
    end  % method

    
    % ---------------------------------------------------------------------
    function SaveScores(self,allScoresCell,sfn)
      % Save prediction scores for a whole experiment.
      % 
      % allScoresCell: cell array, nclassifiers elements, each el is an
      %  allScore (see ScoreFile)
      % sfn: cellstr, full filenames to be saved
      
      % MERGESTUPDATED
      
      nCls = self.nclassifiers;
      assert(iscell(allScoresCell) && numel(allScoresCell)==nCls);
      assert(iscellstr(sfn) && numel(sfn)==nCls);

      if ~self.userHasSpecifiedEverythingFileName
        error('JLabelData:noJabFileNameSpecified', ...
          'A .jab file name must be specified before scores can be saved.');
      end
      
      sf = ScoreFile;
      sf.jabFileNameAbs = self.everythingFileNameAbs;
      sf.version = self.version;
      timestamp = self.classifierTS;
      behnames = self.classifiernames;
      assert(isequal(nCls,numel(timestamp),numel(behnames)));
      
      for iCls = 1:nCls
        sf.allScores = allScoresCell{iCls};
        sf.behaviorName = behnames{iCls};
        sf.timestamp = timestamp(iCls);
        save(sf,sfn{iCls});
      end
    end
    
    
    % ---------------------------------------------------------------------
    function scores = GetValidatedScores(obj,expi,flies,T0,T1,clsIdx)
      % clsIdx: array of classifier indices. Defaults to 1:obj.nclassifiers.
      % scores: numel(clsIdx)-by-(T1-T0+1) array
      
      %MERGESTUPDATED
      
      if nargin<4 || isempty(T0)
        T0 = max(obj.GetTrxFirstFrame(expi,flies));
      end
      if nargin<5 || isempty(T1)
        T1 = min(obj.GetTrxEndFrame(expi,flies));
      end
      if nargin<6
        clsIdx = 1:obj.nclassifiers;
      end
      
      n = T1-T0+1;
      off = 1-T0;
      nCls = numel(clsIdx);
      scores = zeros(nCls,n); % AL: nan instead?
      
      for i = 1:nCls
        iCls = clsIdx(i);
        wd = obj.windowdata(iCls);
        if ~isempty(wd.scores_validated)
          idxcurr = obj.FlyNdx(expi,flies,iCls) & wd.t>=T0 & wd.t<=T1;
          scores(i,wd.t(idxcurr)+off) = wd.scores_validated(idxcurr);
        end
      end      
    end
 
  
    % ---------------------------------------------------------------------
    function [scores,predictions] = GetLoadedScores(obj,expi,flies,T0,T1,clsIdx)
      %MERGESTUPDATED
      if nargin<4 || isempty(T0)
        T0 = max(obj.GetTrxFirstFrame(expi,flies));
      end
      if nargin<5 || isempty(T1)
        T1 = min(obj.GetTrxEndFrame(expi,flies));
      end
      if nargin<6
        clsIdx = 1:obj.nclassifiers;
      end
      
      [scores,predictions] = obj.GetScoresCore(expi,flies,...
        'loaded','loaded_valid','loaded_pp',T0,T1,clsIdx);
    end
    
    
    % ---------------------------------------------------------------------
    function [scores,predictions] = GetPostprocessedScores(obj,expi,flies,T0,T1,clsIdx)
      %MERGESTUPDATED
      if nargin<4 || isempty(T0)
        T0 = max(obj.GetTrxFirstFrame(expi,flies));
      end
      if nargin<5 || isempty(T1)
        T1 = min(obj.GetTrxEndFrame(expi,flies));
      end
      if nargin<6
        clsIdx = 1:obj.nclassifiers;
      end
      
      if obj.HasCurrentScores
        [scores,predictions] = obj.GetScoresCore(expi,flies,...
          'cur','cur_valid','cur_pp',T0,T1,clsIdx);
      else
        [scores,predictions] = obj.GetLoadedScores(expi,flies,T0,T1,clsIdx);
      end      
    end
        
    
    % ---------------------------------------------------------------------
    function scores = GetOldScores(obj,expi,flies,clsIdx)
      if nargin<4
        clsIdx = 1:obj.nclassifiers;
      end
      
      %MERGESTUPDATED
      
      T0 = max(obj.GetTrxFirstFrame(expi,flies));
      T1 = min(obj.GetTrxEndFrame(expi,flies));
      scores = obj.GetScoresCore(expi,flies,'old','old_valid','old_valid',T0,T1,clsIdx);
      % Using arbitrary/random field old_valid for fldPP; this arg only
      % used for computing predictions, which we are not using
    end

    
    % ---------------------------------------------------------------------
    function [stats,flyStats] = GetFlyStats(obj,expi,flyNum)
      % Calculates statistics such as number of labeled bouts, predicted bouts
      % and change in scores.
      %
      % stats: scalar struct, overall stats for exp/fly 
      % flyStats: nclassifier-by-1 struct array classifier-specific stats
            
      % MERGESTUPDATED
      
      assert(isscalar(flyNum));
      
      obj.SetStatus('Computing stats for %s, target %d',obj.expnames{expi},flyNum);
      
      stats = struct();
      flyStats = cell2struct(cell(0,obj.nclassifiers),{});
      
      % General stats
      stats.endframe = obj.endframes_per_exp{expi}(flyNum);
      stats.firstframe = obj.firstframes_per_exp{expi}(flyNum);
      stats.trajLength = stats.endframe-stats.firstframe+1;      
      if obj.hassex
        if obj.hasperframesex
          sexfrac = obj.GetSexFrac(expi,flyNum);
          stats.sexfrac = round(100*sexfrac.M);
        else
          stats.sexfrac = 100*strcmpi(obj.GetSex(expi,flyNum),'M');
        end
      else
        stats.sexfrac = [];
      end
      
      % Compile frames/labels for this fly
      obj.StoreLabelsAndPreLoadWindowData();

      lblIdx2ClsIdx = obj.iLbl2iCls;
      clsIdx2LblIdx = obj.iCls2iLbl;      
      curbouts = zeros(1,0); % bout counter, one element per bout, values are CLASSIFIER indices (in 1:nclassifier)
      curts = zeros(1,0); % all labeled frame indices for this fly (may contain repeats)
      curlabels = zeros(1,0); % label vector for curts; values are LABEL indices (in 1:2*nclassifier)
      lblsExp = obj.labels(expi);
      lblsExpFlies = lblsExp.flies;
      if isequal(lblsExpFlies,[])
        lblsExpFlies = lblsExpFlies(:);
      end
      [tf,iFly] = ismember(flyNum,lblsExpFlies,'rows');
      if tf
        nBouts = numel(lblsExp.t0s{iFly});
        for iBout = 1:nBouts
          t0 = lblsExp.t0s{iFly}(iBout);
          t1 = lblsExp.t1s{iFly}(iBout);
          name = lblsExp.names{iFly}{iBout};
          lblIdx = find(strcmp(name,obj.labelnames));
          assert(isscalar(lblIdx));
          numFrames = t1-t0;
          
          curts(1,end+1:end+numFrames) = t0:(t1-1);
          curlabels(1,end+1:end+numFrames) = lblIdx;
          curbouts(1,end+1) = lblIdx2ClsIdx(lblIdx); %#ok<AGROW>
        end
      end
      
      % General label stats
      nCls = obj.nclassifiers;
      stats.nbouts = numel(curbouts);
      stats.posframes = nnz(curlabels<=nCls);
      stats.negframes = nnz(curlabels>nCls); 
      stats.totalframes = numel(curts); % with multiple classifiers, this could exceed number of frames in track 
      
      for iCls = 1:nCls
        % label stats
        tmp = clsIdx2LblIdx{iCls};
        posLblIdx = tmp(1);
        negLblIdx = tmp(2);

        tmp = struct();
        tmp.nbouts = nnz(curbouts==iCls);        
        tmp.posframes = nnz(curlabels==posLblIdx);
        tmp.negframes = nnz(curlabels==negLblIdx);
        tmp.totalframes = tmp.posframes + tmp.negframes;        
        prefix = fif(obj.gtMode,'gt_','');
        flds = fieldnames(tmp);
        for f = flds(:)', f=f{1}; %#ok<FXSET>
          flyStats(iCls).([prefix f]) = tmp.(f);
        end
        
        pd = obj.predictdata{expi}{flyNum}(iCls);
        
        if ~isempty(pd.loaded_valid) && pd.loaded_valid(1)
          idxcurr = pd.loaded_valid;
          flyStats(iCls).nscoreframes_loaded = nnz(idxcurr);
          flyStats(iCls).nscorepos_loaded = nnz(pd.loaded(idxcurr)>0);
          flyStats(iCls).nscoreneg_loaded = nnz(pd.loaded(idxcurr)<0);
        else
          flyStats(iCls).nscoreframes_loaded = [];
          flyStats(iCls).nscorepos_loaded = [];
          flyStats(iCls).nscoreneg_loaded = [];
        end
        
        tfCls = curlabels==posLblIdx | curlabels==negLblIdx;
        curtsCls = curts(tfCls);
        curlabelsCls = curlabels(tfCls);
        assert(numel(unique(curtsCls))==numel(curtsCls),...
          'For a given classifier each frame may be labeled at most once.');
        curNdx = pd.cur_valid;
        if any(curNdx) && ~isempty(obj.classifier)
          
          % Ignore labels that don't have predicted scores.
          haveScores = curNdx(curtsCls - obj.GetFirstFrames(expi,flyNum)+1); 
          curtsCls(~haveScores) = [];
          curlabelsCls(~haveScores) = [];
          
          if ~isempty(curlabelsCls)
            curWScores = pd.cur(curtsCls - obj.GetFirstFrames(expi,flyNum)+1);
            curPosMistakes = nnz( curWScores<0 & curlabelsCls==posLblIdx );
            curNegMistakes = nnz( curWScores>0 & curlabelsCls==negLblIdx );
          else
            curPosMistakes = [];
            curNegMistakes = [];
          end
          
          curScores = pd.cur(curNdx);          
          flyStats(iCls).nscoreframes = nnz(curNdx);
          flyStats(iCls).nscorepos = nnz(curScores>0);
          flyStats(iCls).nscoreneg = nnz(curScores<0);
          flyStats(iCls).errorsPos = curPosMistakes;
          flyStats(iCls).errorsNeg = curNegMistakes;
        else
          flyStats(iCls).nscoreframes = [];
          flyStats(iCls).nscorepos = [];
          flyStats(iCls).nscoreneg = [];
          flyStats(iCls).errorsPos = [];
          flyStats(iCls).errorsNeg = [];
        end
        
        flyStats(iCls).one2two = [];
        flyStats(iCls).two2one = [];
        if ~isempty(obj.classifier_old{iCls})
          curNdx = pd.old_valid & pd.cur_valid;
          if nnz(curNdx)
            flyStats(iCls).one2two = nnz(pd.cur(curNdx)<0 & pd.old(curNdx)>0);
            flyStats(iCls).two2one = nnz(pd.cur(curNdx)>0 & pd.old(curNdx)<0);
          end
        end
        
        flyStats(iCls).validatedErrorsPos = [];
        flyStats(iCls).validatedErrorsNeg = [];
        if ~isempty(obj.windowdata(iCls).scores_validated)
          curNdx = obj.FlyNdx(expi,flyNum,iCls);
          if nnz(curNdx)
            curScores = obj.windowdata(iCls).scores_validated(curNdx);
            curLabels = obj.windowdata(iCls).labelidx_new(curNdx);
            assert(all(curLabels==posLblIdx | curLabels==negLblIdx));
            
            curPosMistakes = nnz( curScores(:)<0 & curLabels(:)==posLblIdx );
            curNegMistakes = nnz( curScores(:)>0 & curLabels(:)==negLblIdx );
            
            flyStats(iCls).validatedErrorsPos = curPosMistakes;
            flyStats(iCls).validatedErrorsNeg = curNegMistakes;
          end
        end
        
        % MAYANK_JAN15_2016: Do we need gtsuggestions if its normal mode?
        if obj.IsGTMode,
          flyStats(iCls).gt_suggestion_frames = nnz(obj.GetGTSuggestionIdx(expi,flyNum));
        end
        
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
      
      obj.ClearStatus();
    end

   
    % ---------------------------------------------------------------------
    function scores = NormalizeScores(obj,scores,clsIdx)
      % scores: numel(clsIdx)-by-nsamp score array
      % clsIdx: vector of classifier indices labeling rows of scores. 
      % Defaults to 1:obj.nclassifiers.
      % 
      % Normalize the given scores, using the scoreNorm value in self.
      % Seems like it might make sense to add an option to all the
      % methods that get scores out of the JLabelData option, so that
      % callers can request normalized scores, and then make this a private
      % function.  --ALT, Apr 19, 2013   
      %
      % Effect: scores normalized. Each row scores(i,:) normalized by
      % obj.windowdata(clsIdx(i)).scoreNorm. After normalization, scores
      % will be in range [-1,1].
      %
      % Side effect: obj.windowdata(:).scoreNorm initialized if necessary,
      % using obj.windowdata(:).X and current classifiers
      
      %MERGESTUPDATED
      
      if ~exist('clsIdx','var')
        clsIdx = 1:obj.nclassifiers;
      end
      nCls = numel(clsIdx);
      assert(size(scores,1)==nCls);

      for i = 1:nCls
        iCls = clsIdx(i);
        cls = obj.classifier{iCls};
        
        if isempty(obj.windowdata(iCls).scoreNorm) || isnan(obj.windowdata(iCls).scoreNorm)
          % Need to initialize scoreNorm
          
          if isempty(obj.windowdata(iCls).X) || isempty(cls)
            % Can't compute scoreNorm

            % Old comment:
            % Just return the unaltered scores in this case
            % Note that these unaltered scores can have elements outside of
            % [-1,+1].
            % Is this really what we want in this case?
            
            if obj.isST && prctile(abs(scores(iCls,:)),80)>0,
              % not the best way but the only I can think of.
              % MK 20170605
              obj.windowdata(iCls).scoreNorm = prctile(abs(scores(iCls,:)),80);
            else
              obj.windowdata(iCls).scoreNorm = nan;
            end
          else
            wScores = myBoostClassify(obj.windowdata(iCls).X,cls);
            obj.windowdata(iCls).scoreNorm = prctile(abs(wScores),80);
          end
        end

        scoreNorm = obj.windowdata(iCls).scoreNorm;
        tfSm = scores(i,:)<-scoreNorm;
        tfLg = scores(i,:)>scoreNorm;
        scores(i,tfSm) = -scoreNorm;
        scores(i,tfLg) = scoreNorm;
        scores(i,:) = scores(i,:)/scoreNorm;
        % isnan(scoreNorm) => scores(i,:) is nan
      end
    end
    
    
    % ---------------------------------------------------------------------
    function has = HasCurrentScores(obj)
      % has: scalar logical. If true, at least one
      % experiment/fly/classifier has a valid predicted score (for at least
      % one frame)
      
      %MERGEST UPDATED
      
      has = false;
      for expi = 1:obj.nexps
        for flies = 1:obj.nflies_per_exp(expi)
          for ibeh = 1:obj.ntimelines
            if any(obj.predictdata{expi}{flies}(ibeh).cur_valid)
              has = true;
              return;
            end
          end
        end
      end
    end
    
    
    % ---------------------------------------------------------------------
    function params = GetPostprocessingParams(obj)
      params = obj.postprocessparams;
    end
    
    
    % ---------------------------------------------------------------------
    function [success,msg] = SetPostprocessingParams(obj,params)
      assert(iscell(params) && numel(params)==obj.nclassifiers);
      obj.postprocessparams = params;
      [success,msg] = obj.ApplyPostprocessing();
    end
    
    
    % ---------------------------------------------------------------------
    function blen = GetPostprocessedBoutLengths(obj,iCls)
      % blen: row vector, bout length for classifier iCls over all
      % exps/flies
      
      %MERGESTUPDATED
      
      blen = zeros(1,0);
           
      if obj.HasCurrentScores()
        % For predicted scores.
        for endx = 1:obj.nexps
          for flies = 1:obj.nflies_per_exp(endx)
            pd = obj.predictdata{endx}{flies}(iCls);
            idx = find(pd.cur_valid); % AL20141215 
            ts = pd.t(idx);            
            [sortedts,idxorder] = sort(ts);
            % AL 20141215: added +1 to next line, cf ApplyPostProcessing
            gaps = find((sortedts(2:end) - sortedts(1:end-1))>1)+1; 
            gaps = [1;gaps';numel(ts)+1];
            for ndx = 1:numel(gaps)-1
              % loop over 'contiguous' segments of time
              curidx = idx(idxorder(gaps(ndx):gaps(ndx+1)-1)); % indices into pd.t, pd.cur_valid, pd.cur_pp for current time segment
              assert(isequal(size(pd.t),size(pd.cur_valid),size(pd.cur_pp)));
              posts = pd.cur_pp(curidx); 
              labeled = bwlabel(posts);
              aa = regionprops(labeled,'Area');  %#ok
              blen = [blen [aa.Area]];  %#ok
            end
          end
        end        
      else        
        % For loaded scores.
        for endx = 1:obj.nexps
          for flies = 1:obj.nflies_per_exp(endx)
            pd = obj.predictdata{endx}{flies}(iCls);
            curidx = pd.loaded_valid;
            curt = pd.t(curidx);
            if any(curt(2:end)-curt(1:end-1) ~= 1)
              warning('JLabelData:bouts','Scores are not in order');
              return;
            end
            posts = pd.loaded_pp(curidx);
            labeled = bwlabel(posts);
            aa = regionprops(labeled,'Area');  %#ok
            blen = [blen [aa.Area]];  %#ok
          end
        end        
      end
    end

    
    % ---------------------------------------------------------------------
    function [labels,labeledscores,allScores,scoreNorm] = GetAllLabelsAndScores(obj,iCls)
      % labels: current labels for classifier iCls, row vector of -1/1 for no-beh/beh resp
      % labeledscores: row vector same size as labels
      % allScores: row vector, predicted scores for iCls over all exps/flies
      % scoreNorm: scalar
      
      %MERGESTUPDATED
      
      wd = obj.windowdata(iCls);      
      if isempty(wd.exp)
        labels = zeros(1,0);
        labeledscores = zeros(1,0);
      else
        curNdx = wd.labelidx_cur~=0;
        origlabels = wd.labelidx_cur(curNdx);
        
        lblPosNeg = obj.iCls2iLbl{iCls};
        assert(all(origlabels==lblPosNeg(1) | origlabels==lblPosNeg(2)));
        labels = ((origlabels==lblPosNeg(1))-0.5)*2;
        labeledscores = myBoostClassify(wd.X(curNdx,:),obj.classifier{iCls});
      end
      
      allScores = zeros(1,0);
      for expi = 1:obj.nexps
        for flies = 1:obj.nflies_per_exp(expi)
          pd = obj.predictdata{expi}{flies}(iCls);
          curidx = pd.cur_valid;
          allScores = [allScores pd.cur(curidx)]; %#ok<AGROW>
        end
      end
      scoreNorm = wd.scoreNorm;
    end
    
    
    %     % ---------------------------------------------------------------------
%     function expStats = GetExpStats(obj,expi)
%       % Calculates statistics such as number of labeled bouts, predicted bouts
%       % and change in scores.
%       
%       expStats.name = obj.expnames{expi};
%       expStats.nflies = obj.nflies_per_exp(expi);
%       expStats.nlabeledbouts = obj.labelstats(expi).nbouts_labeled;
%       expStats.nlabeledflies = obj.labelstats(expi).nflies_labeled;
%       
%       
%       if ~isempty(obj.predictdata.exp==expi)
%         expid = obj.predictdata.exp==expi;
%         expStats.nscoreframes = nnz(expid);
%         expStats.nscorepos = nnz(obj.predictdata.loaded(expid)>0);
% %         if ~isempty(obj.predictdata.classifierfilenames) && ...
% %             numel(obj.predictdata.classifierfilenames)>=expi
% %           expStats.classifierfilename = obj.predictdata.classifierfilenames{expi};
% %         else
% %           expStats.classifierfilename = '';
% %         end
%       else
%         expStats.nscoreframes = [];
%         expStats.nscorefrac = [];
%         expStats.classifierfilename = '';
%       end
%       
%     end

    %    % ---------------------------------------------------------------------
%    function InitPostprocessparams(obj)
%      % AL no callsites
%      obj.postprocessparams.method = 'Hysteresis';
%      obj.postprocessparams.hystopts(1) = struct('name','High Threshold','tag','hthres','value',0);
%      obj.postprocessparams.hystopts(2) = struct('name','Low Threshold','tag','lthres','value',0);
%      obj.postprocessparams.filtopts(1) = struct('name','Size','tag','size','value',1);
%      obj.postprocessparams.blen = 1;
%    end

  end
  
  methods (Access=private)
    
    function [scores,predictions] = GetScoresCore(obj,expi,flies,...
        fld,fldValid,fldPP,T0,T1,clsIdx)
      % clsIdx: array of classifier indices. Defaults to 1:obj.nclassifiers.
      %
      % scores: numel(clsIdx)-by-(T1-T0+1) 
      % predictions: numel(clsIdx)-by-(T1-T0+1). 1/2 for 
      % behavior/no-behavior, resp. fldPP only needed to compute 
      % predictions.
      
      if ~exist('clsIdx','var')
        clsIdx = 1:obj.nclassifiers;
      end

      n = T1-T0+1;
      off = 1-T0;
      pdArr = obj.predictdata{expi}{flies};
      nCls = numel(clsIdx);

      scores = zeros(nCls,n); % AL: nan instead?
      predictions = zeros(nCls,n); % AL: nan instead?
      for i = 1:nCls
        iTL = clsIdx(i);
        pd = pdArr(iTL);
        idxcurr = pd.(fldValid) & pd.t>=T0 & pd.t<=T1;
        scores(i,pd.t(idxcurr)+off) = pd.(fld)(idxcurr);
        predictions(i,pd.t(idxcurr)+off) = 2-pd.(fldPP)(idxcurr);
      end
    end
  end
  
  %% Save/Load/Import/Export
  
  methods 
    
    % ---------------------------------------------------------------------
    function openJabFile(self, ...
                         fileNameAbs, ...
                         groundTruthingMode, ...
                         varargin)
      % openJabFile(self,fileNameAbs,gtMode,p1,v1,...)
      %
      % Optional PVs:
      %   * macguffin. Macguffin object, jab contents corresponding to
      %   fileNameAbs. This option exists because the caller may want to
      %   modify the project before opening. If this argument is provided,
      %   self.needsave is set to true.
      %   * originalExpDirNames. see below.
      %   * substituteExpDirNames. etc
      % 
      % original/substituteExpDirNames should be cell arrays of the same 
      % length, each element a string giving an absolute path to an 
      % experiment directory.  Each element of originalExpDirNames should
      % be an experiment directory in the .jab file, and the corresponding
      % element of substituteExpDirNames gives an experiment dir name to be
      % used in place of the original one. This is to enable the user to
      % manually locate exp dir names that are missing. Elements of
      % substituteExpDirNames may be empty (eg '') to indicate that the
      % experiment is not to be opened/included in the project. Whether
      % these experiment dir names are treated as normal exp dir names or
      % ground-truthing exp dir names depends on groundTruthingMode.      
                       
      %MERGESTUPDATED
      
      [macguffin,...
       originalExpDirs,...
       substituteExpDirs] = myparse(varargin,...
       'macguffin',[],...
       'originalExpDirs',cell(0,1),...
       'substituteExpDirs',cell(0,1));
                       
      assert(numel(originalExpDirs)==numel(substituteExpDirs));
      macgufProvided = ~isempty(macguffin);
      if macgufProvided
        assert(isa(macguffin,'Macguffin'));
      else
        self.SetStatus('Loading %s',fileNameAbs);
        macguffin = loadAnonymous(fileNameAbs);
        if isstruct(macguffin)
          if isfield(macguffin,'fromAPT') && macguffin.fromAPT
            macguffin = Macguffin(macguffin,macguffin.aptInfo);

            % keep a copy of info as they will get removed when Macguffin
            % is called.
%             behaviorName=macguffin.behaviors.names;
%             movieFileName=macguffin.file.moviefilename;
%             movieIndexFileName=macguffin.file.movieindexfilename;
%             trackFileName=macguffin.file.trxfilename;
%             scoreFileName=macguffin.file.scorefilename;
%             macguffin.behaviors.names=behaviorName;
%             macguffin.file.moviefilename=movieFileName;
%             macguffin.file.movieindexfilename=movieIndexFileName;
%             macguffin.file.trxfilename=trackFileName;
%             macguffin.file.scorefilename=scoreFileName;
          else
            macguffin = Macguffin(macguffin);
          end
        end
        self.ClearStatus();
      end
      macguffin.modernize(true);
      
      self.gtMode = groundTruthingMode;
      
      % Do the substiutions, if any
      substitutionsMade = false;
      if groundTruthingMode
        expDirNames = macguffin.gtExpDirNames;
        labels = macguffin.gtLabels;
      else
        expDirNames = macguffin.expDirNames;
        labels = macguffin.labels;
      end
      newExpDirNames = cell(1,0);
      newLabels = Labels.labels(0);
      removedDirs = [];
      newExpNumbers = []; count = 0;
      for i = 1:length(expDirNames)
        expDirName = expDirNames{i};
        j = whichstr(expDirName,originalExpDirs);
        if isempty(j)
          newExpDirNames{end+1} = expDirNames{i};  %#ok
          newLabels(end+1) = labels(i);  %#ok
          count = count+1;
          newExpNumbers(i) = count;
        else
          if ~isempty(substituteExpDirs{j})
            newExpDirNames{end+1} = substituteExpDirs{j};  %#ok
            newLabels(end+1) = labels(i);  %#ok
            count = count+1;
            newExpNumbers(i) = count;
          else
            removedDirs(end+1) = i;
            newExpNumbers(i) = 0;
            % Empty new name for this experiment; experiment will not be added
          end
          substitutionsMade = true;
        end
      end
      if groundTruthingMode
        macguffin.gtExpDirNames = newExpDirNames;
        macguffin.gtLabels = newLabels;
        
        % Update GTSuggestions
        % Not tested -- Mayank Dec 18 2015
        if numel(removedDirs)>0
          if ~isempty(macguffin.gtSuggestions.randomGTSuggestions)
            % only remove dirs on which there were random suggestions.
            nexp_suggestions = length(macguffin.gtSuggestions.randomGTSuggestions);
            curdirs = removedDirs;
            curdirs(curdirs>nexp_suggestions) = [];
            macguffin.gtSuggestions.randomGTSuggestions(curdirs) = [];
          end
          if ~isempty(macguffin.gtSuggestions.loadedGTSuggestions) 
            nexp_suggestions = length(macguffin.gtSuggestions.loadedGTSuggestions);
            curdirs = removedDirs;
            curdirs(curdirs>nexp_suggestions) = [];
            macguffin.gtSuggestions.loadedGTSuggestions(curdirs) = [];
          end
          if ~isempty(macguffin.gtSuggestions.balancedGTSuggestions)
            expidx = [macguffin.gtSuggestions.balancedGTSuggestions.exp];
            toremove = ismember(expidx,removedDirs);
            macguffin.gtSuggestions.balancedGTSuggestions(toremove) = [];
            for ndx = 1:numel(macguffin.gtSuggestions.balancedGTSuggestions)
              macguffin.gtSuggestions.balancedGTSuggestions(ndx).exp = ...
                newExpNumbers(macguffin.gtSuggestions.balancedGTSuggestions(ndx).exp);
            end
          end
        end
        
      else
        macguffin.expDirNames = newExpDirNames;
        macguffin.labels = newLabels;
      end
      
      % Set the JLD to match the Macguffin
      self.setMacguffin(macguffin,true);
      
      
      % Store file-related stuff
      self.thereIsAnOpenFile = true;
      self.everythingFileNameAbs = fileNameAbs;
      self.userHasSpecifiedEverythingFileName = true;
      self.needsave = macgufProvided  || substitutionsMade;
      self.defaultpath = fileparts(fileNameAbs);
     
      % initialize the status table describing what required files exist
      [success,msg] = self.UpdateStatusTable();
      if ~success,
        error('JLabelData:unableToUpdateStatusTable',msg);
      end
      
      % Load windowdata if appropriate
      cs = macguffin.classifierStuff;
      assert(isequal(self.nclassifiers,numel(cs),numel(self.windowdata)));
      perframeNdx = find(strcmp('perframedir',self.filetypes));
      for iCls = 1:self.nclassifiers
        if ~substitutionsMade && ...
            ~self.IsGTMode() && ... 
            self.loadwindowdata(iCls) && ...
            isprop(cs(iCls),'windowdata') && ...
            isstruct(cs(iCls).windowdata) && ...
            ~isempty(cs(iCls).windowdata) && ...
            ~isempty(cs(iCls).savewindowdata) && ...
            ~isempty(cs(iCls).windowdata) && ...
            cs(iCls).savewindowdata
          
          % determine whether to load windowdata for this classifier
          tfLoadWinData = true;
          if self.isInteractive
            isPerframeNewer = false;
            for ndx = 1:self.nexps
              if self.classifierTS(iCls) < self.filetimestamps(ndx,perframeNdx);
                isPerframeNewer = true;
                expnamenewer = self.expnames{ndx};
                perframeTS = self.filetimestamps(ndx,perframeNdx);
                break;
              end
            end
            
            if isPerframeNewer
              qstr{1} = sprintf('One of the perframe files (Generated on %s) ',...
                datestr(perframeTS));
              qstr{end+1} = sprintf('is newer than the classifier (Trained on %s)',datestr(self.classifierTS(iCls))); %#ok<AGROW>
              qstr{end+1} = sprintf('for the experiment %s.',expnamenewer); %#ok<AGROW>
              qstr{end+1} = ' Still load the windowdata stored in the jab file?'; %#ok<AGROW>
              res = questdlg(qstr, ...
                'Load Window Data?', ...
                'Yes','No', ...
                'No');
              tfLoadWinData = strcmpi(res,'Yes');
            end
          end
          % If ~self.isInteractive, or classifiers newer than all PF
          % dirs, then tfLoadWinData true by default
          
          if tfLoadWinData
            oldScoreNorm = self.windowdata(iCls).scoreNorm;
            oldfeaturenames = self.windowdata(iCls).featurenames;
            self.windowdata(iCls) = cs(iCls).windowdata;
            if isempty(self.windowdata(iCls).scoreNorm) && ~isempty(oldScoreNorm)
              self.windowdata(iCls).scoreNorm = oldScoreNorm;
            end
            if isempty(self.windowdata(iCls).featurenames) && ~isempty(oldfeaturenames)
              self.windowdata(iCls).featurenames = oldfeaturenames;
            end
          end
        end
      end
    end
    
    
    % ---------------------------------------------------------------------
    function openJabFileNoExps(self, ...
        fileNameAbs, ...
        groundTruthingMode)
      
      self.gtMode = groundTruthingMode;
      
      %
      % Open the file
      %
      
      % load the file
      macguffin = loadAnonymous(fileNameAbs);
      % if we get here, file was read successfully
      if isstruct(macguffin)
        if isfield(macguffin,'fromAPT') && macguffin.fromAPT
          macguffin = Macguffin(macguffin,macguffin.aptInfo);
        else
          macguffin = Macguffin(macguffin);
        end
      end
      macguffin.modernize(true);
      
      % Set the JLD to match the Macguffin
      self.setMacguffin(macguffin,false);
      
      % Store file-related stuff
      self.thereIsAnOpenFile = true;
      self.everythingFileNameAbs = fileNameAbs;
      self.userHasSpecifiedEverythingFileName = true;
      self.needsave = false;
      fileDirPathAbs = fileparts(fileNameAbs);
      self.defaultpath = fileDirPathAbs;
      
      % initialize the status table describing what required files exist
      [success,msg] = self.UpdateStatusTable();
      if ~success,
        error('JLabelData:unableToUpdateStatusTable',msg);
      end
    end  % method


    % ---------------------------------------------------------------------
    function newJabFile(obj,macguf)
      % Only called by ProjectSetup/new project creation. 
      % IMPORTANT: macguf has type Macguffin, but it is not a properly
      % initialized Macguffin object. It is semi-initialized object 
      % originating from ProjectSetup for the purposes of initialization, 
      % hence the various massaging here.
            
      obj.gtMode = false;
                     
      obj.setMacguffin(macguf);
      
      % Make up filename
      try
        realbehnames = Labels.verifyBehaviorNames(macguf.behaviors.names);
        fileNameRel = [sprintf('%s_',realbehnames{1:end-1}) realbehnames{end} '.jab'];
      catch excp
        if isequal(excp.identifier,'Macguffin:mainBehaviorNotDefined')
          fileNameRel = 'untitled.jab';
        else
          rethrow(excp);
        end
      end  
      fileNameAbs = fullfile(obj.defaultpath,fileNameRel);

      % Set other file-related instance vars
      obj.thereIsAnOpenFile = true;
      obj.everythingFileNameAbs = fileNameAbs;
      obj.userHasSpecifiedEverythingFileName = false;
      obj.needsave = true; % b/c new file
      
      % initialize the status table describing what required files exist
      [success,msg] = obj.UpdateStatusTable();
      if ~success,
        error(msg);
      end      
    end

    
    % ---------------------------------------------------------------------
    function closeJabFile(self)
      % The list of things we want to be persistent
      listOfPersistentSlots={'defaultpath' ...
                             'expdefaultpath' ...
                             'setstatusfn' ...
                             'clearstatusfn' ...
                             'cacheSize' ...
                             'version' ...
                             'perframeGenerate' ...
                             'perframeOverwrite' ...
                             'isInteractive'}';
      
      % Save the things we want to persist after closing the file
      nPersistentSlots=length(listOfPersistentSlots);
      for i=1:1:nPersistentSlots
        thisSlot=listOfPersistentSlots{i};
        evalString=sprintf('%s=self.%s;',thisSlot,thisSlot);
        eval(evalString);
      end
      
      % Nuke the site from orbit.  It's the only way to be sure.
      self.initialize();
      
      % Re-load the things we want to persist
      for i=1:1:nPersistentSlots
        thisSlot=listOfPersistentSlots{i};
        evalString=sprintf('self.%s=%s;',thisSlot,thisSlot);
        eval(evalString);
      end
    end
    
    
    % ---------------------------------------------------------------------
    function saveJabFile(self,fileNameAbs)
      fileNameRel = fileNameRelFromAbs(fileNameAbs);
      self.SetStatus(sprintf('Saving to %s...',fileNameRel));
      % Extract the structure that will be saved in the everything file
      macguffin = self.getMacguffin();
     
      try
        if exist(fileNameAbs,'file')
          % Warn if overwriting with updated/modernized jab
          oldQ = loadAnonymous(fileNameAbs);
          oldVer = oldQ.version;
          curVer = self.version; 
          if Jab.formatChangedBetweenVersions(oldVer,curVer)
            queststr = sprintf('The file ''%s'' will be saved in the latest format and will not be readable by some older versions of JAABA.',fileNameRel);
            btn = questdlg(queststr,...
              'Update File Format','OK, Continue','Cancel','OK, Continue');
            switch btn
              case 'OK, Continue'
              otherwise
                self.ClearStatus();
                return;
            end
          end
          
          % back up old jab anyway
          backupFileNameAbs = [fileNameAbs,'~'];
          [success,message,identifier] = copyfile(fileNameAbs,backupFileNameAbs); %#ok
          if ~success,
            backupFileNameRel = fileNameRelFromAbs(backupFileNameAbs);
            warning('JLabelData:unableToCreateBackup', ...
                  'Unable to create backup file %s.',backupFileNameRel);
          end
        end
        
        old = warning('query','MATLAB:structOnObject');
        warning('off','MATLAB:structOnObject'); % turn off annoying warning
        macguffinStruct = struct(macguffin);
        warning(old); % restore annoying warning
        saveAnonymous(fileNameAbs,macguffinStruct);
      catch excp
        self.ClearStatus();
        rethrow(excp);
      end
      
      % Do follow-up book-keeping
      self.everythingFileNameAbs=fileNameAbs;
      self.userHasSpecifiedEverythingFileName=true;      
      self.needsave=false;
      fileDirPathAbs=fileparts(fileNameAbs);
      self.defaultpath=fileDirPathAbs;      
      self.ClearStatus();
    end  % method
    
    
    % ---------------------------------------------------------------------
    function importClassifier(self,fileNameAbs)
      %MERGEST SEEMSOK 
      
      macguffin = loadAnonymous(fileNameAbs);

      self.setScoreFeatures(macguffin.scoreFeatures);
      self.setFeatureSublexicon(macguffin.featureLexicon, ...
                                macguffin.featureLexiconName, ...
                                macguffin.sublexiconPFNames);
      self.setWindowFeaturesParams(macguffin.windowFeaturesParams);
      
%       % Generate the necessary files now, so that any problems occur now.
%       for iExp=1:self.nexps
%         [success,msg]=self.GenerateScoreFeaturePerframeFiles(iExp);
%         if ~success,
%           error('JLabelData:unableToGenerateScoreFeaturePerframeFile',msg);
%         end
%         self.UpdateStatusTable('perframedir',iExp);
%         allPerframeFilesExist=self.fileexists(iExp,whichstr('perframedir',self.filetypes));
%         if ~allPerframeFilesExist , 
%           [success,msg]=self.GeneratePerframeFilesExceptScoreFeatures(iExp);
%           if ~success,
%             error('JLabelData:unableToGeneratePerframeFile',msg);
%           end
%         end
%       end
      
      % Load the classifier proper, training params, etc.
      self.setClassifierStuff(macguffin.classifierStuff);

%       % do this to prompt loading of windowdata
%       force=true;
%       self.setCurrentTarget(self.expi,self.flies,force);
      
      self.needsave=true;
    end
    
    
    % ---------------------------------------------------------------------
    function cs = getClassifierStuff(self)
      
      %MERGEST UPDATED
        
      % make sure current labels are committed
      self.StoreLabelsForCurrentAnimal();
      
      cs = ClassifierStuff.empty(0,1);
      for iCls = self.nclassifiers:-1:1
        csArgs = { ...
            'type',self.classifiertype{iCls}, ...
            'params',self.classifier{iCls}, ...
            'trainingParams',self.classifier_params{iCls}, ...
            'timeStamp',self.classifierTS(iCls), ...
            'confThresholds',self.confThresholds(iCls,:), ...
            'scoreNorm',self.windowdata(iCls).scoreNorm, ...
            'postProcessParams',self.postprocessparams{iCls}, ...
            'featureNames',self.windowdata(iCls).featurenames,...
            'savewindowdata',self.savewindowdata(iCls),...
            'predictOnlyCurrentFly',self.predictOnlyCurrentFly(iCls),...
            'selFeatures',self.selFeatures(iCls),...
            };
        if self.savewindowdata(iCls) && ~self.IsGTMode()
          csArgs(end+1:end+2) = {'windowdata',self.windowdata(iCls)};
        end          
        cs(iCls,1) = ClassifierStuff(csArgs{:});
      end
    end
    
    
    % ---------------------------------------------------------------------
    function setClassifierStuff(self,classifierStuff)
      % classifierStuff: array of ClassifierStuffs
      %
      % ALTODO: Seems like classifier-related properties can be
      % consolidated/cleaned up. Why not just use ClassifierStuff?
      
      %MERGESTUPDATED
      
      classifierStuff.modernize();
      
      nrealbeh = self.ntimelines;
      assert(numel(classifierStuff)==nrealbeh);
      self.classifiertype = cell(1,nrealbeh);
      self.classifier = cell(1,nrealbeh);
      self.classifier_params = cell(1,nrealbeh);
      self.classifierTS = nan(1,nrealbeh);
      self.confThresholds = nan(nrealbeh,2);
      self.postprocessparams = cell(1,nrealbeh);
      self.savewindowdata = false(1,nrealbeh);
      for iBeh = 1:nrealbeh
        cs = classifierStuff(iBeh);
        self.classifiertype{iBeh} = cs.type;
        self.classifier{iBeh} = cs.params;
        self.classifier_params{iBeh} = cs.trainingParams;
        self.classifierTS(iBeh) = cs.timeStamp;
        self.confThresholds(iBeh,:) = cs.confThresholds;
        self.windowdata(iBeh).scoreNorm = cs.scoreNorm;
        self.postprocessparams{iBeh} = cs.postProcessParams;
        self.savewindowdata(iBeh) = cs.savewindowdata;
        self.selFeatures(iBeh) = cs.selFeatures;
        self.predictOnlyCurrentFly(iBeh) = cs.predictOnlyCurrentFly;
      end
            
      % AL20141122: Why are we setting windowfeaturenames? This only
      % depends on self.curperframefns and self.windowfeaturescellparams.
      % MK20150805: We don't need do this again. It is already being done
      % in setWindowFeaturesParams.
%       self.SetWindowFeatureNames();
      
      % AL20141126: initialization of loadwindowdata is change from earlier
      % behavior
      self.loadwindowdata = true(1,nrealbeh);
      
      % verify windowFeatureNames in classifierStuff
      for iBeh = 1:nrealbeh
        cs = classifierStuff(iBeh);
        % KB 20201208
        % changing featureNames to a dependent variable, try not to create
        % it many times
        featureNames = cs.featureNames;
        if ~isempty(featureNames) && ...
           ~isempty(featureNames{1}) && ...
           ~isequal(featureNames,self.windowdata(iBeh).featurenames)
          warnstr = sprintf('The feature names stored in the jab file don''t match the current feature names. The loaded classifier ''%s'' shouldn''t be used; retrain a new classifier.',...
            self.labelnames{iBeh});
          if self.isInteractive ,
            uiwait(warndlg(warnstr));
          else
            warning(warnstr) ;  %#ok<SPWRN>
          end
          self.loadwindowdata(iBeh) = false;
        end
      end
            
      self.trainstats = cell(1,nrealbeh);

      % Update the window data near the labels
%       [success,msg] = self.PreLoadPeriLabelWindowData();
%       if ~success,error(msg);end   

      % Move the current predictions out of the way
      self.PredictDataMoveCurToOld();
      
      if ~self.isST
        
        % Set up for fast prediction
        self.FindFastPredictParams();

        % predict for all loaded examples
        self.PredictLoaded();

      end
    end  % setClassifierStuff() method
    
    
    % ---------------------------------------------------------------------
    function everythingParams = getMacguffin(self)
      % Construct the object that will be saved in the everything file
      everythingParams = Macguffin(self);
    end

    
    % ---------------------------------------------------------------------
    function setMacguffin(obj,everythingParams,loadexps)
      % This initializes the JLabelData object based on the contents of
      % everythingParams
  
      % Deal with arguments
      if ~exist('loadexps','var')
        loadexps=true;
      end
      
      before=obj.copy();  % make a copy of the object, in case something goes wrong
      try
        % Note sure what to do here---Macguffin class doesn't have a perframe
        % property at present
        if isfield(everythingParams.extra,'perframe'),
          if isfield(everythingParams.extra.perframe,'params') && isstruct(everythingParams.extra.perframe.params),
            pf_fields = fieldnames(everythingParams.extra.perframe.params);
            for ndx = 1:numel(pf_fields),
              everythingParams.featureLexicon.perframe_params.(pf_fields{ndx}) = everythingParams.extra.perframe.params.(pf_fields{ndx});
            end
          end
          if isfield(everythingParams.extra.perframe,'landmarkParams'),
            obj.landmark_params = everythingParams.extra.perframe.landmarkParams;
          end
        end  % isfield(basicParams,'perframe'),

        % feature config file
        obj.setFeatureSublexicon(everythingParams.featureLexicon, ...
                                 everythingParams.featureLexiconName, ...
                                 everythingParams.sublexiconPFNames);
%         if isequal(everythingParams.featureLexiconName,'custom')
%           %obj.setFeatureLexiconAndTargetSpeciesCustom(everythingParams.featureLexicon, ...
%           %                                            everythingParams.behaviors.type);
%           [success,msg] = obj.setFeatureLexiconAndFLName(everythingParams.featureLexicon,'custom');
%           obj.targettype=everythingParams.behaviors.type;
%         else
%           [success,msg]=obj.setFeatureLexiconAndTargetSpeciesFromFLName(everythingParams.featureLexiconName);
%         end
%         if ~success , 
%           error('JLabelData:unableToSetFeatureSublexicon',msg);
%         end

        % Set the target species
        obj.targettype=everythingParams.behaviors.type;
                
        obj.trxGraphicParams=cookTrxGraphicParams(everythingParams.trxGraphicParams);
        obj.labelcolors=everythingParams.behaviors.labelcolors;
        obj.unknowncolor=everythingParams.behaviors.unknowncolor;

        %
        % load in the rest of the stuff, depending on the fields present
        %

        % read in behavior names
        if isfield(everythingParams.behaviors,'names'),
          obj.labelnames = everythingParams.behaviors.names;
          assert(iscell(obj.labelnames));
          %obj.nbehaviors = numel(obj.labelnames);
        else
          obj.labelnames = {'Behavior','None'};
        end
        assert(numel(obj.labelcolors)==3*numel(obj.labelnames));

        [obj.ntimelines,obj.iLbl2iCls,obj.iCls2iLbl] = Labels.determineNumTimelines(obj.labelnames);

        % Do APT Stuff before the filetype list is used for everything else
        if everythingParams.fromAPT
          if ~isempty(everythingParams.aptInfo.trkfilename),
            [success1,msg] = obj.SetTrkFileName(everythingParams.aptInfo.trkfilename);
            if ~success1
              error('JLabelData:unableToSetTrkFileName', ...
                    msg);
            end
          end
          obj.aptInfo = everythingParams.aptInfo;
          obj.fromAPT = true;
        else
          % remove trk from file types
          ftypes = obj.filetypes;
          for ndx = numel(ftypes):-1:1
            if strcmp(ftypes{ndx},'trk')
              ftypes(ndx)=[];
            end
          end
          obj.filetypes = ftypes;
          obj.fileexists = false(0,numel(obj.filetypes));
          obj.aptInfo = struct;
          obj.fromAPT = false;
        end
        
        obj.stInfo = everythingParams.stInfo;
        obj.stFeatures = everythingParams.stFeatures;
        
        if isfield(everythingParams.file,'moviefilename'),
          if isfield(everythingParams.file,'movieindexfilename'),
            [success1,msg] = obj.SetMovieFileName(everythingParams.file.moviefilename,...
              everythingParams.file.movieindexfilename);
          else
            [success1,msg] = obj.SetMovieFileName(everythingParams.file.moviefilename);
          end
            
          if ~success1,
            error('JLabelData:unableToSetMovieFileName', ...
                  msg);
          end
        end
        if isfield(everythingParams.file,'trxfilename'),
          [success1,msg] = obj.SetTrxFileName(everythingParams.file.trxfilename);
          if ~success1,
            error('JLabelData:unableToSetTrxFileName', ...
                  msg);
          end
        end
        if isfield(everythingParams.file,'scorefilename'),
          scorefilename = everythingParams.file.scorefilename;
        else
          scorefilename = {ScoreFile.defaultScoreFilename(obj.labelnames{1})};
        end
        [success1,msg] = obj.setScoreFileName(scorefilename);
        if ~success1,
          error('JLabelData:unableToSetScoreFileName', ...
                msg);
        end
        if isfield(everythingParams.file,'perframedir'),
          [success1,msg] = obj.SetPerFrameDir(everythingParams.file.perframedir);
          if ~success1,
            error('JLabelData:unableToSetPerframeDirName', ...
                  msg);
          end
        end
        if isfield(everythingParams.file,'clipsdir') && ~isempty(everythingParams.file.clipsdir),
          [success1,msg] = obj.SetClipsDir(everythingParams.file.clipsdir);
          if ~success1,
            error('JLabelData:unableToSetClipsDirName', ...
                  msg);
          end
        end
                
        obj.stfeatures = 'features.mat'; % temporary hardcode
        [success1,msg] = obj.UpdateStatusTable('stfeatures');
        if ~success1
          error('JLabelData:unableToSetSTFeaturesName',msg);
        end

        if isfield(everythingParams.extra,'usePastOnly'),
          obj.usePastOnly = everythingParams.extra.usePastOnly;
        end

        if isprop(everythingParams,'gtSuggestions') && ...
            isfield(everythingParams.gtSuggestions,'GTSuggestionMode')
          obj.GTSuggestionMode = everythingParams.gtSuggestions.GTSuggestionMode;
          obj.randomGTSuggestions = everythingParams.gtSuggestions.randomGTSuggestions;
          obj.thresholdGTSuggestions = everythingParams.gtSuggestions.thresholdGTSuggestions ;
          obj.loadedGTSuggestions = everythingParams.gtSuggestions.loadedGTSuggestions;
          obj.balancedGTSuggestions = everythingParams.gtSuggestions.balancedGTSuggestions ;
        end

  %       if isfield(basicParams,'scoreFeatures') ,
  %         obj.scoreFeatures = basicParams.scoreFeatures;
  %         nScoreFeaturess=length(basicParams.scoreFeatures);
  %         scoreFeaturesPFNames=cell(nScoreFeaturess,1);
  %         for i = 1:nScoreFeaturess ,
  %           [~,pfName] = fileparts(obj.scoreFeatures(i).scorefilename);
  %           scoreFeaturesPFNames{i} = pfName;
  %         end
  %         obj.allperframefns=[obj.allperframefns ; ...
  %                             scoreFeaturesPFNames];
  %       end  % if isfield(basicParams,'scoreFeatures'),

        % Re-load the perframe feature signals, since the PFFs may have changed
        obj.loadPerframeData(obj.expi,obj.flies);

  %       % initialize the post-processing parameters
  %       obj.InitPostprocessparams();

        % initialize everything else
        obj.setWindowFeaturesParams(everythingParams.windowFeaturesParams);
        % AL 201506 order important; setWindowFeaturesParams provides initialization of eg .predictblocks which
        % is necessary should setAllLabels->SetExpDirs require a rollback
        if loadexps,
          obj.setAllLabels(everythingParams);
        end        
        obj.setScoreFeatures(everythingParams.scoreFeatures);
        obj.setClassifierStuff(everythingParams.classifierStuff);
        
      catch excp
        % If there's a problem, restore the object to its original state.
        obj.setToValue(before);
        rethrow(excp);
      end
        
    end  % method    

  
    % ---------------------------------------------------------------------
    function clearClassifierProper(self)
      % Reset the classifier to a blank slate
      
      self.classifier = struct('dim',{}, ...
                             'error',{}, ...
                             'dir',{}, ....
                             'tr',{}, ...
                             'alpha',{});  % 0x1 struct array
      nCls = self.nclassifiers;
      self.classifier = repmat({self.classifier},1,nCls); % ALTODO scattered classifier initialization: JLD, ProjectSetup
      self.classifier_old = self.classifier;
      self.classifierTS = zeros(1,nCls); 
      for iCls = 1:nCls
        self.windowdata(iCls).scoreNorm = nan;
      end
      self.PredictDataInvalidate();
      self.UpdatePredictedIdx(); % update cached predictions for current target
      self.needsave = true;
    end
    
    
    % ---------------------------------------------------------------------
    function basicParams = getBasicParamsStruct(obj)
      basicParams=struct();
      basicParams.featureLexiconName=obj.featureLexiconName;
        basicParams.scoreFeatures=obj.scoreFeatures;
      subdialectPFNames=obj.allperframefns;
      nScoreFeaturess=length(obj.scoreFeatures);
      sublexiconPFNames=subdialectPFNames(1:end-nScoreFeaturess);  
      basicParams.sublexiconPFNames=sublexiconPFNames;

      %assert(numel(basicParams.behaviors.names)==numel(basicParams.behaviors.labelcolors));
      %ALTODO labelcolors one less currently
      %assert(strcmp(obj.labelnames{end},'None'));
      basicParams.behaviors.names=obj.labelnames;
      basicParams.behaviors.labelcolors=obj.labelcolors;
      basicParams.behaviors.unknowncolor=obj.unknowncolor;
      
      basicParams.file.moviefilename=obj.moviefilename;
      basicParams.file.trxfilename=obj.trxfilename;
      basicParams.file.scorefilename=obj.scorefilename;
      %basicParams.scoresinput=obj.scoreFeatures;
      basicParams.trxGraphicParams=obj.trxGraphicParams;
    end
    

    function setsavewindowdata(self,value)
      assert(numel(self.savewindowdata)==self.nclassifiers);
      assert(isscalar(value) || numel(value)==self.nclassifiers);
      self.savewindowdata(:) = value;
      self.needsave = true;      
    end
  
    function setdofeatureselection(self)
      assert(numel(self.selFeatures)==self.nclassifiers);
      for ndx = 1:numel(self.selFeatures),
        self.selFeatures(ndx).do = true;
      end
    end
     
    function setusefeatureselection(self,value)
      assert(numel(self.selFeatures)==self.nclassifiers);
      assert(isscalar(value) || numel(value)==self.nclassifiers);
      for ndx = 1:numel(self.selFeatures),
        self.selFeatures(ndx).use = value;
        if value,
          self.selFeatures(ndx).do = true;
        end
      end
      self.needsave = true;      
    end
    
    function setPredictOnlyCurFly(self,value)
      self.predictOnlyCurrentFly(:) = value;
      self.needsave = true;
    end
    
    function value = getPredictOnlyCurFly(self)
      if self.nclassifiers >= 1
        if ~all(self.predictOnlyCurrentFly==self.predictOnlyCurrentFly(1))
          warndlg('Flag for Predicting only on current fly should be same for all classifiers','Flag mismatch');
        end
        value = self.predictOnlyCurrentFly(1);
      else
        value = false;
      end
    end
    
    function setClassifierParams(self,params)
      %MERGESTUPDATED
      
      assert(iscell(params) && numel(params)==self.nclassifiers);
      oldParams = self.classifier_params;      
      
      for iCls = 1:self.nclassifiers
        prm0 = oldParams{iCls};
        prm1 = params{iCls};
        if prm0.numBins~=prm1.numBins
          self.windowdata(iCls).binVals = [];
          self.windowdata(iCls).bins = [];
        end
      end
      self.classifier_params = params(:)';
    end
    
    % AL 20141210 appears unused, and see below for a commented (identical?) method  
%     % ---------------------------------------------------------------------
%     function SaveCurScores(self,expi,sfn)
%     % Saves the current scores to a file.
%         
%       if nargin < 3
%         sfn = self.GetFile('scores',expi);
%       end
%     
%       if ~self.HasCurrentScores(),
%         %uiwait(warndlg('No scores to save'));
%         return
%       end
%       
%       allScores = struct('scores',{{}},'tStart',[],'tEnd',[],...
%                          'postprocessed',{{}},'postprocessparams',[]);
%       scores_valid = true;
%       pdExp = self.predictdata{expi};
%       nFly = self.nflies_per_exp(expi);
%       firstFrms = self.firstframes_per_exp{expi};
%       endFrms = self.endframes_epr_exp{expi};      
%       assert(iscell(pdExp) && numel(pdExp)==nFly);
%       assert(isequal(nFly,numel(firstFrms),numel(endFrms)));
%       
%       for fly = 1:nFly
%         
%       
%       end
%       
%       if ~scores_valid,
%         % uiwait(warndlg(['Scores have not been computed for all the frames for experiment ' ...
%         %  '%s. Cannot save the scores.'],self.expnames{expi}));
%         % return;
%         error('JLabelData.scoresHaveNotBeenComputed', ...
%               ['Scores have not been computed for all the frames of experiment ' ...
%                '%s. Cannot save the scores.'],self.expnames{expi});  %#ok
%       end
%       allScores.postprocessedparams = self.postprocessparams;
%       allScores.scoreNorm = self.windowdata.scoreNorm;
%       self.SaveScores(allScores,sfn);      
%     end  % method

%     % ---------------------------------------------------------------------
%     function [success,msg] = SetClassifierFileName(self,classifierfilename,varargin)
%     % [success,msg] = SetClassifierFileName(obj,classifierfilename)
%     % Sets the name of the classifier file. If the classifier file exists, 
%     % it loads the data stored in the file. This involves removing all the
%     % experiments and data currently loaded, setting the config file,
%     % setting all the file names set in the config file, setting the
%     % experiments to be those listed in the classifier file, clearing all
%     % the previously computed window data and computing the window data for
%     % all the labeled frames. 
%       [classifierlabels,doreadconfigfile] = ...
%         myparse(varargin,...
%                 'classifierlabels',false,...
%                 'doreadconfigfile',true);
%       success = false;
%       msg = '';
%       self.classifierfilename = classifierfilename;
%       if ~isempty(classifierfilename) && exist(classifierfilename,'file'),
%         classifierParams = load(self.classifierfilename);
%         [success,msg]= ...
%           self.setClassifierParamsOld(classifierParams, ...
%                                       'classifierfilename',classifierfilename, ...
%                                       'classifierlabels',classifierlabels, ...
%                                       'doreadconfigfile',doreadconfigfile);
%       end
%     end  % method

    
%     % ---------------------------------------------------------------------
%     function [success,msg] = setClassifierParamsOld(obj, ...
%                                                     classifierParams, ...
%                                                     varargin)
%       % [success,msg] = SetClassifierFileName(obj,classifierfilename)
%       % Sets the name of the classifier file. If the classifier file exists,
%       % it loads the data stored in the file. This involves removing all the
%       % experiments and data currently loaded, setting the config file,
%       % setting all the file names set in the config file, setting the
%       % experiments to be those listed in the classifier file, clearing all
%       % the previously computed window data and computing the window data for
%       % all the labeled frames.
%       
%       [classifierlabels,doreadconfigfile,classifierfilename] = ...
%         myparse(varargin,...
%                 'classifierlabels',false, ...
%                 'doreadconfigfile',true, ...
%                 'classifierfilename',0);
%     
%       success = false;  %#ok
%       msg = '';  %#ok
%       
%       obj.classifierfilename = classifierfilename;
% 
%       if ischar(classifierParams.labelfilename) && ~strcmp(classifierParams.labelfilename,obj.labelfilename),
%         success = false;
%         msg = sprintf(['Label file name specified for the project (%s) don''t match' ...
%           ' the label file name used to train the classifier (%s). Not loading the classifier'],...
%           obj.labelfilename,classifierParams.labelfilename);
%         return;
%       end
% 
%       if ischar(classifierfilename)
%         obj.SetStatus('Loading classifier from %s',obj.classifierfilename);
%       else
%         obj.SetStatus('Loading classifier...');
%       end
% 
%       % remove all experiments
%       obj.RemoveExpDirs(1:obj.nexps);
% 
%       if doreadconfigfile,
%         % set config file
%         %     if ~strcmp(obj.configfilename,'configfilename'),
%         %       obj.SetConfigFileName(classifierParams.configfilename);
%         %     end
% 
%         % set movie
%         [success,msg] = obj.SetMovieFileName(classifierParams.moviefilename);
%         if ~success,error(msg);end
% 
%         % trx
%         [success,msg] = obj.SetTrxFileName(classifierParams.trxfilename);
%         if ~success,error(msg);end
% 
%         % labelPreLoad
%         [success,msg] = obj.SetLabelFileName(classifierParams.labelfilename);
%         if ~success,error(msg);end
% 
%         % perframedir
%         [success,msg] = obj.SetPerFrameDir(classifierParams.perframedir);
%         if ~success,error(msg);end
% 
%         % clipsdir
%         [success,msg] = obj.SetClipsDir(classifierParams.clipsdir);
%         if ~success,error(msg);end
%       end
% 
%       % featureparamsfilename
% %           [success,msg] = obj.SetFeatureParamsFileName(classifierParams.featureparamsfilename);
% %           if ~success,error(msg);end
% 
%       % load actual window features params instead of filename.
%       if all( isfield(classifierParams,{'windowfeaturesparams','windowfeaturescellparams',...
%                                   'basicFeatureTable','maxWindowRadiusCommonCached'}))
% 
%         classifierParams.windowfeaturesparams = JLabelData.convertTransTypes2Cell(classifierParams.windowfeaturesparams);
%         classifierParams.windowfeaturescellparams = JLabelData.convertParams2CellParams(classifierParams.windowfeaturesparams);
%         if ~( isequal(obj.windowfeaturesparams,classifierParams.windowfeaturesparams) && ...
%                 isequal(obj.maxWindowRadiusCommonCached,classifierParams.maxWindowRadiusCommonCached)),
%             str = sprintf('Window feature parameters in the configuration file');
%             str = sprintf('%s\ndo not match the parameters saved in the classifier',str);
%             str = sprintf('%s\nUsing parameters stored in the classifier file',str);
%             uiwait(warndlg(str));
%             obj.setWindowFeaturesParams(classifierParams.windowfeaturesparams,...
%                                      classifierParams.basicFeatureTable,...
%                                      classifierParams.maxWindowRadiusCommonCached);
%         end
%       end
% 
%       if ~isfield(classifierParams,'featurenames')
%         feature_names = {};
%         for j = 1:numel(obj.curperframefns),
%           fn = obj.curperframefns{j};
%           [~,feature_names_curr] = ComputeWindowFeatures([0,0],...
%             obj.windowfeaturescellparams.(fn){:});
%           feature_names_curr = cellfun(@(x) [{fn},x],feature_names_curr,'UniformOutput',false);
%           feature_names = [feature_names,feature_names_curr]; %#ok<AGROW>
%         end
%         obj.windowdata.featurenames = feature_names;
%       else
%         obj.windowdata.featurenames = classifierParams.featurenames;
%       end
% 
% 
%       % rootoutputdir
% %           [success,msg] = obj.SetRootOutputDir(classifierParams.rootoutputdir);
% %           if ~success,error(msg); end
% 
%       % set experiment directories
%       if classifierlabels && isfield(classifierParams,'labels'),
%         [success,msg] = obj.SetExpDirs(classifierParams.expdirs,classifierParams.outexpdirs,...
%           classifierParams.nflies_per_exp,classifierParams.sex_per_exp,classifierParams.frac_sex_per_exp,...
%           classifierParams.firstframes_per_exp,classifierParams.endframes_per_exp);
%         if ~success,error(msg); end
%         obj.labels = classifierParams.labels;
%         [obj.labelidx,obj.t0_curr,obj.t1_curr] = obj.GetLabelIdx(obj.expi,obj.flies);
%         obj.labelidx_off = 1 - obj.t0_curr;
%         [success,msg] = obj.PreLoadPeriLabelWindowData();
%         if ~success,error(msg); end
%         obj.labelsLoadedFromClassifier = true;
%       else
%         if classifierlabels,
%           uiwait(warndlg('The classifier file didn''t have any labels. Loading the current labels'));
%         end
%         [success,msg] = obj.SetExpDirs(classifierParams.expdirs,classifierParams.outexpdirs,...
%           classifierParams.nflies_per_exp,classifierParams.sex_per_exp,classifierParams.frac_sex_per_exp,...
%           classifierParams.firstframes_per_exp,classifierParams.endframes_per_exp);
%         if ~success,error(msg); end
%       end
%       [success,msg] = obj.UpdateStatusTable();
%       if ~success, error(msg); end
% 
%       % update cached data
% %           obj.windowdata = struct('X',[],'exp',[],'flies',[],'t',[],...
% %             'labelidx_cur',[],'labelidx_new',[],'featurenames',{{}},...
% %             'predicted',[],'predicted_probs',[],'isvalidprediction',[]);
%       [success,msg] = obj.PreLoadPeriLabelWindowData();
%       if ~success,error(msg);end
% 
%       obj.classifier = classifierParams.classifier;
%       obj.classifiertype = classifierParams.classifiertype;
%       obj.classifierTS = classifierParams.classifierTS;
%       obj.windowdata.scoreNorm = classifierParams.scoreNorm;
%       obj.confThresholds = classifierParams.confThresholds;
%       if isfield(classifierParams,'postprocessparams')
%         obj.postprocessparams = classifierParams.postprocessparams;
%       end
% 
%       paramFields = fieldnames(classifierParams.classifier_params);
%       for ndx = 1:numel(paramFields)
%         obj.classifier_params.(paramFields{ndx}) = classifierParams.classifier_params.(paramFields{ndx});
%       end
%       % predict for all loaded examples
%       obj.PredictLoaded();
% 
%       % set labelidx_cur
%       obj.SetTrainingData(classifierParams.trainingdata);
% 
% %           if strcmp(obj.classifiertype,'boosting'),
% %             [obj.windowdata.binVals, obj.windowdata.bins] = findThresholds(obj.windowdata.X);
% %           end
% 
%       % make sure inds is ordered correctly
%       if ~isempty(obj.classifier),
%         switch obj.classifiertype,
% 
%           case 'ferns',
%             waslabeled = obj.windowdata.labelidx_cur ~= 0;
%             obj.classifier.inds = obj.predict_cache.last_predicted_inds(waslabeled,:);
% 
%         end
%       end
% 
%       % clear the cached per-frame, trx data
%       obj.ClearCachedPerExpData();
% 
% %         catch ME,
% %           errordlg(getReport(ME),'Error loading classifier from file');
% %         end
% 
%       obj.ClearStatus();
%       obj.classifierfilename = classifierfilename;
%       obj.FindFastPredictParams();
%     end  % setClassifierParamsOld() method

%     % ------------------------------------------------------------------------
%     function basicParams=basicParamsFromMacguffin(everythingParams)
%       basicParams=struct();
%       basicParams.featureLexiconName=everythingParams.featureLexiconName;
%       basicParams.featureLexicon=everythingParams.featureLexicon;
%       %basicParams.scoreFeatures=everythingParams.scoreFeatures;
%       basicParams.sublexiconPFNames=everythingParams.sublexiconPFNames;
%       basicParams.behaviors=everythingParams.behaviors;  % need the animal type, in case featureLexiconName is 'custom'
%       basicParams.behaviors.names=everythingParams.behaviors.names(1);  % just want the first one
%       basicParams.file=everythingParams.file;
%       basicParams.trxGraphicParams=everythingParams.trxGraphicParams;
%       basicParams.landmarkParams=everythingParams.landmarkParams;
%     end    

  end
  
  
  %% Evaluating performance
  
  methods

    % ---------------------------------------------------------------------
    function [success,msg,crossErrorCell] = CrossValidate(obj,varargin)
    % Cross validate on bouts.
    %
    % success: nclassifier-by-1 logical
    % msg: nclassifier-by-1 cellstr
    % crossErrorCell: nclassifier-by-1 cell
      
    %MERGESTUPDATED
    
      nCls = obj.nclassifiers;
      success = false(nCls,1);
      msg = repmat({''},nCls,1);
      crossErrorCell = cell(nCls,1);
      
      obj.StoreLabelsAndPreLoadWindowData();      
      
      tmpsuccess = obj.PreLoadPeriLabelWindowData();
      if ~tmpsuccess        
        return;
      end
      
      [setidx,byexp] = myparse(varargin,'setidx',[],'byexp',false);
      
      for iCls = 1:nCls
        
        wd = obj.windowdata(iCls);
        islabeled = wd.labelidx_new~=0 & wd.labelidx_imp;
        if ~any(islabeled)
          msg{iCls} = 'No Labeled Data';
          continue;
        end
        
        if ~strcmp(obj.classifiertype{iCls},'boosting'); 
          msg{iCls} = 'Not boosting classifier';
          continue;
        end        
        
        obj.SetStatus('Cross-validating classifier ''%s'' for %d examples...',...
          obj.classifiernames{iCls},nnz(islabeled));
        
        obj.UpdateBoostingBins(iCls);
        
        bouts = obj.getLabeledBouts(iCls);
        
        if byexp && isempty(setidx)
          [~,~,setidx] = unique(wd.exp);
        end
        
        % Don't use unimportant labels
        labels = wd.labelidx_new;
        labels(~wd.labelidx_imp) = 0;
        iLbl = obj.iCls2iLbl{iCls};
        iLblPos = iLbl(1);
        iLblNeg = iLbl(2);
        labels012 = Labels.labelVec2label012(labels,iLblPos,iLblNeg);
                
        [success(iCls),msg{iCls},crossScores] = ...
          obj.crossValidateBout(iCls,labels012,bouts,setidx);        
        if ~success(iCls)
          continue;          
        end
                
        assert(isrow(crossScores));
        obj.windowdata(iCls).scores_validated = zeros(numel(islabeled),1);
        obj.windowdata(iCls).scores_validated(islabeled) = crossScores(1,:);
        
        % AL: various versions of label-vectors here and above presumably
        % redundant
        labelsLabeled = wd.labelidx_new(islabeled);
        labelsImpLabeled = wd.labelidx_imp(islabeled);
        assert(all(labelsImpLabeled==1));
        labelsLabeled12 = Labels.labelVec2label012(labelsLabeled,iLblPos,iLblNeg);        
        assert(all(labelsLabeled12==1 | labelsLabeled12==2));
        % modLabels: 1==positive+important; 3==negative+important
        modLabels = 2*labelsLabeled12 - labelsImpLabeled;
        
        %crossError=zeros(1,size(crossScores,1));
        %nSomething = size(crossScores,1);
        %crossError = struct('numbers',[],'frac',[]);
        crossError = obj.createConfMat(iCls,crossScores,modLabels);
        
        waslabeled = false(numel(islabeled),1);
        waslabeled(1:numel(wd.labelidx_old)) = wd.labelidx_old~=0 & wd.labelidx_imp;
        oldSelect = waslabeled(islabeled);
        oldScores = crossScores(oldSelect);
        tfWasIsLbled = waslabeled(:)&islabeled(:);
        labelsCur = wd.labelidx_cur(tfWasIsLbled);
        labelsImp = wd.labelidx_imp(tfWasIsLbled);
        assert(all(labelsImp==1));
        labelsCur12 = Labels.labelVec2label012(labelsCur,iLblPos,iLblNeg);
        
        oldLabels = 2*labelsCur12 - labelsImp;
        oldError = obj.createConfMat(iCls,oldScores,oldLabels);
        crossError.oldNumbers = oldError.numbers;
        crossError.oldFrac = oldError.frac;
        
        crossErrorCell{iCls} = crossError;
      end
      
      obj.ClearStatus();
    end

    
    % ---------------------------------------------------------------------    
    function [success,msg,scores] = ...
        crossValidateBout(obj,iCls,labels,bouts,setidx)
      %
      % Let windowdata have length nsamp, ie nsamp = numel(obj.windowdata(iCls).t)
      %
      % iCls: classifier index
      % labels: vector of length nsamp. Values are 1/2 for pos/neg labels, 
      %   resp. Conceptually like obj.windowdata(iCls).labelidx; can 
      %   contain 0 (?).
      % bouts: see getLabeledBouts(). Used only for generating setidx when
      %   necessary.
      % setidx: optional, integer grouping vector, length nsamp. Values in 
      %   range 1:k. Partitions windowdata for k-fold crossvalidation. If 
      %   not provided, setidx is generated using cvpartition.
      %
      % scores: vector of classification scores. Length nsamp, except 
      %  indices for labels==0 removed.
      
      %MERGESTUPDATED
      
      data = obj.windowdata(iCls).X;
      binVals = obj.windowdata(iCls).binVals;
      bins = obj.windowdata(iCls).bins;
      params = obj.classifier_params{iCls};
      
      [nsamp,nftrs] = size(data);
      assert(size(binVals,2)==nftrs);
      assert(isvector(labels) && numel(labels)==nsamp);
      assert(all(labels==0 | labels==1 | labels==2)); % AL: are there labels==0?
      assert(size(bouts.ndx,2)==nsamp);      

      modLabels = sign((labels==1)-0.5); % labels==0?

      if ~exist('setidx','var') || isempty(setidx)
        assert(all(bouts.label==1 | bouts.label==2));
        posBouts = bouts.label==1;
        negBouts = ~posBouts;
        nFolds = params.CVfolds;
        
        numPosBouts = nnz(posBouts);
        numNegBouts = nnz(negBouts);        
        if numPosBouts<nFolds || numNegBouts<nFolds
          scores = zeros(1,nsamp);
          scores(labels==0) = [];
          success = false;          
          if numPosBouts<nFolds
            msg = 'Too few bouts of behavior to do cross-validation';
          end
          if numNegBouts<nFolds
            msg = 'Too few bouts of not-behavior to do cross-validation';
          end          
          return;
        end
        
        cvpart = cvpartition(posBouts,'kfold',nFolds);
        
        setidx = nan(1,nsamp);
        for iFold = 1:nFolds
          iBoutTest = test(cvpart,iFold); % should contain both pos and neg bouts
          iBoutTest = find(iBoutTest);
          testNdx = false(1,nsamp);
          for iB = iBoutTest(:)'
            testNdx = testNdx | bouts.ndx(iB,:);
          end
          assert(all(isnan(setidx(testNdx))),...
            'Overlapping bouts encountered for classifier %d.',iCls);
          setidx(testNdx) = iFold;
        end
      end
      
      validateattributes(setidx,{'numeric'},{'positive' 'integer' 'vector' 'numel' nsamp});
      k = max(setidx);
      assert(isequal(unique(setidx(:)'),1:k));
      
      % AL20150528: I may have encountered a situation where the 
      % crossvalidation results table did not sum up as expected across
      % rows (with totals corresponding to number of frames labeled etc). 
      % I am including this temporary debug statement so we may check
      % against the results table.
      fprintf(1,'Count of labels (pos/neg/tot): %d/%d/%d\n',...
        nnz(labels==1),nnz(labels==2),numel(labels));
      
      scores = zeros(1,nsamp);
      for bno = 1:k
        curTestNdx = setidx==bno;
        curTrainNdx = setidx~=bno;
        
        curTrainLabels = modLabels(curTrainNdx);
        wt = getWeights(curTrainLabels);
        tt = tic;
        curbins = curTrainNdx;
        [~,curModel] = loglossboostLearnRandomFeatures(data(curTrainNdx,:),curTrainLabels,...
          params.iter,wt,binVals,bins(:,curbins),params);
        tScores = myBoostClassify(data(curTestNdx,:),curModel);
        scores(1,curTestNdx) = tScores;
        etime = toc(tt);
        obj.SetStatus('%d%% cross-validation done.  Time Remaining: %d s ',...
          round(bno/k*100),round((k-bno)*etime));
        drawnow();
      end

      scores(:,labels==0) = [];
      success = true;
      msg = '';
    end
    
    
    % ---------------------------------------------------------------------
    function [curScores,modLabels] = getCurrentScoresForROCCurve(obj,iCls)
      %MERGESTUPDATED
      
      if ~obj.classifierIsPresent()
        error('JLabelData:noClassifier', ...
              'No classifier has been trained to set the confidence thresholds.');
      end
      
      clsfr = obj.classifier{iCls};
      windata = obj.windowdata(iCls);
      if isempty(windata.X),
        error('JLabelData:noWindowData',...
              'No window data exists for the current classifier. Please retrain the classifier');
      end
      labels = windata.labelidx_new;
      iLbls = obj.iCls2iLbl{iCls};
      
      curNdx = labels~=0;
      curLabels = labels(curNdx);
      curLabels12 = Labels.labelVec2label012(curLabels,iLbls(1),iLbls(2));
      modLabels = ((curLabels12==1)-0.5)*2;
      curScores = myBoostClassify(windata.X(curNdx,:),clsfr);
    end
        
    
    % ---------------------------------------------------------------------
    function newError = TestOnNewLabels(obj,iCls)
      %MERGESTUPDATED
      
      obj.StoreLabelsAndPreLoadWindowData();
      [success,msg] = obj.PreLoadPeriLabelWindowData();
      if ~success
        warning(msg);
        return;
      end
            
      newError = struct();
      
      windata = obj.windowdata(iCls);

      prevLabeled = windata.labelidx_cur~=0;
      Nprev = numel(prevLabeled);
      newLabels = windata.labelidx_new~=0;
      tOld = newLabels(1:Nprev);
      tOld(prevLabeled) = false;
      newLabels(1:Nprev) = tOld;
      
      if ~nnz(newLabels)
        fprintf('No new labeled data\n');
        return;
      end
      
      % Find out the index of scores with the same exp, flynum and time as
      % the newly labeled data.
      
      orderedScores = [];
      orderedLabels = [];
      orderedLabels_imp = [];
      nlexp = windata.exp(newLabels);
      nlflies = windata.flies(newLabels);
      nlt = windata.t(newLabels);
      nlLabels = windata.labelidx_new(newLabels);
      iLbls = obj.iCls2iLbl{iCls};
      nlLabels012 = Labels.labelVec2label012(nlLabels,iLbls(1),iLbls(2));
      nlLabels_imp = windata.labelidx_imp(newLabels);
      
      classifierfilename = 'None';
      %setClassifierfilename = 1;
      unNLExp = unique(nlexp);
      for curExp = unNLExp(:)'
        curNLexpNdx = nlexp==curExp;
        unFlies = unique(nlflies(curNLexpNdx));
        for curFly = unFlies(:)'
          tfExpFly = nlexp==curExp & nlflies==curFly;
          curT = nlt(tfExpFly);
          curLabels = nlLabels012(tfExpFly);
          curLabels_imp = nlLabels_imp(tfExpFly);
          
          pd = obj.predictdata{curExp}{curFly}(iCls);
          curScoreNdx = find(pd.cur_valid);
          scoresT = pd.t(curScoreNdx);
          [curValidScoreNdx,loc] = ismember(scoresT,curT);
          if nnz(curValidScoreNdx)~=numel(curT)
            warndlg('Scores are missing for some labeled data');
            newError = struct();
            return;
          end
          
          orderedLabels = [orderedLabels; curLabels(loc(loc~=0))];  %#ok
          orderedLabels_imp = [orderedLabels_imp; curLabels_imp(loc(loc~=0))];  %#ok
          orderedScores = [orderedScores; pd.cur(curScoreNdx(curValidScoreNdx~=0))'];  %#ok
        end
%         if setClassifierfilename,
%           classifierfilename = obj.windowdata.classifierfilenames{curExp};
%           setClassifierfilename = 0;
%         elseif strcmp(classifierfilename,'multiple'),
%         elseif ~strcmp(classifierfilename,obj.windowdata.classifierfilenames{curExp}),
%           classifierfilename = 'multiple';
%         end
          
      end
      
      modLabels = 2*orderedLabels - orderedLabels_imp;      
      newError = obj.createConfMat(iCls,orderedScores,modLabels);
      newError.classifierfilename = classifierfilename;      
    end
        
  end
  
  methods (Access=private)
    
    % ---------------------------------------------------------------------
    function errorRates = createConfMat(obj,iCls,scores,modLabels)
      % iCls: classifier index
      % scores: vector 
      % modLabels: label vector for scores with values as follows:
      %   1: lbl1+important
      %   2: lbl1+notimportant
      %   3: lbl2+important
      %   4: lbl2+notimportant
      % errorRates: .numbers, .frac: 4x3 arrays. rows are labeled:
      %   row1: lbl1+important
      %   row2: lbl1+any, ie lbl1+(important or notimportant)
      %   row3: lbl2+important
      %   row4: lbl2+any
      
      %MERGESTUPDATED
      
      % AL: This method is a little awkward in that the original intent 
      % appears to have been to handle more than two labels/behaviors;
      % for now we are explicitly specifying a classifier index and
      % assuming only two label/behavior types in modLabels.      
      
      assert(isvector(scores) && isvector(modLabels) ...
        && numel(scores)==numel(modLabels));
      
      scoreNorm = obj.windowdata(iCls).scoreNorm;
      if isempty(scoreNorm) || isnan(scoreNorm)
        scoreNorm = 0;
      end
      confThresholds = obj.confThresholds(iCls,:);
    
      NBEH = 2;
      confMat = zeros(2*NBEH,3);
      for ndx = 1:2*NBEH
        if mod(ndx,2)
          curIdx = modLabels==ndx;
        else
          % either ndx or ndx-1; or lblN for both important and
          % nonimportant
          curIdx = modLabels>(ndx-1.5) & modLabels<(ndx+0.5);
        end
        confMat(ndx,1) = nnz( scores(curIdx)>= (confThresholds(1)*scoreNorm)); 
        confMat(ndx,2) = nnz(-scores(curIdx)<  (confThresholds(2)*scoreNorm) & ...
                              scores(curIdx)<  (confThresholds(1)*scoreNorm) );
        confMat(ndx,3) = nnz(-scores(curIdx)>= (confThresholds(2)*scoreNorm));
      end
      
      errorRates = struct();
      errorRates.numbers = confMat;
      errorRates.frac = errorRates.numbers./repmat( sum(errorRates.numbers,2),[1 3]);
    end

  end
  
  %% Show Similar Frames
  
  methods
    
    % ---------------------------------------------------------------------
    function DoFastBagging(obj)
      
      % MERGEST UPDATED
      
      % AL: Looks a lot like DoBagging
      
      assert(obj.nclassifiers==1,'Unsupported for multiple classifiers.');
      ICLS = 1;
      
      obj.StoreLabelsForCurrentAnimal();
      [success,msg] = obj.PreLoadPeriLabelWindowData();
      if ~success, 
        warning(msg);
        return;
      end

      islabeled = obj.windowdata(ICLS).labelidx_new ~= 0;

      if ~any(islabeled),                              return; end
      if ~strcmp(obj.classifiertype{ICLS},'boosting'); return; end
      if isempty(obj.classifier{ICLS}), obj.Train;             end

      if isempty(obj.windowdata(ICLS).binVals),
        obj.windowdata(ICLS).binVals = findThresholds(...
          obj.windowdata(ICLS).X(islabeled,:),...
          obj.classifier_params{ICLS},'deterministic',obj.deterministic);
        obj.windowdata(ICLS).bins = findThresholdBins(obj.windowdata(ICLS).X(islabeled,:),...
          obj.windowdata(ICLS).binVals);
      end
      bins = obj.windowdata(ICLS).bins;

      bmodel = fastBag(obj.windowdata(ICLS).X(islabeled,:),...
        obj.windowdata(ICLS).labelidx_new(islabeled),...
        obj.windowdata(ICLS).binVals,bins,obj.classifier_params{ICLS},obj);
      
      obj.bagModels = bmodel;
%       obj.bagModels = obj.classifier;

      obj.SetStatus('Computing parameters for fast distance computation..');
      % Find the parameters for fast prediction.
      feature_names = obj.windowdata(ICLS).featurenames;
      
      % which features are actually used
      dims = [obj.bagModels(:).dim];


      wfidxcurr = unique(dims,'stable');
      wfs = feature_names(wfidxcurr);
      feature_names = feature_names(dims);
            
      wf2pff = cellfun(@(x)x{1},wfs,'UniformOutput',false);
      [pffs,~,wf2pffidx] = unique(wf2pff);
      
      windowfeaturescellparams = struct;
      for pfi = 1:numel(pffs),
        pf = pffs{pfi};
        wfidx_cur = wf2pffidx==pfi;
        windowfeaturescellparams.(pf) = WindowFeatureName2Params(wfs(wfidx_cur));
      end
      
      classifiers_indexed = obj.bagModels;
      for j = 1:numel(classifiers_indexed),
        classifiers_indexed(j).dim = j;
      end
      
      obj.fastPredictBag.classifier = classifiers_indexed;
      obj.fastPredictBag.windowfeaturescellparams = windowfeaturescellparams;
      obj.fastPredictBag.wfs = feature_names;
      obj.fastPredictBag.pffs = pffs;
      obj.fastPredictBag.ts = obj.classifierTS(ICLS);
      obj.fastPredictBag.tempname = tempname;
 
      features_names_sel = {};
      for ndx = 1:numel(pffs)
        [~,curf] = ComputeWindowFeatures([0,0],...
          windowfeaturescellparams.(pffs{ndx}){:});
        feature_names_curr = cellfun(@(x) [{pffs{ndx}},x],curf,'UniformOutput',false);  %#ok
        features_names_sel = [features_names_sel,feature_names_curr]; %#ok<AGROW>
      end      
      
      ttt = tic;
      wfidx = nan(1,numel(feature_names));
      matched = false(1,numel(dims));
      for j = 1:numel(feature_names),
        if matched(j), continue, end
        
        idxcurr = find(WindowFeatureNameCompare(feature_names{j},features_names_sel));
        if numel(idxcurr) ~= 1,
          error('Error matching wfs for classifier with window features computed');
        end
        curidx = dims==dims(j);
        matched(curidx) = true;
        wfidx(curidx) = idxcurr;
        if (mod(j,100)==0)
          telapsed = toc(ttt);
          obj.SetStatus('Indexing the feature names ... %d%% done: Time Remaining:%ds',...
            round(nnz(matched)/numel(matched)*100),round(telapsed/nnz(matched)*(numel(matched)-nnz(matched))));
        end
      end
  
%       ttt = tic;
%       wfidx = nan(1,numel(feature_names));
%       for j = 1:numel(feature_names),
%         
%         idxcurr = find(WindowFeatureNameCompare(feature_names{j},features_names_sel));
%         if numel(idxcurr) ~= 1,
%           error('Error matching wfs for classifier with window features computed');
%         end
%         wfidx(j) = idxcurr;
%         if(mod(j,100)==0)
%           telapsed = toc(ttt);
%           obj.SetStatus('Indexing the feature names ... %d%% done: Time Remaining:%ds',...
%             round(j/numel(feature_names)*100),round(telapsed/j*(numel(feature_names)-j)));
%         end
%       end

      obj.fastPredictBag.wfidx = wfidx;
      
      obj.fastPredictBag.dist = cell(1,obj.nexps);
      obj.ClearStatus();      
    end
    

    % ---------------------------------------------------------------------
    function [success,msg] = SetCurrentFlyForBag(obj,exp,fly,t)
      % Initializes/sets .fastPredictBag.
      
      % MERGEST UPDATED
      
      assert(obj.nclassifiers==1,'Unsupported for multiple classifiers.');
      ICLS = 1;
      
      obj.fastPredictBag.curexp = exp;
      obj.fastPredictBag.fly = fly;
      obj.fastPredictBag.t = t;
      
      [success,msg,t0,t1,X] = obj.ComputeWindowDataChunk(exp,fly,t,'center',true,ICLS);
      
      curX = X((t0:t1)==t,:);
      curF = zeros(1,numel(obj.bagModels));
      for ndx = 1:numel(obj.bagModels);
        curWk = obj.bagModels(ndx);
        dd = curX(:,curWk.dim)*curWk.dir;
        tt = curWk.tr*curWk.dir;
        curF(ndx) = sign( (dd>tt) - 0.5) * curWk.alpha;
      end
      obj.fastPredictBag.curF = curF;
      obj.fastPredictBag.dist = {};
      obj.fastPredictBag.trainDist = {};
    end
    
    function [success,msg] = UnsetCurrentFlyForBag(obj)
      % Clears .fastPredictBag.
      
      %MERGEST OK
      
      success = true;
      msg = '';
      obj.fastPredictBag.curexp = [];
      obj.fastPredictBag.fly = [];
      obj.fastPredictBag.t = [];
      obj.fastPredictBag.dist = {};
      obj.fastPredictBag.trainDist = {};
      obj.fastPredictBag.curF = [];
    end

    
    % ---------------------------------------------------------------------
    function dist = GetDistance(obj,expi,flies)
      % MERGEST UPDATED
      
      assert(obj.nclassifiers==1,'Unsupported for multiple classifiers.');
      ICLS = 1;
      
      if obj.HasDistance(expi,flies)
        obj.ComputeBagDistanceTraining();
        distNorm = median(obj.fastPredictBag.trainDist);
        dist = obj.fastPredictBag.dist{expi}{flies};
        dist = dist/distNorm;
        dist(dist>1) = 1;
      elseif ~isempty(obj.bagModels) && ~isempty(obj.fastPredictBag.curexp),
        obj.ComputeBagDistanceTraining();
        distNorm = median(obj.fastPredictBag.trainDist);
        T0 = obj.firstframes_per_exp{expi}(flies);
        T1 = obj.endframes_per_exp{expi}(flies);
        dist = nan(1,T1-T0+1);
        idx = obj.FlyNdx(expi,flies,ICLS);
        t = obj.windowdata(ICLS).t(idx);
        dist(t-T0+1) = obj.fastPredictBag.trainDist(idx);
        dist = dist/distNorm;
        dist(dist>1) = 1;
      else
        T0 = obj.firstframes_per_exp{expi}(flies);
        T1 = obj.endframes_per_exp{expi}(flies);
        dist = nan(1, T1- T0 +1);
      end      
    end
    
    
    % ---------------------------------------------------------------------
    function [nextT, distT] = NextClosestBagFly(obj,dir,curt,expi,flies,curV,ignore,jumpList,jumpRestrict)
      
      % MERGEST UPDATED
      
      assert(obj.nclassifiers==1,'Unsupported for multiple classifiers.');
      ICLS = 1;

      nextT = []; distT = [];
      if  isempty(obj.fastPredictBag.dist) || ...
          isempty(obj.fastPredictBag.dist{expi}) || ...
         numel(obj.fastPredictBag.dist{expi})<flies || ...
         isempty(obj.fastPredictBag.dist{expi}{flies})
        [success,msg] = obj.ComputeBagDistFly(expi,flies);
        if ~success,
          uiwait(warndlg(msg));
          return;
        end
      end
      dist = obj.fastPredictBag.dist{expi}{flies};
      
      T0 = obj.firstframes_per_exp{expi}(flies);
      T1 = obj.endframes_per_exp{expi}(flies);
      if isempty(curV)
        curV = dist(curt-T0+1);
      end
      
      for ndx = 1:numel(jumpList)
        if jumpList(ndx).exp ~= expi || jumpList(ndx).fly ~= flies
          continue;
        end
        idx = jumpList(ndx).t - T0 + 1 + (-ignore:ignore);
        idx(idx<1) = [];
        idx(idx>numel(dist)) = [];
        dist(idx) = inf;
      end
      
      if strcmp(jumpRestrict,'behavior')
        obj.PredictFast(expi,flies,T0,T1,ICLS);
        obj.ApplyPostprocessing(expi,flies);
        dist(obj.predictdata{expi}{flies}(ICLS).cur < 0) = inf;
      elseif strcmp(jumpRestrict,'none')
        obj.PredictFast(expi,flies,T0,T1,ICLS);
        obj.ApplyPostprocessing(expi,flies);
        dist(obj.predictdata{expi}{flies}(ICLS).cur > 0) = inf;
      end
      
      [nextT,distT] = obj.FindNextClosest(dist,curV,dir);
      nextT = nextT+T0-1;
    end
    
    
    % ---------------------------------------------------------------------
    function [bestnextT,bestdistT,bestfly] = NextClosestBagExp(obj,dir,curt,expi,flies,ignore,jumpList,jumpRestrict)
      % MERGEST OK
      
      bestnextT = []; bestfly = [];
      switch dir
        case 'next'
          bestdistT = inf;
        case 'prev' 
          bestdistT = -inf;
        otherwise
          uiwait(warndlg('Undefined direction'));
          bestdistT = [];
          return;
      end
      
      if isempty(obj.fastPredictBag.dist) || ...
          isempty(obj.fastPredictBag.dist{expi}) || ...
          numel(obj.fastPredictBag.dist{expi})<flies || ...
          isempty(obj.fastPredictBag.dist{expi}{flies})
        [success,msg] = obj.ComputeBagDistFly(expi,flies);
        if ~success,
          uiwait(warndlg(msg));
          return;
        end
      end

      T0 = obj.firstframes_per_exp{expi}(flies);
      curV = obj.fastPredictBag.dist{expi}{flies}(curt-T0+1);
      
      for fly = 1:obj.nflies_per_exp(expi)
        [nextT,distT] = obj.NextClosestBagFly(dir,curt,expi,fly,curV,ignore,jumpList,jumpRestrict);
        if strcmp(dir,'next')
          if distT < bestdistT
            bestdistT = distT;
            bestnextT = nextT;
            bestfly = fly;
          end
        else
          if distT > bestdistT
            bestdistT = distT;
            bestnextT = nextT;
            bestfly = fly;
          end          
        end
      end      
    end
    
    
    % ---------------------------------------------------------------------
    function [nextT, distT, fly, exp ] = ...
        NextClosestBagTraining(obj,jtype,curT,curExp,curFly,ignore,jumpList,jumpRestrict)
      
      %MERGEST UPDATED
      
      assert(obj.nclassifiers==1,'Unsupported for multiple classifiers.');
      ICLS = 1;
      
      nextT = []; distT = []; fly = []; exp = [];
      
      obj.ComputeBagDistanceTraining();
      
      dist = obj.fastPredictBag.trainDist;
      
      curndx = find(obj.FlyNdx(curExp,curFly,ICLS) & obj.windowdata(ICLS).t == curT);
      if isempty(curndx)
        curV = 0;
      else
        curV = dist(curndx);
      end
      
      % Ignore close by bouts to previously seen examples
      for ndx = 1:numel(jumpList)
        idx = obj.FlyNdx(jumpList(ndx).exp,jumpList(ndx).fly,ICLS) & ...
              (abs(obj.windowdata(ICLS).t - jumpList(ndx).t)<=ignore );
        dist(idx) = inf;
      end
      
      labels = obj.windowdata(ICLS).labelidx_new;
      if strcmp(jumpRestrict,'behavior')
        dist(labels ~= 1) = inf;
      elseif strcmp(jumpRestrict,'none')
        dist(labels ~= 2) = inf;
      end
      
      nextNdx = [];
      switch jtype
        case 'next'
          
          if max(dist) > curV,
            tt = dist-curV;
            tt(tt<=0) = inf;
            [~,nextNdx] = min(tt);
            distT = dist(nextNdx);
          end
          
        case 'prev'

          if min(dist) < curV,
            tt = dist-curV;
            tt(tt>=0) = -inf;
            [~,nextNdx] = max(tt);
            distT = dist(nextNdx);
          end
        
      end
      
      if ~isempty(nextNdx)
        fly = obj.windowdata(ICLS).flies(nextNdx);
        exp = obj.windowdata(ICLS).exp(nextNdx);
        nextT = obj.windowdata(ICLS).t(nextNdx);
      end
    end
    
    
    % ---------------------------------------------------------------------
    function DoBagging(obj)
      % MERGEST UPDATED
      
      assert(obj.nclassifiers==1,'Unsupported for multiple classifiers.');
      ICLS = 1;
      
      obj.StoreLabelsForCurrentAnimal();
      [success,msg] = obj.PreLoadPeriLabelWindowData();
      if ~success,
        warning(msg);
        return;
      end
      
      islabeled = obj.windowdata(ICLS).labelidx_new ~= 0;
      
      if ~any(islabeled),                              return; end
      if ~strcmp(obj.classifiertype{ICLS},'boosting'); return; end
      if isempty(obj.classifier{ICLS}), obj.Train;             end
      
      bouts = obj.getLabeledBouts(ICLS);
      
      obj.SetStatus('Bagging the classifier with %d examples...',nnz(islabeled));
      
      obj.windowdata(ICLS).binVals = ...
        findThresholds(obj.windowdata(ICLS).X(islabeled,:),obj.classifier_params{ICLS},'deterministic',obj.deterministic);
      obj.windowdata(ICLS).bins = findThresholdBins(obj.windowdata(ICLS).X,obj.windowdata(ICLS).binVals);
      
      [obj.bagModels,obj.distMat] = ...
        doBaggingBouts(obj.windowdata(ICLS).X, ...
        obj.windowdata(ICLS).labelidx_new,obj, ...
        obj.windowdata(ICLS).binVals,...
        obj.classifier_params{ICLS},bouts);
      
      obj.windowdata(ICLS).distNdx.exp = obj.windowdata(ICLS).exp(islabeled);
      obj.windowdata(ICLS).distNdx.flies = obj.windowdata(ICLS).flies(islabeled);
      obj.windowdata(ICLS).distNdx.t = obj.windowdata(ICLS).t(islabeled);
      obj.windowdata(ICLS).distNdx.labels = obj.windowdata(ICLS).labelidx_new(islabeled);
      
      obj.ClearStatus();
    end
    
    
    % ---------------------------------------------------------------------
    function SimilarFrames(obj,curTime,JLabelHandles)
      
      % MERGEST UPDATED
      
      assert(obj.nclassifiers==1,'Unsupported for multiple classifiers.');
      ICLS = 1;
      
      if isempty(obj.frameFig) || ~ishandle(obj.frameFig),
        obj.InitSimilarFrames(JLabelHandles),
      end
      
      distNdx = find( (obj.windowdata(ICLS).distNdx.exp == obj.expi) & ...
        (obj.windowdata(ICLS).distNdx.flies == obj.flies) & ...
        (obj.windowdata(ICLS).distNdx.t == curTime),1);
      
      windowNdx = find( (obj.windowdata(ICLS).exp == obj.expi) & ...
        (obj.windowdata(ICLS).flies == obj.flies) & ...
        (obj.windowdata(ICLS).t == curTime),1);
      
      if isempty(distNdx) % The example was not part of the training data.
        outOfTraining = 1;
        [~,~,t0,~,curX] = obj.ComputeWindowDataChunk(obj.expi,obj.flies,curTime,[],[],ICLS);

        curX = curX(curTime-t0+1,:);
        curD = zeros(1,length(obj.bagModels)*length(obj.bagModels{1}));
        count = 1;
        for bagNo = 1:length(obj.bagModels)
          curModel = obj.bagModels{bagNo};
          for j = 1:length(curModel)
            curWk = curModel(j);
            dd = curX(curWk.dim)*curWk.dir;
            tt = curWk.tr*curWk.dir;
            curD(count) = (dd>tt)*curWk.alpha;
            count = count+1;
          end
        end
      else
        outOfTraining = 0;
        curD = obj.distMat(distNdx,:);
      end

      % Compute the distance 
      diffMat = zeros(size(obj.distMat));
      for ndx = 1:size(diffMat,2);
        diffMat(:,ndx) = abs(obj.distMat(:,ndx)-curD(ndx));
      end
      dist2train = nanmean(diffMat,2)*200;
      [rr,rrNdx] = sort(dist2train,'ascend');
      
      if ~outOfTraining
        rr = rr(2:end); %#ok
        curEx = rrNdx(1); rrNdx = rrNdx(2:end); %#ok
      else
        curEx = [];  %#ok
      end
      
      % Find 5 closest pos and neg examples.
      % This looks complicated then it should be.
      % DEBUG: find values of actual labels 
     
      trainLabels = obj.windowdata(ICLS).distNdx.labels;
      allPos = rrNdx(trainLabels(rrNdx)>1.5);
      allNeg = rrNdx(trainLabels(rrNdx)<1.5);      
      
      curP = zeros(1,5);
      curN = zeros(1,5);
      count = 0;
      for ex = allPos'
        if count>4; break; end;
        isClose = 0;
        if ~outOfTraining && ...
          obj.windowdata(ICLS).exp(windowNdx) == obj.windowdata(ICLS).distNdx.exp(ex) && ...
           obj.windowdata(ICLS).flies(windowNdx) == obj.windowdata(ICLS).distNdx.flies(ex) && ...
           abs( obj.windowdata(ICLS).t(windowNdx) - obj.windowdata(ICLS).distNdx.t(ex) )<5,
           continue; 
        end
        
        for used = curP(1:count)
          if obj.windowdata(ICLS).distNdx.exp(used) == obj.windowdata(ICLS).distNdx.exp(ex) && ...
             obj.windowdata(ICLS).distNdx.flies(used) == obj.windowdata(ICLS).distNdx.flies(ex) && ...
             abs( obj.windowdata(ICLS).distNdx.t(used) - obj.windowdata(ICLS).distNdx.t(ex) )<5,
             isClose = 1;
             break;
          end
        end
        
        if isClose; continue; end;
        count = count+1;
        curP(count) = ex;
      end
      
      count = 0;
      for ex = allNeg'
        if count>4; break; end;
        isClose = 0;
        if ~outOfTraining && ...
          obj.windowdata(ICLS).exp(windowNdx) == obj.windowdata(ICLS).distNdx.exp(ex) && ...
           obj.windowdata(ICLS).flies(windowNdx) == obj.windowdata(ICLS).distNdx.flies(ex) && ...
           abs(obj.windowdata(ICLS).t(windowNdx) - obj.windowdata(ICLS).distNdx.t(ex))<5,
           continue; 
        end
        
        for used = curN(1:count)
          if obj.windowdata(ICLS).distNdx.exp(used) == obj.windowdata(ICLS).distNdx.exp(ex) && ...
             obj.windowdata(ICLS).distNdx.flies(used) == obj.windowdata(ICLS).distNdx.flies(ex) && ...
             abs(obj.windowdata(ICLS).distNdx.t(used) - obj.windowdata(ICLS).distNdx.t(ex))<5,
             isClose = 1; 
             break; 
          end
        end
        
        if isClose; continue; end;
        count = count+1;
        curN(count) = ex;
      end
      
      varForSSF.curFrame.expNum = obj.expi;
      varForSSF.curFrame.flyNum = obj.flies;
      varForSSF.curFrame.curTime = curTime;
      
      for k = 1:4
        varForSSF.posFrames(k).expNum = obj.windowdata(ICLS).distNdx.exp(curP(k));
        varForSSF.posFrames(k).flyNum = obj.windowdata(ICLS).distNdx.flies(curP(k));
        varForSSF.posFrames(k).curTime = obj.windowdata(ICLS).distNdx.t(curP(k));
        varForSSF.negFrames(k).expNum = obj.windowdata(ICLS).distNdx.exp(curN(k));
        varForSSF.negFrames(k).flyNum = obj.windowdata(ICLS).distNdx.flies(curN(k));
        varForSSF.negFrames(k).curTime = obj.windowdata(ICLS).distNdx.t(curN(k));
      end
      showSimilarFrames('setFrames',obj.frameFig,varForSSF);
    end  % method
    
  end
  
  methods (Access=private)
    
    % ---------------------------------------------------------------------
    function [success,msg,dist] = ComputeBagFeatures(obj,curexp,curfly,curF)
      % Use the fast feature computation to find the bag features.
      
      % MERGESTOK
      
      success = true; msg = '';
      featureFileName = sprintf('%s_%s_%d',obj.fastPredictBag.tempname,obj.expnames{curexp},curfly);
      if exist(featureFileName,'file'),
        load(featureFileName,'dX'),
        dist = sum(abs(dX - repmat(curF,[size(dX,1) 1]))); %#ok<NODEF>
        return;
      end
      
      obj.SetStatus(sprintf('Computing distance to fly %d in exp:%s',curfly,obj.expnames{curexp}));
      perframeInMemory = ~isempty(obj.flies) && obj.IsCurFly(curexp,curfly);
      perframefile = obj.GetPerframeFiles(curexp);
      
      T0 = obj.GetTrxFirstFrame(curexp,curfly);
      T1 = obj.GetTrxEndFrame(curexp,curfly);
      
      perframedata_cur = obj.perframedata;
      windowfeaturescellparams = obj.fastPredictBag.windowfeaturescellparams;
      pffs = obj.fastPredictBag.pffs;
      allperframefns = obj.allperframefns;
      
      X_all = [];
      for t0 = T0:(2*obj.predictwindowdatachunk_radius):T1
        t1 = min(T1,t0+2*obj.predictwindowdatachunk_radius-1)-T0+1;
        
        % for the parfor loop.
        x_curr_all = cell(1,numel(pffs));
        X = []; fnames = {};
        parfor j = 1:numel(pffs),
          
          fn = pffs{j};
          
          ndx = find(strcmp(fn,allperframefns));
          if perframeInMemory,
            perframedata = perframedata_cur{ndx};  %#ok
          else
            perframedata = load(perframefile{ndx});  %#ok
            perframedata = perframedata.data{curfly(1)};  %#ok
          end
          
          t11 = min(t1,numel(perframedata));
          [x_curr,curf] = ...
            ComputeWindowFeatures(perframedata,...
            windowfeaturescellparams.(fn){:},'t0',t0-T0+1,'t1',t11);  %#ok
          fnames{j} = curf;
          if t11 < t1,
            x_curr(:,end+1:end+t1-t11) = nan;
          end
          
          x_curr_all{j} = single(x_curr);
        end
        
        allFeatures = {};
        for j = 1:numel(pffs),
          fn = pffs{j};
          x_curr = x_curr_all{j};
          % add the window data for this per-frame feature to X
          nold = size(X,1);
          nnew = size(x_curr,2);
          if nold > nnew,
            warning(['Number of examples for per-frame feature %s does not '...
              'match number of examples for previous features'],fn);
            x_curr(:,end+1:end+nold-nnew) = nan;
          elseif nnew > nold && ~isempty(X),
            warning(['Number of examples for per-frame feature %s does not '...
              'match number of examples for previous features'],fn);
            X(end+1:end+nnew-nold,:) = nan;
          end
          X = [X,x_curr']; %#ok<AGROW>
          jj = fnames{j};
          feature_names_curr = cellfun(@(x) [{pffs{j}},x],jj,'UniformOutput',false);  %#ok
          
          allFeatures = [allFeatures,feature_names_curr];  %#ok
        end
        
        X_all = [X_all;X];  %#ok
      end
      
      X_all = X_all(:,obj.fastPredictBag.wfidx);
      
      dX = zeros(size(X_all,1),numel(obj.fastPredictBag.classifier));
      for ndx = 1:numel(obj.fastPredictBag.classifier)
        curWk = obj.fastPredictBag.classifier(ndx);
        dd = X_all(:,curWk.dim)*curWk.dir;
        tt = curWk.tr*curWk.dir;
        dX(:,ndx) = sign( (dd>tt)-0.5 )*curWk.alpha;
        
      end
      save(featureFileName,'dX');
      dist = sum(abs(dX - repmat(curF,[size(dX,1) 1])),2);
      obj.ClearStatus();
    end
    
    
    % ---------------------------------------------------------------------
    function [success,msg] = ComputeBagDistFly(obj,expi,fly)
      % MERGESTOK
      [success,msg,dist] = obj.ComputeBagFeatures(expi,fly,obj.fastPredictBag.curF);
      obj.fastPredictBag.dist{expi}{fly} = dist;
    end
    
    
    % ---------------------------------------------------------------------
    function [success,msg] = ComputeBagDistExp(obj,curexp,curfly,curt,expi)  %#ok
      % AL: no callsites?
      for ndx = 1:numel(obj.nflies_per_exp(expi))
        [success,msg,dist] = obj.ComputeBagFeatures(expi,ndx,obj.fastPredict.curF);
        obj.fastPredictBag.dist{expi}{ndx} = dist;
      end
      
    end
    
    
    % ---------------------------------------------------------------------
    function [nextT,distT] = FindNextClosest(obj,dist,curV,dir) %#ok
      % MERGEST OK
      
      nextT = []; distT = [];
      switch dir
        case 'next'
          
          if max(dist) > curV,
            tt = dist-curV;
            tt(tt<=0) = inf;
            [~,nextT] = min(tt);
            distT = dist(nextT);
          end
          
        case 'prev'
          
          if min(dist) < curV,
            tt = dist-curV;
            tt(tt>=0) = -inf;
            [~,nextT] = max(tt);
            distT = dist(nextT);
          end
          
      end
      
    end
    
    
    % ---------------------------------------------------------------------
    function has = HasDistance(obj,expi,flies)
      % MERGEST OK
      has = ~(isempty(obj.bagModels) || ...
        isempty(obj.fastPredictBag.dist) || ...
        numel(obj.fastPredictBag.dist)<expi || ...
        isempty(obj.fastPredictBag.dist{expi}) || ...
        numel(obj.fastPredictBag.dist{expi})<flies || ...
        isempty(obj.fastPredictBag.dist{expi}{flies}));
    end
    
    
    % ---------------------------------------------------------------------
    function ComputeBagDistanceTraining(obj)
      % MERGEST UPDATED
      
      assert(obj.nclassifiers==1,'Unsupported for multiple classifiers.');
      ICLS = 1;
      
      if isempty(obj.fastPredictBag.trainDist)
        dX = zeros(size(obj.windowdata(ICLS).X,1),numel(obj.bagModels));
        for ndx = 1:numel(obj.bagModels)
          curWk = obj.bagModels(ndx);
          dd = obj.windowdata(ICLS).X(:,curWk.dim)*curWk.dir;
          tt = curWk.tr*curWk.dir;
          dX(:,ndx) = sign( (dd>tt)-0.5 )*curWk.alpha;
        end
        
        dist = sum(abs(dX - repmat(obj.fastPredictBag.curF,[size(dX,1) 1])),2);
        
        obj.fastPredictBag.trainDist = dist;
      end
    end
    
    
    % ---------------------------------------------------------------------
    function InitSimilarFrames(obj,HJLabel)
      obj.frameFig = showSimilarFrames;
      showSimilarFrames('SetJLabelData',obj.frameFig,obj,HJLabel);
      showSimilarFrames('CacheTracksLabeled',obj.frameFig);
      showSimilarFrames('add_prep_list', obj.frameFig);
    end
    
  end  
  
  %% Ground truthing
  
  methods 

    % ---------------------------------------------------------------------
    function [success,msg] = setGTSuggestionMode(obj,modeString,varargin)
      % Sets the ground-truth suggestion mode. This determines how things
      % returned by GetGTSuggestionIdx() are calculated.
      
      %MERGEST OK
      
      switch modeString,
        case 'Random'
          [success,msg] = obj.SuggestRandomGT(varargin{:});
          obj.needsave = true;
        case 'Balanced'
          [success,msg] = obj.SuggestBalancedGT(varargin{:});
          obj.needsave = true;
        case 'Imported'
          [success,msg] = obj.SuggestLoadedGT(varargin{:});
          obj.needsave = true;
        case 'Threshold'
          [success,msg] = obj.SuggestThresholdGT(varargin{:});
          obj.needsave = true;
        otherwise
          error('JLabelData:noSuchGTSuggestionMode', ...
                'Internal error: No such ground-truth suggestion mode');
      end
    end


    % ---------------------------------------------------------------------
    function avgBoutLen = GetAverageLoadedPredictionBoutLength(obj)
      %MERGEST OK
      assert(obj.nclassifiers==1,'Only supported for single classifiers.');
      
      if ~obj.HasLoadedScores(1)
        error('JLabelData:noLoadedScores','No scores have been loaded.');
      end
      
      blen = [];
      for endx = 1:obj.nexps
        for flies = 1:obj.nflies_per_exp(endx)
          curidx = obj.predictdata{endx}{flies}.loaded_valid;
          posts = obj.predictdata{endx}{flies}.loaded(curidx)>0;
          labeled = bwlabel(posts);
          aa = regionprops(labeled,'Area'); %#ok
          blen = [blen [aa.Area]]; %#ok
        end
      end
      avgBoutLen = mean(blen);
      
    end
    
    
    % ---------------------------------------------------------------------
    function SaveSuggestionGT(obj,filename)
      
      %MERGEST OK
      
      assert(obj.nclassifiers==1,'Only supported for single classifiers.');

      fid = fopen(filename,'w');
      switch obj.GTSuggestionMode
        case 'Random'
          for expi = 1:obj.nexps
            for fly = 1:obj.nflies_per_exp(expi)
              start = obj.randomGTSuggestions{expi}(fly).start;
              last = obj.randomGTSuggestions{expi}(fly).end;
              fprintf(fid,'exp:%s,fly:%d,start:%d,end:%d\n',obj.expnames{expi},fly,start,last);
            end
          end
        case 'Threshold'
          for expi = 1:obj.nexps
            if obj.predictdata{expi}{1}.loaded_valid(1)
              for fly = 1:obj.nflies_per_exp(expi)
                T1 = obj.GetTrxEndFrame(expi,fly);
                suggestedidx = zeros(1,T1);
                suggestedidx(obj.predictdata{expi}{fly}.t) = ...
                  obj.NormalizeScores(obj.predictdata{expi}{fly}.loaded) > ...
                  -obj.thresholdGTSuggestions;
                [t0s,t1s] = get_interval_ends(suggestedidx);
                for ndx = 1:numel(t0s)
                  if t1s(ndx)< T0, continue; end
                  fprintf(fid,'exp:%s,fly:%d,start:%d,end:%d\n',obj.expnames{expi}.fly,t0s(ndx),t1s(ndx));
                end
              end
            end
          end
        case 'Balanced'
          for ndx = 1:numel(obj.balancedGTSuggestions)
%             if obj.balancedGTSuggestions(ndx).exp ~= expi
%               continue;
%             end
            start = obj.balancedGTSuggestions(ndx).start;
            last = obj.balancedGTSuggestions(ndx).end;
            fprintf(fid,'exp:%s,fly:%d,start:%d,end:%d\n',...
              obj.expnames{obj.balancedGTSuggestions(ndx).exp},...
              obj.balancedGTSuggestions(ndx).flies,start,last);
          end
        case 'Imported'
          for expi = 1:obj.nexps
            for fly = 1:obj.nflies_per_exp(expi)
              start = obj.loadedGTSuggestions{expi}(fly).start;
              last = obj.loadedGTSuggestions{expi}(fly).end;
              if numel(start) == 1 && last < start,
                continue;
              end
              for ndx = 1:numel(start)
                fprintf(fid,'exp:%s,fly:%d,start:%d,end:%d\n',obj.expnames{expi},fly,start(ndx),last(ndx));
              end
            end
          end
      end
      fclose(fid);

    end

    
    % ---------------------------------------------------------------------
    function suggestedidx = GetGTSuggestionIdx(obj,expi,flies,T0,T1)

      % MERGEST OK
      
      assert(obj.nclassifiers==1,'Only supported for single classifiers.');

      % Get the indices of gt suggestion  
      if nargin<4
        T0 = obj.GetTrxFirstFrame(expi,flies);
        T1 = obj.GetTrxEndFrame(expi,flies);
      end
      n = T1-T0+1;
      off = 1 - T0;
      
      if isempty(obj.GTSuggestionMode)
        suggestedidx = false(1,n);
        return;
      end
      
      suggestedidx = false(1,n);
      
      switch obj.GTSuggestionMode,
        case 'Random'
          if numel(obj.randomGTSuggestions)<expi || isempty(obj.randomGTSuggestions{expi}),
            suggestedidx = false(1,n);
            return;
          end
          start = obj.randomGTSuggestions{expi}(flies).start;
          last = obj.randomGTSuggestions{expi}(flies).end;
          range = start+off:last+off;
          selIdx = range(range>0);
          suggestedidx(selIdx) = true;
        
        case 'Imported'
          if numel(obj.loadedGTSuggestions)<expi || isempty(obj.loadedGTSuggestions{expi}),
            suggestedidx = false(1,n);
            return;
          end
          suggestedidx = false(1,n);
          for ndx = 1:numel(obj.loadedGTSuggestions{expi}(flies).start)
            start = obj.loadedGTSuggestions{expi}(flies).start(ndx);
            last = obj.loadedGTSuggestions{expi}(flies).end(ndx);
            range = start+off:last+off;
            selIdx = range(range>0);
            suggestedidx(selIdx) = true;
          end
          suggestedidx(n+1:end) = [];

        case 'Threshold'
          if (obj.predictdata{expi}{flies}.loaded_valid(1))
            suggestedidx = ...
              obj.NormalizeScores(obj.predictdata{expi}{flies}.loaded) > ...
              -obj.thresholdGTSuggestions;
          elseif any(obj.predictdata{expi}{flies}.cur_valid)
            idxcurr = obj.predictdata{expi}{flies}.cur_valid;
            suggestedidx = false(size(idxcurr));
            suggestedidx(idxcurr) = obj.NormalizeScores(obj.predictdata{expi}{flies}.cur(idxcurr)) > ...
              -obj.thresholdGTSuggestions;
          end
          
        case 'Balanced'
          suggestedidx = false(1,n);
          for ndx = 1:numel(obj.balancedGTSuggestions)
            if obj.balancedGTSuggestions(ndx).exp ~= expi ||...
              obj.balancedGTSuggestions(ndx).flies ~= flies
              continue;
            end
            start = obj.balancedGTSuggestions(ndx).start;
            last = obj.balancedGTSuggestions(ndx).end;
            if start>T1 || last <T0, continue ;end
            start = max(start,T0); last = min(last,T1);
            range = start+off:last+off;
            selIdx = range(range>0);
            suggestedidx(selIdx) = true;
          end
      end
      
    end

    
    % ---------------------------------------------------------------------
    function crossError = GetGTPerformance(obj)
      % Computes classifier performance on the GT data.
      
      %MERGEST UPDATED
      
      assert(obj.nclassifiers==1,'Only supported for single classifiers.');
      assert(obj.gtMode,'JLabelData:wrongMode',...
              'Can only call GetGTPerformance() in ground-truthing mode.');

      obj.StoreLabelsAndPreLoadWindowData();
      
      hasloaded = obj.HasLoadedScores(1);
      if ~hasloaded
        for expi = 1:obj.nexps
          for i = 1:size(obj.labels(expi).flies,1)
            flies = obj.labels(expi).flies(i,:);
            labels_curr = obj.GetLabels(expi,flies);
%             ts = [];            
            for j = 1:numel(labels_curr.t0s)
%               ts = [ts,labels_curr.t0s(j):(labels_curr.t1s(j)-1)]; %#ok<AGROW>
              obj.PredictFast(expi,flies,labels_curr.t0s(j),labels_curr.t1s(j)-1,1);
            end
            
            obj.ApplyPostprocessing(expi,flies);
            
            % assumes that if have any loaded score for an experiment we
            % have scores for all the flies and for every frame.
%             [success1,msg] = obj.PreLoadWindowData(expi,flies,ts);
%             if ~success1,
%               warndlg(msg);
%               return;
%             end
          end
        end
%        obj.PredictLoaded();
      end
      
      gt_scores =[];
      gt_labels = [];
      
      for expi = 1:obj.nexps,
        for i = 1:size(obj.labels(expi).flies,1),
          
          flies = obj.labels(expi).flies(i,:);
          labels_curr = obj.GetLabels(expi,flies);
          
          % Find the important labels
          labels_imp = [];
          for j = 1:numel(labels_curr.imp_t0s),
            t0 = labels_curr.imp_t0s(j);
            t1 = labels_curr.imp_t1s(j);
            labels_imp = [labels_imp t0:t1-1];  %#ok
          end
          
          for j = 1:numel(labels_curr.t0s),
            t0 = labels_curr.t0s(j);
            t1 = labels_curr.t1s(j);
            
            curLabel = 2*repmat(find(strcmp(labels_curr.names{j},obj.labelnames)),1,t1-t0); 
            % single classifier: 2 for beh, 4 for no-beh
            curLabel(ismember(t0:t1-1,labels_imp)) = curLabel(ismember(t0:t1-1,labels_imp))-1;
            % now, 1/2/3/4 <-> behImp/behNotImp/noBehImp/noBehNotImp
            
            gt_labels = [gt_labels curLabel]; %#ok
            
            if hasloaded,
              idx = obj.predictdata{expi}{flies}.t(:) >=t0 & obj.predictdata{expi}{flies}.t(:) <t1;
              ts = obj.predictdata{expi}{flies}.t(idx);
              scores = obj.predictdata{expi}{flies}.loaded_pp(idx)-0.5;
              [check,ndxInLoaded] = ismember(t0:(t1-1),ts);
              if any(check==0), warndlg('Loaded scores are missing scores for some loaded frames'); end
              gt_scores = [gt_scores scores(ndxInLoaded)];  %#ok
            else
              idx = obj.predictdata{expi}{flies}.t(:) >=t0 & obj.predictdata{expi}{flies}.t(:) <t1;
              ts = obj.predictdata{expi}{flies}.t(idx);
              scores = obj.predictdata{expi}{flies}.cur_pp(idx)-0.5;
               [check,ndxInLoaded] = ismember(t0:(t1-1),ts);
              if any(check==0), warndlg('calculated scores are missing for some labeled frames'); end
              gt_scores = [gt_scores scores(ndxInLoaded)];  %#ok
            end
          end
          
        end
      end

      crossError = obj.createConfMat(1,gt_scores,gt_labels);
    end
    
    
    % ---------------------------------------------------------------------
    function gtMode = IsGTMode(obj)
      gtMode = obj.gtMode;
    end
    
  end
  
  methods (Access=private) 
    
    % ---------------------------------------------------------------------
    function [success,msg] = SuggestRandomGT(obj,perfly,perexp)
      % Set obj.randomGTSuggestions, obj.GTSuggestionMode
      %
      % This is currently public, but seems like maybe it should be private,
      % and get called when callers actually query
      % self.randomGTSuggestions, or something.  -- ALT, Apr 19, 2013
      
      % MERGEST OK
      
      success = false;
      msg = '';
      
      % Do nothing if we already have suggestions with the same settings for
      % all the experiments.
      if numel(obj.randomGTSuggestions)<numel(obj.nflies_per_exp)
        recompute = true;
      else
        recompute = false;
        for endx = 1:obj.nexps
          prevperexp = 0;
          
          for fndx = 1:obj.nflies_per_exp(endx)
            if isempty(obj.randomGTSuggestions{endx}(fndx).start);
              continue;
            end
            prevperfly = obj.randomGTSuggestions{endx}(fndx).end - ...
              obj.randomGTSuggestions{endx}(fndx).start+1;
            if (prevperfly ~= perfly)
              recompute = true;
              break; % AL: Ideally breaks out of both loops, not just inner?
            end
            prevperexp = prevperexp+1;
          end
          
          if prevperexp ~= perexp
            recompute = true;
          end
        end
      end
      if ~recompute
        obj.GTSuggestionMode = 'Random';
        return;
      end
      
      for endx = 1:obj.nexps
        % AL: appears to be the wrong initialization, don't we want
        % repmat(...,1,obj.nflies_per_exp(endx)?
        obj.randomGTSuggestions{endx} = repmat(struct('start',[],'end',[]),1,perexp);
        
        validflies = find( (obj.endframes_per_exp{endx} - ...
          obj.firstframes_per_exp{endx})>perfly );
        if numel(validflies)<perexp
          msg = sprintf('Experiment %s does not have %d flies with more than %d frames',...
            obj.expnames{endx},perexp,perfly);
          success = false;
          return;
        end
        permuteValid = validflies(randperm(numel(validflies)));
        randFlies = permuteValid(1:perexp);
        
        for fndx = 1:obj.nflies_per_exp(endx)
          if any(fndx == randFlies)
            first = obj.firstframes_per_exp{endx}(fndx);
            last = obj.endframes_per_exp{endx}(fndx);
            suggestStart = first + round( (last-first-perfly)*rand(1));
            obj.randomGTSuggestions{endx}(fndx).start = suggestStart;
            obj.randomGTSuggestions{endx}(fndx).end = suggestStart+perfly-1;
          else
            obj.randomGTSuggestions{endx}(fndx).start = [];
            obj.randomGTSuggestions{endx}(fndx).end = [];
          end
        end
      end
      success = true;
      obj.GTSuggestionMode = 'Random';
      obj.needsave = true;
    end
    
    
    % ---------------------------------------------------------------------
    function [success,msg] = SuggestBalancedGT(obj,intsize,numint)
      % Set obj.balancedGTSuggestions, obj.GTSuggestionMode
      %
      % Suggest frames such that half the suggested frames are predicted
      % positive and half are negative.
      % Excludes frames that have normal labels.
      
      % MERGEST OK
      
      assert(obj.nclassifiers==1,'Unsupported for multiple classifiers.');
      
      success = true; msg = '';
      
      if ~obj.HasLoadedScores(1)
        uiwait(warndlg('No scores have been loaded. Load precomputed scores to use this'));
      end
      
      jabparams = load(obj.everythingFileNameAbs,'-mat');
      assert(numel(jabparams.x.classifierStuff)==1,'Unsupported for multiple classifiers.');
      frames2exclude = cell(1,obj.nexps);
      for ndx = numel(jabparams.x.expDirNames)
        matchingGtexp = find(strcmp(jabparams.x.expDirNames{ndx},obj.expdirs));
        if isempty(matchingGtexp), continue; end
        frames2exclude{matchingGtexp} = jabparams.x.labels(ndx);
      end
      
      numpos = 0;
      numneg = 0;
      for expi = 1:obj.nexps,
        for flies = 1:obj.nflies_per_exp(expi)
          numpos = numpos + nnz(obj.predictdata{expi}{flies}.loaded>0);
          numneg = numneg + nnz(obj.predictdata{expi}{flies}.loaded<0);
        end
      end
      poswt = numneg/(numneg+numpos);
      negwt = numpos/(numneg+numpos);
      
      int = struct('exp',[],'flies',[],'tStart',[],'wt',[]);
      obj.balancedGTSuggestions = {};
      for endx = 1:obj.nexps
        obj.SetStatus('Couting predictions on experiment%d',endx);
        for flies = 1:obj.nflies_per_exp(endx)
          if ~obj.predictdata{endx}{flies}.loaded_valid(1),
            msg = sprintf('No Scores have been loaded for %s, cannot suggest intervals for ground truthing\n',...
              obj.expnames{endx});
            success = false;
            return;
          end
          
          curt = obj.predictdata{endx}{flies}.t;
          if numel(curt)<intsize; continue; end
          numT = numel(curt)-intsize+1;
          int.exp(1,end+1:end+numT) = endx;
          int.flies(1,end+1:end+numT) = flies;
          int.tStart(1,end+1:end+numT) = curt(1:end-intsize+1);
          curwt = (obj.predictdata{endx}{flies}.loaded<0)*negwt +(obj.predictdata{endx}{flies}.loaded>0)*poswt ;
          
          cumwt = cumsum(curwt);
          sumwt = cumwt(intsize+1:end)-cumwt(1:end-intsize);
          sumwt = [cumwt(intsize) sumwt]; %#ok<AGROW>
          
          if ~isempty(frames2exclude{endx}) && any(frames2exclude{endx}.flies == flies)
            
            fndx = find(frames2exclude{endx}.flies == flies);
            labeledF = false(size(curwt));
            for bndx = 1:numel(frames2exclude{endx}.t0s{fndx})
              tStart = frames2exclude{endx}.off(fndx);
              curt0 = max(1, frames2exclude{endx}.t0s{fndx}(bndx)-intsize + tStart );
              curt1 = min(numel(sumwt),frames2exclude{endx}.t1s{fndx}(bndx)+intsize +tStart);
              labeledF(curt0:curt1-1) = true;
            end
            
            sumwt(labeledF) = 0;
          end
          
          int.wt(1,end+1:end+numT) = sumwt;
          
        end
      end
      
      obj.balancedGTSuggestions = [];
      for ndx = 1:numint
        obj.SetStatus('Finding interval %d to label',ndx);
        
        % weight sampling was off by 1
        % fixed 20140331 by KB
        
        % old sampling
        %cumwt = cumsum(int.wt)/sum(int.wt);
        cumwt = cumsum([0,int.wt(1:end-1)])/sum(int.wt);
        intlocs = rand;
        locsSel = find(cumwt<=intlocs,1,'last');
        
        if isempty(locsSel), locsSel = numel(cumwt); end
        expi = int.exp(locsSel);
        flies = int.flies(locsSel);
        tStart = int.tStart(locsSel);
        obj.balancedGTSuggestions(ndx).start = tStart;
        obj.balancedGTSuggestions(ndx).end = tStart+intsize-1;
        obj.balancedGTSuggestions(ndx).exp = expi;
        obj.balancedGTSuggestions(ndx).flies = flies;
        
        % Removing intervals that overlap
        overlap = int.exp == int.exp(locsSel) & ...
          int.flies == int.flies(locsSel) & ...
          abs( int.tStart-int.tStart(locsSel))<=intsize;
        int.exp(overlap) = [];
        int.flies(overlap) = [];
        int.tStart(overlap) = [];
        int.wt(overlap) = [];
      end
      
      obj.GTSuggestionMode = 'Balanced';
      obj.needsave = true;
    end
    
    
    % ---------------------------------------------------------------------
    function [success,msg] = SuggestLoadedGT(obj,filename)
      %MERGEST OK
      
      success = false; %#ok<NASGU>
      msg = '';
      fid = fopen(filename);
      if fid < 0,
        msg = sprintf('Could not open file %s for reading',filename);
      end
      dat = cell(0,4);
      missingexps = {};
      while true,
        s = fgetl(fid);
        if ~ischar(s),
          break;
        end
        s = strtrim(s);
        m = regexp(s,'^exp:(.*),fly:(.*),start:(.*),end:(.*)','tokens','once');
        if isempty(m),
          continue;
        end
        if isempty(m{3}) || isempty(m{4}),
          continue;
        end
        curexp = find(strcmp(m{1},obj.expnames));
        if isempty(curexp)
          if ~any(strcmp(m{1},missingexps))
            missingexps{end+1} = m{1}; %#ok<AGROW>
          end
          continue;
        end
        dat{end+1,1} = curexp; %#ok<AGROW>
        dat(end,2:4) = cellfun(@(x) str2double(x),m(2:end),'UniformOutput',false);
      end
      
      if ~isempty(missingexps)
        expstring = '';
        for ii = 1:numel(missingexps),
          expstring = [expstring ' ' missingexps{ii}]; %#ok<AGROW>
        end
        uiwait(warndlg(sprintf(...
          'Experiments:%s are not currently loaded. Not loading the GT suggestions for these experiments',...
          expstring)));        
      end
      
      dat = cell2mat(dat);
      fclose(fid);
      expi = dat(:,1); fly = dat(:,2); t0s = dat(:,3); t1s = dat(:,4);
      for curexpi = 1:obj.nexps
        for ndx = 1:obj.nflies_per_exp(curexpi)
          loc = ismember(fly,ndx) & ismember(expi,curexpi);
          if ~any(loc),
            obj.loadedGTSuggestions{curexpi}(ndx).start = 1;
            obj.loadedGTSuggestions{curexpi}(ndx).end = 0;
          else
            obj.loadedGTSuggestions{curexpi}(ndx).start = t0s(loc);
            obj.loadedGTSuggestions{curexpi}(ndx).end = t1s(loc);
          end
        end
      end
      obj.GTSuggestionMode = 'Imported';
      obj.needsave = true;
      success = true;
    end

    
    % ---------------------------------------------------------------------
    function SuggestThresholdGT(obj,threshold)
      obj.thresholdGTSuggestions = threshold;
      obj.GTSuggestionMode = 'Threshold';
    end
    
  end
    
  %% Misc
  
  methods
    
    function obj = JLabelData(varargin)
    % obj = JLabelData(configParams,...)
    %
    % constructor: first input should be the config params. All other
    % inputs are optional. if configfilename is not input, throws an error. 
    % 
    % optional inputs: 
    %
    % TODO: debug this
    % override stuff set in the config file: 
    %
    % moviefilename, trxfilename, perframedir, clipsdir: names of
    % files within experiment directories: 
    % featureparamsfilename: file containing feature parameters
    % rootoutputdir: in case we don't want to write to the experiment
    % directory, we will mirror the experiment directory structure in the
    % rootoutputdir this can be the same as the input root directory
    %
    % defaultpath: default location to look for experiments
    % setstatusfn: handle to function that inputs sprintf-like input and
    % outputs the corresponding string to a status bar.
    % clearstatusfn: handle to function that returns the status to the
    % default string
    % classifierfilename: name of classifier file to save/load classifier from
 
      % Initialize the object
      obj.initialize();
    
      % args should be key-value pairs
      if mod(numel(varargin),2) ~= 0,
        error('JLabelData:oddNumberOfArgsToConstructor',  ...
              'Number of inputs to JLabelData constructor must be even');
      end

      % parse arguments into keywords and corresponding values
      keys = varargin(1:2:end);
      values = varargin(2:2:end);     
      
      % Set the function to be called when the SetStatus method is invoked
      i = find(strcmpi(keys,'setstatusfn'),1);
      if isempty(i),
        obj.setstatusfn = @disp;  % do-nothing function
      else
        obj.setstatusfn = values{i};
      end

      % Set the function to be called when the ClearStatus method is
      % invoked
      i = find(strcmpi(keys,'clearstatusfn'),1);
      if isempty(i),
        obj.clearstatusfn = @()([]);  % do-nothing function
      else
        obj.clearstatusfn = values{i};
      end
                          
      % default path
      i = find(strcmpi(keys,'defaultpath'),1);
      if ~isempty(i),
        [success,msg] = obj.SetDefaultPath(values{i});
        if ~success,
          warning(msg);
        end
      end
      
      % isInteractive
      i = find(strcmpi(keys,'isInteractive'),1);
      if ~isempty(i),
        obj.isInteractive=values{i};
      end
      
      % cacheSize
      % Not used for anything as of Apr 30, 2013 --ALT
      i = find(strcmpi(keys,'cacheSize'),1);
      if ~isempty(i),
        obj.cacheSize = values{i};
      end
      
      % Set the JAABA version
      try 
        
        if isdeployed,
          verfilename = deployedRelative2Global('version.txt');
        else
          verfilename = 'version.txt';
        end
        vid = fopen(verfilename,'r');
        vv = textscan(vid,'%s');
        fclose(vid);
        obj.version = vv{1}{1};
      catch ME
        warning('Cannot detect JAABA Version. Setting it to 0.0. Error:\n%s',getReport(ME));  
        obj.version = '0.0';
      end
      
      % % initialize the status table describing what required files exist
      % [success,msg] = obj.UpdateStatusTable();
      % if ~success,
      %   error(msg);
      % end      
    end  % constructor method

  % Status display    

    % ---------------------------------------------------------------------
    function SetStatusFn(obj,statusfn)
      obj.setstatusfn = statusfn;
    end

    
    % ---------------------------------------------------------------------
    function SetClearStatusFn(obj,clearfn)
      obj.clearstatusfn = clearfn;
    end   
    

    % ---------------------------------------------------------------------
    function SetStatus(obj,varargin)
      % SetStatus(obj,<sprintf-like arguments>)
      % Update an associated status text according to the input sprintf-like
      % arguments.
      
      if isempty(obj.setstatusfn),
        fprintf(varargin{:});
        fprintf('\n');
      else
        obj.setstatusfn(sprintf(varargin{:}));
        drawnow;
      end
      %       allF = findall(0,'type','figure');
      %       jfigNdx = find(strcmp(get(allF,'name'),'JAABA'));
      %       jfig = allF(jfigNdx);
      %       if ~isempty(jfig),
      %         set(jfig,'pointer','watch');
      %       end
    end  % method
    
    
    % ---------------------------------------------------------------------
    function ClearStatus(obj)
      % ClearStatus(obj)
      % Return an associated status text to the default.
      
      if ~isempty(obj.clearstatusfn),
        obj.clearstatusfn();
        drawnow;
      end
    end  % method

    
    
    % ---------------------------------------------------------------------
    function [success,msg] = SetDefaultPath(obj,defaultpath)
    % [success,msg] = SetDefaultPath(obj,defaultpath)
    % sets the default path to load experiments from. only checks for
    % existence of the directory.
      
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
    

    function [success,msg] = SetExpDefaultPath(obj,defaultpath)
    % [success,msg] = SetDefaultPath(obj,defaultpath)
    % sets the default path to load experiments from. only checks for
    % existence of the directory.
      
      success = false;
      msg = '';
      
      if ischar(defaultpath),
        
        if ~isempty(defaultpath) && ~exist(defaultpath,'file'),
          msg = sprintf('defaultpath directory %s does not exist',defaultpath);
          return;
        end
          
        obj.expdefaultpath = defaultpath;
        success = true;
      end

    end
    
    
    % ---------------------------------------------------------------------
    function tf = classifierIsPresent(obj)
      % For multiclassifier projects, returns true only if all classifiers
      % present      
      
      cls = obj.classifier;
      tf = ~isempty(cls) && ~any(cellfun(@isempty,cls));
    end
    
    
    % ---------------------------------------------------------------------
    function someExperimentIsCurrent = getSomeExperimentIsCurrent(self)
      if self.thereIsAnOpenFile
        nExp = self.nexps;
        someExperimentIsCurrent = (1<=self.expi) && ...
                                  (self.expi<=nExp);
      else
        someExperimentIsCurrent = false;
      end
    end
    
   
    % ---------------------------------------------------------------------
    function setToValue(self,jld)
      % Makes self into a clone of jld.  I.e. sets all the independent
      % properties of self to have the same _value_ as the corresponding
      % property of jld.
      mc=metaclass(self);
      propertyList=mc.PropertyList;
      for i=1:length(propertyList)
        property=propertyList(i);
        if ~property.Dependent
          self.(property.Name)=jld.(property.Name);
        end
      end
    end  % method    
        
    
    function isRandom = getColorAssignment(self)
      if isfield(self.trxGraphicParams,'assignment') && ...
          strcmp(self.trxGraphicParams.assignment,'static')
        isRandom = false;
      else
        isRandom = true;
      end
      
    end   

    
    % ---------------------------------------------------------------------
    function setNeedSave(obj)
      % Get the behavior name, a string
      obj.needsave = true;
    end
    
    
    % ---------------------------------------------------------------------
    function delete(~)
    end
    
    
  end
  
  methods % more private
    
    % ---------------------------------------------------------------------
    function initialize(self)
      % Set all properties to the state they should have just after the
      % JLabelData object is created, before any .jab file has been loaded
      % or a new file is created.  The state of the object after calling
      % this function does not depend at all on the state prior to calling
      % it.  Nothing is spared.  Scorched earth.
      % Also: No attempt is done to close open files, or anything like
      % that.  It just sets the properties to the values they should have
      % initially.
      self.expi = 0;
      self.flies = [];    
      self.targettype = 'fly';
      self.trx = {};
      self.perframedata = {};
      self.featureLexiconName='';
      self.featureLexicon=[];
      self.windowdata = WindowData.windowdata(0);
      self.selFeatures = SelFeatures.createEmpty();
      self.predictdata = {};
      self.predictblocks = Predict.predictblocks(0);
      self.fastPredict = Predict.fastPredict(0);
      self.windowdatachunk_radius = 100;
      self.predictwindowdatachunk_radius = 2000;
      self.labels = Labels.labels(0); 
      self.labelidx = struct('vals',[],'imp',[],'timestamp',[]);
      self.labelidx_off = 0;
      self.t0_curr = 0;
      self.t1_curr = 0;
      self.predictedidx = [];
      self.scoresidx = [];  
      self.scoresidx_old = [];
      self.scoreTS = [];  
      %self.erroridx = [];
      self.labelnames = {};
      %self.nbehaviors = 0;
      self.ntimelines = 0;
%       self.labelstats = struct('nflies_labeled',{},'nbouts_labeled',{});
      self.perframe_params = {};
      self.landmark_params = struct;  % scalar struct with no fields
      self.classifiertype = 'boosting';
      self.classifier = [];
      self.classifier_old = [];
      self.lastFullClassifierTrainingSize = 0;
      self.classifierTS = zeros(1,0); 
      self.trainstats = cell(1,0);
      self.classifier_params = cell(1,0);
      self.filetypes = {'movie','trx','perframedir','clipsdir','scores','stfeatures','trk'};
      self.moviefilename = 0;
      self.trxfilename = 0;
      self.scorefilename = 0;
      self.perframedir = 0;
      self.clipsdir = 0;
      self.scores = 0;
      self.stfeatures = 0;
      self.trkfilename = 0;
      %self.openmovie = false;
      % self.ismovie = false;
      self.expdirs = {};
      self.nflies_per_exp = [];
      self.firstframes_per_exp = {};
      self.endframes_per_exp = {};
      self.frac_sex_per_exp = {};
      self.sex_per_exp = {};
      self.hassex = false;
      self.hasperframesex = false;
      self.defaultpath = '';
      self.expdefaultpath = '';
      self.windowfeaturesparams = cell(0,1);
      self.windowfeaturescellparams = cell(0,1);
      self.allperframefns = {};
      self.curperframefns = {};
      self.perframeunits = {};
      self.scoreFeatures = struct('classifierfile',{}, ...
                                  'ts',{}, ...
                                  'scorefilename',{});
      self.fileexists = false(0,numel(self.filetypes));
      self.filetimestamps = nan(0,numel(self.filetypes));
      self.allfilesexist = true;
      self.filesfixable = true;
      self.perframeGenerate = [];
      self.perframeOverwrite = [];
      self.arenawarn = true;
      self.hasarenaparams = [];
      self.setstatusfn = '';
      self.clearstatusfn = '';
      self.frameFig = [];
      self.distMat = [];
      self.bagModels = {};
      self.fastPredictBag = ...
        struct('classifier',[], ...
               'windowfeaturescellparams',[],...
               'wfs',[], ...
               'pffs',[], ...
               'ts',[], ...
               'tempname',[], ...
               'curF',[], ...
               'dist',[], ...
               'trainDist',[]);
      self.confThresholds = zeros(0,2);
      self.doUpdate = true;
      self.predictOnlyCurrentFly = false(1,0);

      self.gtMode = [];
      self.randomGTSuggestions = {};
      self.thresholdGTSuggestions = [];
      self.loadedGTSuggestions = {};
      self.balancedGTSuggestions = {};
      self.GTSuggestionMode = '';
      self.cacheSize = 4000;
      self.postprocessparams = [];
      self.version = '';
      self.otherModeLabelsEtc = struct('expDirNames',{cell(1,0)}, ...
                                       'labels',{struct([])});
      self.trxGraphicParams=[];
      self.labelcolors = [];
      self.unknowncolor = [0 0 0];
      self.isInteractive=true;
      self.thereIsAnOpenFile=false;
      self.everythingFileNameAbs='';
      self.userHasSpecifiedEverythingFileName=false;
      self.needsave=false;
      self.savewindowdata = false(1,0);
      self.loadwindowdata = true(1,0);
      
      self.usePastOnly = false;
      
    end  % method
    
    % ---------------------------------------------------------------------
    function ClearCachedPerExpData(obj)
    % ClearCachedPerExpData(obj)
    % Clears all cached data for the currently loaded experiment
      obj.unsetCurrentTarget();
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
      %obj.erroridx = [];
    end
    
  end
  
  methods (Static)
    
    function params = convertTransTypes2Cell(params)
      % Convert the trans_types field into cell type
      if ~isstruct(params), return; end
      fnames = fieldnames(params);
      for ndx = 1:numel(fnames)
        if isstruct(params.(fnames{ndx}))
          params.(fnames{ndx}) = JLabelData.convertTransTypes2Cell(params.(fnames{ndx}));
        end
      end
      if isfield(params,'trans_types')
        if ~iscell(params.trans_types)
          params.trans_types = {params.trans_types};
        end
        x = cellfun(@isempty,params.trans_types);
        params.trans_types(x)=[];
      end
    end
    
    
    % ---------------------------------------------------------------------
    function cellparams = convertParams2CellParams(params)
      % Convert the windowFeatureParams structure, which stores a 'list' of
      % per-frame features and the parameters that determine how each is
      % converted to a set of window features, to it's corresponding
      % cell-based form, which is apparently sometimes useful.
      cellparams = struct;
      fns1 = fieldnames(params);
      for i1 = 1:numel(fns1),
        fn1 = fns1{i1};
        fns2 = fieldnames(params.(fn1));
        cellparams.(fn1) = {};
        feature_types = {};
        for i2 = 1:numel(fns2),
          fn2 = fns2{i2};
          if ~isstruct(params.(fn1).(fn2)),
            cellparams.(fn1)(end+1:end+2) = {fn2,params.(fn1).(fn2)};
          else
            cellparams.(fn1)(end+1:end+2) = {[fn2,'_params'],struct2paramscell(params.(fn1).(fn2))};
            feature_types{end+1} = fn2; %#ok<AGROW>
          end
        end
        cellparams.(fn1)(end+1:end+2) = {'feature_types',feature_types};
      end
    end  % method
  
  end
   

end
