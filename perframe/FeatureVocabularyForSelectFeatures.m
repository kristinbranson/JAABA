classdef FeatureVocabularyForSelectFeatures < handle
  
  % instance variables
  properties
%     jld  % instance of JLabelData, a ref
%     scoreFeatures  % the scores-as-input structure, detailing the files that
%                    % hold classifiers to be used as inputs, and their
%                    % per-frame feature names
%     windowFeatureParams
%    featureLexicon  % the dictionary of possible per-frame features
    wfParamsFromAmount  % scalar structure holding the window feature parameters
                        % that go with each of the possible "amount" levels.
                        % e.g. 'normal', 'more', 'less'.
    pfCategoriesFromName  % scalar structure.  Each field is the name of a 
                          % per-frame feature, and the value is a cell
                          % array of strings, telling what categories that
                          % per-frame feature is in.  Most PFs are in only
                          % one category, but some are in more than one.
    pfCategoryNames  % a cell array of strings, holding the names of all the
                     % possible per-frame feature categories
    subdialectPFNames
      % A cell array of strings, containing the names of all the 
      % to-be-calculated per-frame features.  This includes the 
      % scores-as-input PFs.  I like to call this the "subdialect".
    defaultPFTransTypesFromName  
      % A structure, each field is a per-frame feature, and each value is a cell 
      % array of strings holding the default translation types
      % for that per-frame feature  (formerly transType)
    vocabulary
      % a structure where each field is a per-frame feauture, and the value
      % holds the window feature parameters for that PF.
      % Should have as many fields as there are per-frame features in the
      % to-be-calculated vocabulary (the subdialect).  vocabulary.(perframeFeatureName).enabled is a boolean
      % indicating whether that per-frame feature is in the vocabulary.
      % (formerly data)
  end  % properties
  
  % things that look like instance vars from the outside, but
  % which are calculated on-the-fly
  properties (Dependent=true)
    wfAmounts
  end
  
  % class constants
  properties (Access=public,Constant=true)
    % Window features come in different types.  Below is a list of the
    % different types.  (formerly windowComp)
    wfTypes = {'mean','min','max','hist','prctile',...
               'change','std','harmonic','diff_neighbor_mean',...
               'diff_neighbor_min','diff_neighbor_max','zscore_neighbors'};
    % Window features have parameters.  Below are the names of those parameters.
    % These parameters exist for each WF type.  (They should arguably be
    % called the "standard" parameters, or the "common" parameters,
    % because they exist for all WF types.  Whereas the extra parameters,
    % below, are each particular to one WF type.
    % (Formerly winParams.)
    wfParamNames = {'max_window_radius','min_window_radius','nwindow_radii', ...
                    'trans_types','window_offsets'};
%     % These are the default values for the window feature parameters.  (Formerly defaultWinParams.)
%     % They are only used if the lexicon fails to specify pre-set
%     % window-feature amounts (e.g. normal, more, less)
%     defaultWFParams = {10,1,3,{'none'},0};
    % A Window feature type can sometimes have an extra parameter.  Below are the names 
    % of those parameters, one for each WF type.  WF types without an extra param have the 
    % empty string for their entry.  (Formerly winextraParams.)                  
    wfExtraParamNames = {'','','','hist_edges','prctiles','change_window_radii', ...
                         '','num_harmonic','','','',''};
%     % Below are the default values for the extra window feature parameters.  (Formerly winextraDefaultParams.)              
%     % They are only used if the lexicon fails to specify pre-set
%     % window-feature amounts (e.g. normal, more, less)
%     defaultWFExtraParams = {[],[],[],[],[-400000 0 40000],[5 10 30 50 70 90 95],[1 3],...
%                             [],2,[],[],[],[]};
  end  % class constants
  
  % class methods
  methods (Access=public,Static=true)
    % ---------------------------------------------------------------------
    function wfParamsFromAmount=computeWFParamsFromAmount(featureLexicon)
      % Specify some constants we need
      % Below are the default values for the window feature parameters.  (Formerly defaultWinParams.)
      % They are only used if the lexicon fails to specify pre-set
      % window-feature amounts (e.g. normal, more, less)
      defaultWFParams = {10,1,3,{'none'},0};
      % Below are the default values for the extra window feature parameters.  (Formerly winextraDefaultParams.)              
      % They are only used if the lexicon fails to specify pre-set
      % window-feature amounts (e.g. normal, more, less)
      defaultWFExtraParams = {[],[],[],[-400000 0 40000],[5 10 30 50 70 90 95],[1 3],...
                              [],2,[],[],[],[]};
    
      % Get some class constants we'll need
      wfTypes = ...
        FeatureVocabularyForSelectFeatures.wfTypes;
      wfParamNames = ...
        FeatureVocabularyForSelectFeatures.wfParamNames;
      wfExtraParamNames = ...
        FeatureVocabularyForSelectFeatures.wfExtraParamNames;
      % Read the default parameters for different categories.
      wfAmounts = fieldnames(featureLexicon.defaults);
      if isempty(wfAmounts) 
        % If no preset amounts have been specified in the lexicon file,
        % fall back on the built-in defaults.
        wfAmount='default';
        wfParams = struct;
        % % Set the parameters for the WFs of type 'default'
        % wfParams.default.enabled = true;
        % for ndx = 1:numel(wfParamNames)
        %   wfParamName = wfParamNames{ndx};
        %   wfParams.default.values.(wfParamName) = defaultWFParams{ndx};
        % end
        %wfParams.default.values.sanitycheck = false;
        % Set the parameters for the other WF types
        for ndx = 1:numel(wfTypes)
          wfType = wfTypes{ndx};
          wfParams.(wfType).enabled = false;
          for paramIndex = 1:numel(wfParamNames)
            wfParamName = wfParamNames{paramIndex};
            wfParams.(wfType).values.(wfParamName) = defaultWFParams{paramIndex};
          end
          if ~isempty(wfExtraParamNames{ndx})
            extraParamName = wfExtraParamNames{ndx};
            wfParams.(wfType).values.(extraParamName) = defaultWFExtraParams{ndx};
          end
        end
        wfParamsFromAmount.(wfAmount) = wfParams;
      else
        % The usual case --- preset amounts are defined in the lexicon.
        for wfAmountIndex = 1:numel(wfAmounts)
          wfAmount=wfAmounts{wfAmountIndex};
          wfParams = struct;
          wfParamsThisAmount = featureLexicon.defaults.(wfAmount);
          % % Set the parameters for the WFs of type 'default'
          % wfParams.default.enabled = true;
          % for ndx = 1:numel(wfParamNames)
          %   wfParamName = wfParamNames{ndx};
          %   wfParams.default.values.(wfParamName) = wfParamsThisAmount.(wfParamName);
          % end
          % wfParams.default.values.sanitycheck = false;
          % Set the parameters for the other WF types
          for ndx = 1:numel(wfTypes)
            % Copy the standard WF params for this WF type
            wfType = wfTypes{ndx};
            for paramIndex = 1:numel(wfParamNames)
              wfParamName = wfParamNames{paramIndex};
              wfParams.(wfType).values.(wfParamName) = wfParamsThisAmount.(wfParamName);
            end
            % Copy the extra WF parameter, if present for this WF type
            if ~isempty(wfExtraParamNames{ndx})
              extraParamName = wfExtraParamNames{ndx};
              if isfield(wfParamsThisAmount,wfType) && isfield(wfParamsThisAmount.(wfType),extraParamName)
                % Use a value from the preset, if present
                wfParams.(wfType).values.(extraParamName) = wfParamsThisAmount.(wfType).(extraParamName);
              else
                % If not present, fall back to the default default
                wfParams.(wfType).values.(extraParamName) = defaultWFExtraParams{ndx};
              end
            end
            % If this WF type is included in the lexicon preset, enable the WF type
            wfParams.(wfType).enabled=isfield(wfParamsThisAmount,wfType);
            % If this WF type is included in the lexicon preset, check for any WF
            % parameter names besides the usual standard and extra ones,
            % and include them in the processed preset parameters.
            % This seems weird.  Maybe we can leave it out?
            % if isfield(wfParamsThisAmount,wfType)
            %   wfParamNamesThisAmount = fieldnames(wfParamsThisAmount.(wfType));
            %   for dndx = 1:numel(wfParamNamesThisAmount)
            %     if ismember(wfParamNamesThisAmount{dndx},wfParamNames) ,
            %       wfParamName = wfParamNamesThisAmount{dndx};
            %       wfParams.(wfType).values.(wfParamName) = wfParamsThisAmount.(wfType).(wfParamName);
            %     end
            %   end
            % end
          end
          wfParamsFromAmount.(wfAmount) = wfParams;
        end
      end
    end  % method    
    
    
    % ---------------------------------------------------------------------
    function defaultPFTransTypesFromName= ...
        defaultPFTransTypesFromNameFromFeatureLexicon(featureLexicon,scoreFeatures)
      % Calculate the structure that maps PF names to transformation types.
      
      % Get the list of all the per-frame features in the lexicon,
      % including the scores-as-inputs ones
      lexiconPFNames=fieldnames(featureLexicon.perframe);
      scorePFNames = {scoreFeatures(:).scorefilename};
      subdialectPFNames=[lexiconPFNames;scorePFNames];

      % Make one entry in the result structure per PF
      defaultPFTransTypesFromName = struct;
      for i = 1:length(subdialectPFNames)
        pfName = subdialectPFNames{i};
        if ismember(pfName,scorePFNames) ,
          transTypesThis = {'none'};
        else  
          transTypesThisRaw = featureLexicon.perframe.(pfName).trans_types;
          if ischar(transTypesThisRaw)
            transTypesThis = {transTypesThisRaw};
          else
            transTypesThis = transTypesThisRaw;
          end  
        end
        defaultPFTransTypesFromName.(pfName)=transTypesThis;
      end  % for
    end  % method
    
    
    % ---------------------------------------------------------------------
    function pfCategoriesFromName= ...
        pfCategoriesFromNameFromFeatureLexicon(featureLexicon,scoreFeatures)
      % Calculate the structure that maps a PF name to the PF categories
      % that it belongs to.
      
      % Get the list of all the per-frame features in the lexicon,
      % including the scores-as-inputs ones
      lexiconPFNames=fieldnames(featureLexicon.perframe);
      scorePFNames = {scoreFeatures(:).scorefilename};
      subdialectPFNames=[lexiconPFNames;scorePFNames];

      % Make one entry in the result structure per PF
      pfCategoriesFromName = struct;
      for i = 1:length(subdialectPFNames)
        pfName = subdialectPFNames{i};
        if ismember(pfName,scorePFNames) ,
          pfCategoriesThis = {'scores'};
        else  
          pfCategoriesThisRaw = featureLexicon.perframe.(pfName).type;
          if ischar(pfCategoriesThisRaw)
            pfCategoriesThis = {pfCategoriesThisRaw};
          else
            pfCategoriesThis = pfCategoriesThisRaw;
          end  
        end
        pfCategoriesFromName.(pfName)=pfCategoriesThis;
      end  % for
    end  % method
    
    
    % ---------------------------------------------------------------------
    function pfCategoryNames=pfCategoryNamesFromFeatureLexicon(featureLexicon)
      fallpf = fieldnames(featureLexicon.perframe);
      pfCategoryNames = {};
      for pfndx = 1:numel(fallpf)
        curpf = fallpf{pfndx};
        curtypes  = featureLexicon.perframe.(curpf).type; 
        if ischar(curtypes)
          curT = curtypes;
          if ~any(strcmp(pfCategoryNames,curT))
            pfCategoryNames{end+1} = curT;  %#ok
          end
        else    
          for tndx = 1:numel(curtypes)
            curT = curtypes{tndx};
            if ~any(strcmp(pfCategoryNames,curT))
              pfCategoryNames{end+1} = curT;  %#ok
            end
          end
        end
      end
      if ~any(strcmp(pfCategoryNames,'scores')),
        pfCategoryNames{end+1} = 'scores';
      end
    end  % method
    
  end  % class methods
  
  
  % -----------------------------------------------------------------------  
  % public instance methods
  methods
    % ---------------------------------------------------------------------
    function self=FeatureVocabularyForSelectFeatures(featureLexicon, ...
                                                     scoreFeatures, ...
                                                     subdialectPFNames, ...
                                                     windowFeatureParams, ...
                                                     maxWindowRadiusCommon)
      % Constructor.  jld an instance of JLabelData, maxWindowRadiusCommon
      % a value to be set for max_window_radius for all PFs, all WF types,
      % and all preset WF amounts.
      
      % Set a bunch of things that are derived from the constructor input,
      % but stay constant for the life of the object
      self.wfParamsFromAmount= ...
        FeatureVocabularyForSelectFeatures.computeWFParamsFromAmount(featureLexicon);
      self.defaultPFTransTypesFromName= ...
        FeatureVocabularyForSelectFeatures.defaultPFTransTypesFromNameFromFeatureLexicon(featureLexicon, ...
                                                                                         scoreFeatures);
      self.pfCategoriesFromName= ...
        FeatureVocabularyForSelectFeatures.pfCategoriesFromNameFromFeatureLexicon(featureLexicon, ...
                                                                                  scoreFeatures);
      self.pfCategoryNames= ...
        FeatureVocabularyForSelectFeatures.pfCategoryNamesFromFeatureLexicon(featureLexicon);
      
      % Populate the list of per-frame feature names
      self.subdialectPFNames = subdialectPFNames;
      
      % Populate the feature vocabulary proper
      enabledPFs = fieldnames(windowFeatureParams);
      wfParamNames = self.wfParamNames;
      wfAmounts=self.wfAmounts;
      wfAmount=wfAmounts{1};  % by convention, the first one is the default setting (e.g. 'normal')
      wfParamsAmount=self.wfParamsFromAmount.(wfAmount);
      nToBeCalculatedPFNames=length(subdialectPFNames);
      vocabulary = cell(1,nToBeCalculatedPFNames);
      for pfIndex = 1:nToBeCalculatedPFNames
        pfName = subdialectPFNames{pfIndex};
        vocabulary{pfIndex}.name = pfName;
        if ismember(pfName,enabledPFs) ,
          % if this PF is enabled, copy its params out of
          % windowFeatureParams
          vocabulary{pfIndex}.enabled = true;  
          vocabulary{pfIndex}.sanitycheck = windowFeatureParams.(pfName).sanitycheck;
          % Fill in the fallback values for this PF.  These are used if 
          fallbackValues=wfParamsAmount.default.values;
%           fallbackValues=struct();
%           for wfParamNdx = 1:numel(wfParamNames)
%             wfParamName = wfParamNames{wfParamNdx};
%             fallbackValues.(wfParamName) = ...
%               wfParamsAmount.default.values.(wfParamName);
%             end
%           end
          % Fill for the rest of the window feature types.
          enabledWFTypes = fieldnames(windowFeatureParams.(pfName));
          for wfTypeNdx = 1:numel(self.wfTypes)
            wfType = self.wfTypes{wfTypeNdx};
            if ismember(wfType,enabledWFTypes) ,
              vocabulary{pfIndex}.(wfType) = ...
                wfParamsOfOneType(windowFeatureParams.(pfName).(wfType), ...
                                  fallbackValues, ...
                                  wfParamNames, ...
                                  self.wfExtraParamNames{wfTypeNdx});
            else
              % If the window type is disabled, use default values for its
              % parameters
              vocabulary{pfIndex}.(wfType).enabled = false;
              for wfParamNdx = 1:numel(wfParamNames)
                wfParamName = wfParamNames{wfParamNdx};
                vocabulary{pfIndex}.(wfType).values.(wfParamName) = fallbackValues.(wfParamName);
              end
              if ~isempty(self.wfExtraParamNames{wfTypeNdx})
                extraParamName = self.wfExtraParamNames{wfTypeNdx};
                %vocabulary{pfIndex}.(wfType).values.(extraParamName) = self.defaultWFExtraParams{wfTypeNdx};
                vocabulary{pfIndex}.(wfType).values.(extraParamName) = wfParamsAmount.(wfType).values.(extraParamName);
              end
            end
          end
        else
          % If PF is disabled, use values in the preset
          vocabulary{pfIndex}=wfParamsAmount;
          vocabulary{pfIndex}.enabled = false;
          vocabulary{pfIndex}.sanitycheck = false;
          % Tack on the default translation types, which depend on which
          % per-frame feature this is
          for i=1:length(self.wfTypes)
            wfType=self.wfTypes{i};
            vocabulary{pfIndex}.(wfType).values.trans_types=self.defaultPFTransTypesFromName.(pfName);
          end
        end
      end
      self.vocabulary = vocabulary;

      % If a common max window radius is set, use it
      if exist('maxWindowRadiusCommon','var') && ~isempty(maxWindowRadiusCommon)
        % set the global maxWindowRadius to the given value
        self.setMaxWindowRadiusForAllWFs(maxWindowRadiusCommon);
        self.setMaxWindowRadiusForAllWFAmounts(maxWindowRadiusCommon);
      end
    end  % constructor method
    
    
    % ---------------------------------------------------------------------
    function enableAllPFsInCategory(self,pfCategoryIndex)
      %pfCategory=self.pfCategoryNames{pfCategoryIndex};
      %subdialectPFNames=self.subdialectPFNames;
      %pfCategoriesFromName=self.pfCategoriesFromName;
      pfNamesThisCategory=self.getPFNamesInCategory(pfCategoryIndex);
      % Iterate over the per-frame features, looking for ones that are
      % within the selected category.
      for i = 1:length(pfNamesThisCategory)
        pfName=pfNamesThisCategory{i};
        self.enablePerframeFeature(pfName)
      end
    end  % method

    
    % ---------------------------------------------------------------------
    function setAllPFsInCategoryToWFAmount(self,pfCategoryIndex,wfAmount)
      % This does what it says, but note that it doesn't add or subtract
      % any PFs from the vocabulary.
      %pfCategory=self.pfCategoryNames{pfCategoryIndex};
      %subdialectPFNames=self.subdialectPFNames;
      %pfCategoriesFromName=self.pfCategoriesFromName;
      pfNamesThisCategory=self.getPFNamesInCategory(pfCategoryIndex);
      % Iterate over the per-frame features, looking for ones that are
      % within the selected category.
      for i = 1:length(pfNamesThisCategory)
        pfName=pfNamesThisCategory{i};
        self.setPFToWFAmount(pfName,wfAmount)
      end
    end  % method

    
    % ---------------------------------------------------------------------
    function setPFEnablement(self,pfSpecifier,enabled)
      % Turn on/off the given per-frame feature, but leave the window feature 
      % parameters alone.
      % pfIndicator can be a name or an index.
      if ischar(pfSpecifier)
        % this means pfIndex is really a name
        pfName=pfSpecifier;
        pfIndex=find(strcmp(pfName,self.subdialectPFNames));  
      else
        pfIndex=pfSpecifier;
      end
      self.vocabulary{pfIndex}.enabled = enabled;
    end  % method

    
    % ---------------------------------------------------------------------
    function enablePerframeFeature(self,pfSpecifier)
      % Turn on the given per-frame feature, but leave the window feature 
      % parameters alone.
      % This method is deprecated.
      setPFEnablement(self,pfSpecifier,true);
    end  % method
    
    
    % ---------------------------------------------------------------------
    function disablePerframeFeature(self,pfSpecifier)
      % Turn off the given per-frame feature.  (I.e. remove all it's window 
      % features from the vocabulary)
      % This method is deprecated.
      setPFEnablement(self,pfSpecifier,false);
    end  % method
    
    
    % ---------------------------------------------------------------------
    function setPFToWFAmount(self,pfSpecifier,wfAmount)
      % Set the window feature parameters for the named PF to those
      % specified by the named wfAmount (normal, more, less).  This is
      % orthogonal to whether that PF is enabled or not.  Note that this
      % also sets the transformation types for the PF to the default
      % set.
      if isequal(wfAmount,'custom')
        return
      end
      if ischar(pfSpecifier)
        pfName=pfSpecifier;
        pfIndex=find(strcmp(pfName,self.subdialectPFNames));
      else
        pfIndex=pfSpecifier;
        pfName=self.subdialectPFNames{pfIndex};
      end
      defaultPFTransTypes=self.defaultPFTransTypesFromName.(pfName);
      %self.vocabulary{pfIndex}.enabled = true;
      %pfCategories=self.pfCategoriesFromName.(pfName);
      %pfCategory = pfCategories{1};  % why the first one?
      %categoryNdx = find(strcmp(pfCategory,self.pfCategoryNames));
      wfParamTemplate=self.wfParamsFromAmount.(wfAmount);
      %self.vocabulary{pfIndex}.enabled = true;
      for wfTypeIndex = 1:numel(self.wfTypes)
        wfType = self.wfTypes{wfTypeIndex};
        %wfParamsThisType=wfParamTemplate.(wfType);
        %if ~wfParamsThisType.enabled
        %  self.vocabulary{pfIndex}.(wfType).enabled = false;
        %  continue;
        %end
        self.setWindowFeaturesOfSingleType(pfIndex,wfTypeIndex,wfParamTemplate);
        self.vocabulary{pfIndex}.(wfType).values.trans_types = defaultPFTransTypes;
      end
    end  % method
    
    
    % ---------------------------------------------------------------------
    function setWFTypeEnablement(self,pfName,wfType,enable)
      % For the given per-frame feature name, enable/disable the given
      % window-feature type
      pfIndex=find(strcmp(pfName,self.subdialectPFNames));
      self.vocabulary{pfIndex}.(wfType).enabled = enable;  %#ok
    end
    
    
    % ---------------------------------------------------------------------
    function setWFParam(self,pfName,wfType,wfParamName,newValue)
      % For the given per-frame feature name, window-feature type, and
      % window-feature parameter name, set the parameter to newValue
      if ischar(pfName)
        % the usual case
        pfIndex=find(strcmp(pfName,self.subdialectPFNames));
      else
        % this means pfName is really an index
        pfIndex=pfName;
        %pfName=self.subdialectPFNames{pfIndex};
      end
      if ~ischar(wfType)
        % this means wfType is really an index
        wfIndex=wfType;
        wfType=FeatureVocabulary.wfTypes{wfIndex};
      end      
      self.vocabulary{pfIndex}.(wfType).values.(wfParamName) = newValue;
    end

    
    % ---------------------------------------------------------------------
    function addWFTransformation(self,pfName,wfType,newTransformation)
      % For per-frame feature pfName, window feature type wfType, add a new
      % transformation, newTransformation.  (A transformation is one of 'none',
      % 'flip', 'abs', 'relative'.)
      pfIndex=find(strcmp(pfName,self.subdialectPFNames));
      transformations = self.vocabulary{pfIndex}.(wfType).values.trans_types;     
      if ~any(strcmp(newTransformation,transformations))
        % if it's not in the list, add it to the end
        transformations{end+1}=newTransformation;
        self.vocabulary{pfIndex}.(wfType).values.trans_types=transformations;
      end
    end  % method
    
    
    % ---------------------------------------------------------------------
    function removeWFTransformation(self,pfName,wfType,transformation)
      % For per-frame feature pfName, window feature type wfType, remove a 
      % transformation, transformation.  (A transformation is one of 'none',
      % 'flip', 'abs', 'relative'.)  If transformation is not in the
      % transformation list for that window feature, nothing changes.  If
      % removing the given transformation would leave the transformation
      % list empty, nothing changes.
      pfIndex=find(strcmp(pfName,self.subdialectPFNames));
      transformations = self.vocabulary{pfIndex}.(wfType).values.trans_types;
      toBeRemoved = strcmp(transformation,transformations);
      transformations(toBeRemoved) = [];
      % only commit the change if there is still at least one transformation
      % left
      if ~isempty(transformations)
        self.vocabulary{pfIndex}.(wfType).values.trans_types=transformations;
      end
    end  % method

    
    % ---------------------------------------------------------------------
    function result=isWFTransformation(self,pfName,wfType,transformation)
      % For per-frame feature pfName, window feature type wfType, returns
      % true iff the given transformation is in the transformation list.
      pfIndex=find(strcmp(pfName,self.subdialectPFNames));
      transformations = self.vocabulary{pfIndex}.(wfType).values.trans_types;  %#ok
      result=ismember(transformation,transformations);
    end  % method

    
    % ---------------------------------------------------------------------
    function setMaxWindowRadiusForAllWFAmounts(self,newValue)
      % This sets the maximum window radius for all of the preset window
      % feature amounts, for all window feature types.  Note, however, that
      % it does not change the window feature parameters
      % themselves.  So WF params that were previously set using one of the
      % WF amount pre-sets will retain the old values.
      wfAmounts = self.wfAmounts;
      wfTypes = self.wfTypes;
      for i = 1:numel(wfAmounts)
        wfAmount = wfAmounts{i};
        for j = 1:numel(wfTypes)
          wfType=wfTypes{j};
          wfParams=self.wfParamsFromAmount.(wfAmount);
          if isfield(wfParams,wfType),
            self.wfParamsFromAmount.(wfAmount).(wfType).values.max_window_radius = newValue;
          end
        end
      end
    end  % method
    
    
    % ---------------------------------------------------------------------
    function setMaxWindowRadiusForAllWFs(self,newValue)
      % This sets the maximum window radius for all per-frame features, for
      % all window-feature types.  Thus for all window features.
      pfNames=self.subdialectPFNames;
      wfTypes = self.wfTypes;
      for i = 1:length(pfNames)
        for j = 1:length(wfTypes)
          wfType=wfTypes{j};
          if isfield(self.vocabulary{i},wfType),
            self.vocabulary{i}.(wfType).values.max_window_radius=newValue;
          end
        end
      end
    end  % method
    
    
    % ---------------------------------------------------------------------
    function [thereIsConsensus,consensusValue]=getConsensusMaxWindowRadiusForAllWFs(self)
      % Determines the common value of max_window_radius across all WFs, in
      % all PFs.  thereIsConsensus will be true iff there is a common
      % value.  The common value is returned in consensusValue.  If there
      % is no common value, consensusValue is unspecified.  If there are no
      % WFs at all in the vocabulary, thereIsConsensus will be true, but
      % consensusValue will be empty.
      thereIsConsensus=true;
      consensusValue=[];
      pfNames=self.subdialectPFNames;
      wfTypes = self.wfTypes;
      for i = 1:length(pfNames)
        for j = 1:length(wfTypes)
          wfType=wfTypes{j};
          if isfield(self.vocabulary{i},wfType),
            thisValue=self.vocabulary{i}.(wfType).values.max_window_radius;
            if isempty(consensusValue)
              % If this is the first value we've examined, it becomes the
              % consensus value
              consensusValue=thisValue;
            else
              % if this is not the first value we've examined, compare it
              % to the consensusValue.
              if thisValue~=consensusValue
                % No consensus value---return empty matrix
                thereIsConsensus=false;
                return
              end
            end
          end
        end
      end
    end  % method
    
    
    % ---------------------------------------------------------------------
    function copyWFParams(self,pfNdxFrom,wfTypeNdxFrom,pfNdxTo,wfTypeNdxTo)
      wfTypeFrom = self.wfTypes{wfTypeNdxFrom};
      % something to copy from?
      if ~isfield(self.vocabulary{pfNdxFrom},wfTypeFrom),
        return;
      end
      wfTypeTo = self.wfTypes{wfTypeNdxTo};
      self.vocabulary{pfNdxTo}.(wfTypeTo).enabled = self.vocabulary{pfNdxFrom}.(wfTypeFrom).enabled;
      for i = 1:numel(self.wfParamNames),
        wfParamName = self.wfParamNames{i};
        self.vocabulary{pfNdxTo}.(wfTypeTo).values.(wfParamName) = ...
          self.vocabulary{pfNdxFrom}.(wfTypeFrom).values.(wfParamName);
      end
      if ~isempty(self.wfExtraParamNames{wfTypeNdxTo})
        extraParam = self.wfExtraParamNames{wfTypeNdxTo};
        if isfield(self.vocabulary{pfNdxFrom}.(wfTypeFrom).values,extraParam)
          self.vocabulary{pfNdxTo}.(wfTypeTo).values.(extraParam) = ...
            self.vocabulary{pfNdxFrom}.(wfTypeFrom).values.(extraParam);
        else
          self.vocabulary{pfNdxTo}.(wfTypeTo).values.(extraParam) = '';
        end
      end
    end  % method    
    
    
    % ---------------------------------------------------------------------
    function result=get.wfAmounts(self)
      result=fieldnames(self.wfParamsFromAmount);
    end  % method      

    
    % ---------------------------------------------------------------------
    function result=pfIsInVocabulary(self,pfIndex)
      if ischar(pfIndex)
        % means pfIndex is really a per-frame feature name
        pfName=pfIndex;
        pfIndex=find(strcmp(pfName,self.subdialectPFNames));
      end
      result=self.vocabulary{pfIndex}.enabled;
    end  % method      

    
    % ---------------------------------------------------------------------
    function result=wfTypeIsInVocabulary(self,pfIndex,wfType)
      if ischar(pfIndex)
        % means pfIndex is really a per-frame feature name
        pfName=pfIndex;
        pfIndex=find(strcmp(pfName,self.subdialectPFNames));
      end
      if isnumeric(wfType)
        % means wfType is really a window feature index
        wfIndex=wfType;
        wfType=self.wfTypes{wfIndex};
      end
      result=self.vocabulary{pfIndex}.(wfType).enabled;
    end  % method      
    
    
    % ---------------------------------------------------------------------    
    function windowFeatureParams = getInJLabelDataFormat(self)
      % Converts the feature vocabulary into the format used by JLabelData,
      % returns this.
      % Note that none of the information about window-feature amount
      % presets (e.g. normal, more, less) gets used in computing the
      % output.
      windowFeatureParams = struct;
      for pfIndex = 1:numel(self.subdialectPFNames)
        pfName = self.subdialectPFNames{pfIndex};
        if ~self.vocabulary{pfIndex}.enabled; continue;end
        pfParams = self.vocabulary{pfIndex};
        windowFeatureParams.(pfName).sanitycheck = pfParams.sanitycheck;
        % the default window feature type is always included
        % for wfParamsIndex = 1:numel(self.wfParamNames)
        %   wfParamName = self.wfParamNames{wfParamsIndex};
        %   windowFeatureParams.(pfName).(wfParamName) = pfParams.default.values.(wfParamName);
        % end
        for wfTypeIndex = 1:numel(self.wfTypes)
          wfType = self.wfTypes{wfTypeIndex};
          if pfParams.(wfType).enabled,
            for wfParamsIndex = 1:numel(self.wfParamNames)
              wfParamName = self.wfParamNames{wfParamsIndex};
              windowFeatureParams.(pfName).(wfType).(wfParamName) =...
                  pfParams.(wfType).values.(wfParamName);
            end
            if ~isempty(self.wfExtraParamNames{wfTypeIndex})
              extraParam = self.wfExtraParamNames{wfTypeIndex};
              extraParamVal = pfParams.(wfType).values.(extraParam);
              windowFeatureParams.(pfName).(wfType).(extraParam) = extraParamVal;
            end
          end
        end
      end
    end  % method
    
    
    % ---------------------------------------------------------------------
    function pfNames=getPFNamesInCategory(self,pfCategoryName)
      % For the given per-frame feature category, return a cell array of
      % strings containing the names of all per-frame features in that
      % category.
      if isnumeric(pfCategoryName)
        % this means the category name is really a category index
        pfCategoryIndex=pfCategoryName;
        pfCategoryName=self.pfCategoryNames{pfCategoryIndex};
      end
      pfNames=cell(1,0);
      for iPF = 1:numel(self.subdialectPFNames)
        pfName=self.subdialectPFNames{iPF};
        categoriesThisPF=self.pfCategoriesFromName.(pfName);
        if ismember(pfCategoryName,categoriesThisPF)
          % if the selected category contains this per-frame feature,
          % do stuff
          pfNames{1,end+1}=pfName;  %#ok
        end
      end
    end  % method    
    

    % ---------------------------------------------------------------------
    function level=getPFCategoryLevel(self,pfCategoryName)
      % For the given per-frame feature category, calculate what 'level' of 
      % per-frame features are in the vocabulary.  Returns one of 'none',
      % 'custom', 'all'.
      pfNamesInCategory=self.getPFNamesInCategory(pfCategoryName);
      nPFsInCategory=length(pfNamesInCategory);
      nPFsInCategoryAndVocab=0;
      for i = 1:nPFsInCategory
        pfName=pfNamesInCategory{i};
        if self.pfIsInVocabulary(pfName),
          nPFsInCategoryAndVocab=nPFsInCategoryAndVocab+1;
        end
      end
      if nPFsInCategoryAndVocab==0
        if nPFsInCategoryAndVocab==nPFsInCategory
          level='n/a';
        else
          level='none';
        end
      elseif nPFsInCategoryAndVocab==nPFsInCategory
        level='all';
      else
        level='custom';
      end
    end  % method

    
    % ---------------------------------------------------------------------
    function wfAmount=getWFAmountForPFCategory(self,pfCategoryName)
      % Get the current amount for the given per-frame feature category.  Returns
      % one of 'normal', 'more', 'less', 'custom', and 'n/a' (if there are
      % no PFs in the category.
      pfNamesInCategory=self.getPFNamesInCategory(pfCategoryName);
      nPFsInCategory=length(pfNamesInCategory);
      if (nPFsInCategory==0)
        wfAmount='n/a';
        return
      end
      % Get the amount for the first PF
      pfName=pfNamesInCategory{1};
      wfAmountPutative=self.getWFAmountForPF(pfName);
      % If the first one is 'custom', no need to check the rest
      if isequal(wfAmountPutative,'custom')
        wfAmount='custom';
        return
      end
      % If we get here, wfAmountThisPutative is something other than
      % 'custom'
      % The rest have to match wfAmountPutative, or else the amount is
      % 'custom'
      for i = 2:nPFsInCategory
        pfName=pfNamesInCategory{i};
        wfAmountThis=self.getWFAmountForPF(pfName);
        if ~isequal(wfAmountThis,wfAmountPutative)
          wfAmount='custom';
          return
        end
      end
      % If we get here, they all match wfAmountPutative
      wfAmount=wfAmountPutative;
    end
            
      
    % ---------------------------------------------------------------------
    function wfAmount=getWFAmountForPF(self,pfIndex)
      % get the current amount for the given per-frame feature.  Returns
      % one of 'normal', 'more', 'less', and 'custom'.
      if ischar(pfIndex)
        % means it's really a name
        pfName=pfIndex;
        pfIndex=find(strcmp(pfName,self.subdialectPFNames));  
      else
        % pfIndex is really an index
        pfName=self.subdialectPFNames{pfIndex};
      end
      wfAmount='custom';  % if no match, we default to custom
      defaultPFTransTypes=self.defaultPFTransTypesFromName.(pfName);
      pfParams=self.vocabulary{pfIndex};  % all the WF params for the given PF
      nAmounts=length(self.wfAmounts);
      for j=1:nAmounts
        wfAmountTest=self.wfAmounts{j};
        isMatch=true;
        pfParamsTemplate=self.wfParamsFromAmount.(wfAmountTest);  
          % all the WF params for a single PF, from the pre-set amount
        for wfTypeIndex = 1:length(self.wfTypes)
          wfType = self.wfTypes{wfTypeIndex};
          if isfield(pfParams,wfType) && ~isfield(pfParamsTemplate,wfType),
            % If this window-feature type is present in the vocabulary, but not
            % in the template, then skip to next wfType
            continue
          end
          % check that the enablement in the vocab matches the template
          if pfParams.(wfType).enabled ~= pfParamsTemplate.(wfType).enabled ,
            isMatch=false;
            break
          end
          % For all the regular WF parameter names, check that they're
          % equal
          for i = 1:numel(self.wfParamNames),
            wfParamName = self.wfParamNames{i};
            if isequal(wfParamName,'trans_types')
              % Skip transformation types, we do those below
              continue
            end
            if pfParams.(wfType).values.(wfParamName) ~= ...
               pfParamsTemplate.(wfType).values.(wfParamName) ,
              isMatch=false;
              break
            end
          end
          if ~isMatch, break, end
          % If there is an extra WF parameter for this WF type, copy it over
          % also.
          extraWFParamName = self.wfExtraParamNames{wfTypeIndex};
          if ~isempty(extraWFParamName) && ...
             isfield(pfParamsTemplate.(wfType).values,extraWFParamName) && ...
             ~isequal(pfParams.(wfType).values.(extraWFParamName), ...
                      pfParamsTemplate.(wfType).values.(extraWFParamName)) ,
            isMatch=false;
            break
          end
          % Check the transformation types also
          if ~isequal(unique(pfParams.(wfType).values.trans_types) , ...
                      unique(defaultPFTransTypes) ) ,
            isMatch=false;
            break
          end
        end  % for wfTypeIndex = 1:length(self.wfTypes)
        if isMatch,
          wfAmount=wfAmountTest;
          break  % if match, we're done
        else
          continue  % if no match, try next wfAmount
        end
      end  % for wfAmountTest=self.wfAmounts
    end  % method
      
  end  % public instance methods
  
  
  % -----------------------------------------------------------------------
  % private instance methods
  methods (Access=private)
    % ---------------------------------------------------------------------
    function setWindowFeaturesOfSingleType(self,pfNdx,wfTypeNdx,wfParamTemplate)
      % For the per-frame feature with index pfNdx, sets the window features 
      % with type given by index wfTypeNdx to the amount given in wfParamTemplate.  

      wfType = self.wfTypes{wfTypeNdx};
      %wfParamTemplate=self.wfParamsFromAmount.(wfAmount);
      % Want to do something like:
      %   self.vocabulary{pfNdx}.(wfType)=pfParamTemplate.(wfType)  ,
      % but adapted to the particular PF, I guess ---ALT, Mar 31, 2013
      % something to copy from?
      if isfield(self.vocabulary{pfNdx},wfType) && ~isfield(wfParamTemplate,wfType),
        % If this window-feature type is present in the vocabulary, but not
        % in the template, then disable it in the vocab, and return.
        % (If this ever happens, doesn't it mean the template is broken?
        % And since we control the template, this need never happen.
        % --ALT, Mar 31, 2013)
        self.vocabulary{pfNdx}.(wfType).enabled = false;
        return
      end
      % Set the enablement in the vocab to match the template
      self.vocabulary{pfNdx}.(wfType).enabled = wfParamTemplate.(wfType).enabled;
      % For all the regular WF parameter names, copy them over
      for i = 1:numel(self.wfParamNames),
        wfParamName = self.wfParamNames{i};
        self.vocabulary{pfNdx}.(wfType).values.(wfParamName) = ...
          wfParamTemplate.(wfType).values.(wfParamName);
      end
      % If there is an extra WF parameter for this WF type, copy it over
      % also.
      extraWFParamName = self.wfExtraParamNames{wfTypeNdx};
      if ~isempty(extraWFParamName) 
        if isfield(wfParamTemplate.(wfType).values,extraWFParamName) ,
          value=wfParamTemplate.(wfType).values.(extraWFParamName);
        else
          % If no value is specified, set to empty matrix
          value=[];
        end
        self.vocabulary{pfNdx}.(wfType).values.(extraWFParamName) = ...
          value;
      end
    end  % method    
    
  end  % private methods
end  % classdef
