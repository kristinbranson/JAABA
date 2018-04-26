classdef ClassifierStuff < handle
  % Class that includes the classifier plus classifier metadata, like the
  % classifier training parameters
  
  properties
    type
    params
    timeStamp
    confThresholds % 1x2 vec
    scoreNorm % AL does windowdata.scoreNorm dup .scoreNorm?
    postProcessParams
    trainingParams
    featureNames
    windowdata 
    savewindowdata
    predictOnlyCurrentFly
    selFeatures
  end
  
  methods
    
    function self = ClassifierStuff(varargin)
      
      self.type = 'boosting';
      self.params = struct('dim',cell(0,1), ...
                           'error',[], ...
                           'dir',[], ....
                           'tr',[], ...
                           'alpha',[]);
      self.timeStamp = 0;  % default time stamp (bad idea?)
      self.confThresholds = [0 0];
      self.scoreNorm = nan;  % default time stamp, number of days since Jan 0, 0000 (typically non-integer)
      self.postProcessParams.method = 'Hysteresis';
      self.postProcessParams.hystopts(1) = struct('name','High Threshold','tag','hthres','value',0);
      self.postProcessParams.hystopts(2) = struct('name','Low Threshold','tag','lthres','value',0);
      self.postProcessParams.filtopts(1) = struct('name','Size','tag','size','value',1);
      self.postProcessParams.blen = 1;
      self.trainingParams = ...
        struct('iter',100, ...
               'iter_updates',10, ...
               'numSample',2500, ...
               'numBins',30, ...
               'CVfolds',7, ...
               'baseClassifierTypes','Decision Stumps', ...
               'baseClassifierSelected',1,...
               'nselfeatures',1000);
      self.featureNames = {};
      self.windowdata = WindowData.windowdata(0);
      self.savewindowdata = false;
      self.predictOnlyCurrentFly = false;
      self.selFeatures =  SelFeatures.createEmpty();
      
      % varargin should be key-value pairs
      % Keys should be property names, values should be initialization
      % values for the corresponding property
      keys = varargin(1:2:end);
      values = varargin(2:2:end);
      
      % Process the keys, assigning to any that match properties, and
      % throwing an error if a key matches no property
      propertyNames = properties(self);
      for iKey=1:length(keys)
        key = keys{iKey};
        iProperty = whichstr(key,propertyNames);
        if isempty(iProperty)
          error('Classifier:noPropertyMatchingKey', ...
                'Classifier: No property matching key %s',key);
        else
          propertyName = propertyNames{iProperty};
          value = values{iKey};
          self.(propertyName) = value;
        end
      end  
    end
        
    function tfmodified = modernize(self)
      tfmodified = false;
      for i = 1:numel(self)
        if isempty(self(i).savewindowdata)
          self(i).savewindowdata = false;
          tfmodified = true;
        end
      end
      
      % MK: May 09 2016
      for i = 1:numel(self)
        % AL: prob don't need to check isprop(...,'selFeatures') b/c the
        % property is in the class now -- any/all ClassifierStuffs running
        % with this version should get the prop 
        if ~isprop(self(i),'selFeatures') || isempty(self(i).selFeatures),
          tfmodified = true;
          self(i).selFeatures = SelFeatures.createEmpty();
        else
          [self(i).selFeatures,sfmod] = SelFeatures.modernize(self(i).selFeatures);
          tfmodified = tfmodified || sfmod;
        end
      end
      
      % MK: May 09 2016
      for i = 1:numel(self)
        if ~isfield(self(i).trainingParams,'nselfeatures'),
          tfmodified = true;
          self(i).trainingParams.nselfeatures = 1000;
        end
      end
    
      % MK: May 09 2016
      for i = 1:numel(self)
        % AL: isprop etc
        if ~isprop(self(i),'predictOnlyCurrentFly') || isempty(self(i).predictOnlyCurrentFly),
          tfmodified = true;
          self(i).predictOnlyCurrentFly = false;
        end
      end
      
%       % MK: modernize window data.
%       wdmodified = false;
%       for i = 1:numel(self)
%           [self(i).windowdata,wdm] =  WindowData.modernize(self(i).windowdata);
%           if wdm,
%               wdmodified = true;
%           end
%       end
%       if wdmodified,
%           tfmodified = true;
%       end
      
      % MK: modernize window data.
      % AL: fixed for classifiers without windowdata
      wdmodified = false;
      for i = 1:numel(self)
          if ~isequal(self(i).windowdata,[])
              [self(i).windowdata,wdm] =  WindowData.modernize(self(i).windowdata);
              if wdm,
                  wdmodified = true;
              end
          end
      end
      if wdmodified,
          tfmodified = true;
      end
      
    end
    
  end 
  
end
