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
      self.scoreNorm = 0;  % default time stamp, number of days since Jan 0, 0000 (typically non-integer)
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
      if ~isprop(self,'selFeatures') || isempty(self.selFeatures),
        tfmodified = true;
        for i = 1:numel(self)
          self(i).selFeatures = SelFeatures.createEmpty();
        end
      else
        [self.selFeatures,sfmod] = SelFeatures.modernize(self.selFeatures);
        if sfmod,
          tmodified = true; %#ok<NASGU>
        end
      end
      
      if ~isfield(self.trainingParams,'nselfeatures'),
        tfmodified = true;
        for i = 1:numel(self)
          self(i).trainingParams.nselfeatures = 1000;
        end
      end
      
      if ~isfield(self,'predictOnlyCurrentFly'),
        tfmodified = true;
        for i = 1:numel(self)
          self(i).predictOnlyCurrentFly = false;
        end
      end

    end
        
  end 
  
end
