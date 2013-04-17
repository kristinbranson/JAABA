classdef EverythingParams < handle
  properties
    featureLexiconName
    behaviors
    file
    trxGraphicParams
    labelGraphicParams
    scoreFeatures
    sublexiconPFNames
    windowFeaturesParams
    labels
    gtLabels
    expDirNames
    gtExpDirNames
    classifier
    version='0.5.0'
  end  % properties
  methods (Access=private)
    function initFromFeatureLexiconName(self,featureLexiconName)
      self.featureLexiconName=featureLexiconName;
      self.behaviors.names = {};
      self.file.perframedir = 'perframe';
      self.file.clipsdir = 'clips';
      %self.windowfeatures = struct;
      self.behaviors.labelcolors = [0.7,0,0,0,0,0.7];
      self.behaviors.unknowncolor = [0,0,0];
      self.trxGraphicParams.colormap = 'jet';
      self.trxGraphicParams.colormap_multiplier = 0.7;
      self.trxGraphicParams.extra_linestyle = '-';
      self.trxGraphicParams.extra_marker = '.';
      self.trxGraphicParams.extra_markersize = 12;
      self.labelGraphicParams.colormap = 'line';
      self.labelGraphicParams.linewidth = 3;
      self.file.scorefilename = '';
      self.file.trxfilename = '';
      self.file.moviefilename = '';
      self.scoreFeatures = struct('classifierfile',{},'ts',{},'scorefilename',{});
      featureLexicon=featureLexiconFromFeatureLexiconName(featureLexiconName);
      featureLexiconPFNames = fieldnames(featureLexicon.perframe);
      self.sublexiconPFNames = featureLexiconPFNames;
      self.windowFeaturesParams=[];  % seems wrong---shouldn't it match sublexiconPFNames+scoreFeatures?
      self.labels=[];  % do I need to make a struct array with the right fields?
      self.gtLabels=[];  % do I need to make a struct array with the right fields?
      self.expDirNames={};
      self.gtExpDirNames={};
      self.classifier=[]; % do I need to make a struct with the right fields?
    end  % method      
    function initFromJLabelData(self,jld)
      self.featureLexiconName=jld.featureLexiconName;
      self.featureLexicon=jld.featureLexicon;
      self.scoreFeatures=jld.scoreFeatures;
      subdialectPFNames=jld.allperframefns;
      nScoreFeaturess=length(jld.scoreFeatures);
      sublexiconPFNames=subdialectPFNames(1:end-nScoreFeaturess);
      self.sublexiconPFNames=sublexiconPFNames;
      self.behaviorself.type=jld.targettype;
      self.behaviorself.names=jld.labelnames;
      self.behaviorself.labelcolors=jld.labelcolors;
      self.behaviors.unknowncolor=jld.unknowncolor;
      self.file.moviefilename=jld.moviefilename;
      self.file.trxfilename=jld.trxfilename;
      self.file.scorefilename=jld.scorefilename;
      self.file.clipsdir=jld.clipsdir;                  
      self.file.perframedir=jld.perframedir;                  
      self.labelGraphicParams=jld.labelGraphicParams;
      self.trxGraphicParams=jld.trxGraphicParams;
      self.landmarkParams=jld.landmark_params;
      % Get the labels, put them in self
      if jld.gtMode ,
        self.gtExpDirNames=jld.expdirs;
        self.expDirNames=jld.otherModeLabelsEtc.expDirNames;
      else
        self.expDirNames=jld.expdirs;
        self.gtExpDirNames=jld.otherModeLabelsEtc.expDirNames;
      end
      [self.labels,self.gtLabels]=jld.storeAndGetLabelsAndGTLabels();
      % Get the window feature params, put in self
      self.windowFeaturesParams=jld.windowfeaturesparams;
      % Put the classifier in self
      self.classifier=jld.getClassifier();
    end  % method
  end  % private methods
  methods
    function self=EverythingParams(varargin)
      if length(varargin)==1 && ischar(varargin{1})
        self.initFromFeatureLexiconName(varargin{1});
      elseif length(varargin)==1 && isequal(class(varargin{1}),'JLabelData')
        self.initFromJLabelData(varargin{1})
      else
        error('EverythingParams:badArgumentsToConstructor', ...
              'The arguments to EverythingParams() are no good');
      end
    end  % constructor method
  end
end  % classdef
