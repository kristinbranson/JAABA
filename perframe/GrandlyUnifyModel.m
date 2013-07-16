classdef GrandlyUnifyModel < handle
  properties (GetAccess=public, SetAccess=private)
    projectFileName=''
      % The project file name to be converted.
    scoreFeatureMatFileNames=cell(0,1);
      % The file names of any score feature classifiers in the project file
    scoreFeatureJabFileNames=cell(0,1);
      % The names of the jab files that go with each score feature
      % If the jab file name for a particular score feature is unspecified,
      % that element will be the empty string.
    classifierFileName=''  
      % The classifier file name to be converted.
    gtExpDirNames={}
      % The list of ground-truth experiment directory names.
    iCurrentGTExpDir=[]
      % The index of the currently selected GT experiment dir.  Empty
      % iff gtExpDirNames is empty.
  end
  
  methods
    % ---------------------------------------------------------------------
    function self=GrandlyUnifyModel()
    end  % constructor method

    % ---------------------------------------------------------------------
    function setProjectFileName(self,projectFileName)
      oldProjectFileName=self.projectFileName;
      oldScoreFeatureMatFileNames=self.scoreFeatureMatFileNames;
      oldScoreFeatureJabFileNames=self.scoreFeatureJabFileNames;
      try
        self.projectFileName=projectFileName;
        % Need to read project file, extract score feature names, if any.
        [self.scoreFeatureMatFileNames,self.scoreFeatureJabFileNames]= ...
          GrandlyUnifyModel.extractScoreFeatureNames(projectFileName);
      catch excp
        self.projectFileName=oldProjectFileName;
        self.scoreFeatureMatFileNames=oldScoreFeatureMatFileNames;
        self.scoreFeatureJabFileNames=oldScoreFeatureJabFileNames;
        rethrow(excp);
      end
    end

    % ---------------------------------------------------------------------
    function setScoreFeatureJabFileName(self,matFileName,jabFileName)
      i=whichstr(matFileName,self.scoreFeatureMatFileNames);
      if isempty(i)
        error('GrandlyUnifyModel:noSuchScoreFeatureMatFileName', ...
              sprintf('No score feature mat file name %s.',matFileName));  %#ok
      else
        self.scoreFeatureJabFileNames{i}=jabFileName;
      end
    end

    % ---------------------------------------------------------------------
    function setClassifierFileName(self,fileName)
      self.classifierFileName=fileName;
    end

    % ---------------------------------------------------------------------
    function result=isProjectFileSpecified(self)
      result=~isempty(self.projectFileName);
    end

    % ---------------------------------------------------------------------
    function result=isScoreFeatureJabFileSpecified(self,index)
      if ischar(index)
        matFileName=index;
        i=whichstr(matFileName,self.scoreFeatureMatFileNames);
        if isempty(i)
          error('GrandlyUnifyModel:noSuchScoreFeatureMatFileName', ...
                sprintf('No score feature mat file name %s.',matFileName));  %#ok
        end
      else
        i=index;
      end
      result=~isempty(self.scoreFeatureJabFileNames{i});
    end

    % ---------------------------------------------------------------------
    function result=isClassifierFileSpecified(self)
      result=~isempty(self.classifierFileName);
    end

    % ---------------------------------------------------------------------
    function addGTExpDirName(self,dirName)
      %dirNameAbs=absolutifyFileName(dirName,self.workingDirName);
      self.gtExpDirNames{end+1}=dirName;
      self.iCurrentGTExpDir = length(self.gtExpDirNames);
      %self.workingDirName=fileparts(dirName);
    end

    % ---------------------------------------------------------------------
    function removeCurrentGTExpDirName(self)
      i = self.iCurrentGTExpDir;
      if isempty(i), return; end
      n=length(self.gtExpDirNames);
      % delete the i'th entries
      self.gtExpDirNames(i) = [];
      % update i if we just deleted the n'th element
      if (i==n)
        i=fif(i>1,i-1,[]);
      end
      self.iCurrentGTExpDir=i;
    end  % method

    % ---------------------------------------------------------------------
    function convert(self,jabFileName)
      %jabFileNameAbs=absolutifyFileName(jabFileName,self.workingDirName);
      everythingFileFromOldStyleProjectAndClassifierFiles(...
        jabFileName, ...
        self.projectFileName, ...
        self.classifierFileName, ...
        self.gtExpDirNames, ...
        self.scoreFeatureMatFileNames, ...
        self.scoreFeatureJabFileNames);    
    end  % method
  end  % methods

  methods (Static)
    % ---------------------------------------------------------------------
    function [scoreFeatureMatFileNames,scoreFeatureJabFileNames]=extractScoreFeatureNames(projectFileName)
      projectParams=load(projectFileName,'-mat');
      scoresAsInput=projectParams.scoresinput;
      scoreFeatureMatFileNames={scoresAsInput.classifierfile}';
      nScoreFeatures=length(scoresAsInput);
      scoreFeatureJabFileNames=cell(nScoreFeatures,1);
      for i=1:nScoreFeatures
        scoreFeatureJabFileNames{i}='';
      end
    end % method
  end  % class methods
end  % classdef
