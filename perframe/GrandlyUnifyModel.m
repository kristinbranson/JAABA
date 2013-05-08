classdef GrandlyUnifyModel < handle
  properties (GetAccess=public, SetAccess=private)
    projectFileName=''
      % The absolute path of the project file name to be converted.
    classifierFileName=''  
      % The absolute path of the classifier file name to be converted.
    gtExpDirNames={}
      % The list of ground-truth experiment directories.  Absolute paths,
      % again.
    iCurrentGTExpDir=[]
      % The index of the currently selected GT experiment dir.  Empty
      % iff gtExpDirNamesAbs is empty.
    workingDirName='';
      % the absolute path of the default dir used for file choosers, etc.
  end
  methods
    function self=GrandlyUnifyModel()
      self.workingDirName=pwd();
    end  % constructor method
    function setProjectFileName(self,projectFileName)
      projectFileNameAbs=absolutifyFileName(projectFileName,self.workingDirName);
      self.projectFileName=projectFileNameAbs;
      self.workingDirName=fileparts(projectFileNameAbs);
    end
    function setClassifierFileName(self,fileName)
      fileNameAbs=absolutifyFileName(fileName,self.workingDirName);
      self.classifierFileName=fileNameAbs;
      self.workingDirName=fileparts(fileNameAbs);
    end
    function result=isProjectFileSpecified(self)
      result=~isempty(self.projectFileName);
    end
    function result=isClassifierFileSpecified(self)
      result=~isempty(self.classifierFileName);
    end
    function addGTExpDirName(self,dirName)
      dirNameAbs=absolutifyFileName(dirName,self.workingDirName);
      self.gtExpDirNames{end+1}=dirNameAbs;
      self.iCurrentGTExpDir = length(self.gtExpDirNames);
      self.workingDirName=fileparts(dirName);
    end
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
    function convert(self,jabFileName)
      jabFileNameAbs=absolutifyFileName(jabFileName,self.workingDirName);
      everythingFileFromOldStyleProjectAndClassifierFiles(...
        jabFileNameAbs, ...
        self.projectFileName, ...
        self.classifierFileName, ...
        self.gtExpDirNames);    
    end  % method
  end  % methods
end  % classdef
