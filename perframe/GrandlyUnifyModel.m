classdef GrandlyUnifyModel < handle
  properties
    projectFileName=''
      % The absolute path of the project file name to be converted.
    classifierFileName=''  
      % The absolute path of the classifier file name to be converted.
    gtExpDirNames={}
      % The list of ground-truth experiment directories.  Absolute paths,
      % again.
    iCurrentGTExpDir=[]
      % The index of the currently selected GT experiment dir.  Empty
      % iff gtExpDirNames is empty.      
  end    
  properties (Dependent)
    projectFileSpecified
    classifierFileSpecified
  end
  methods
    function self=GrandlyUnifyModel()
    end  % constructor method
    function set.projectFileName(self,newValue)
      self.projectFileName=newValue;
    end
    function set.classifierFileName(self,newValue)
      self.classifierFileName=newValue;
    end
    function result=get.projectFileSpecified(self)
      result=~isempty(self.projectFileName);
    end
    function result=get.classifierFileSpecified(self)
      result=~isempty(self.classifierFileName);
    end
    function addGTExpDirName(self,fileNameAbs)
      self.gtExpDirNames{end+1}=fileNameAbs;
      self.iCurrentGTExpDir = length(self.gtExpDirNames);
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
    end
    
  end
end
