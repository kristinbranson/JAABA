classdef GrandlyUnifyController < handle
  properties
    model  % a ref to the model
    view  % a ref to the view
    workingDirName
  end
  methods
    function self=GrandlyUnifyController()
      self.model=GrandlyUnifyModel();
      self.view=GrandlyUnifyView(self.model,self);
      self.workingDirName=pwd();
    end  % constructor method

    % ---------------------------------------------------------------------
    function chooseProjectFileButtonPressed(self)
      [fileNameRel,pathName] = ...
        uigetfile({'*.mat','Matlab MAT files (*.mat)'}, ...
                  'Choose JAABA Project File', ...
                  self.workingDirName);
      if fileNameRel == 0; return; end;
      fileNameAbs = fullfile(pathName,fileNameRel);
      try
        self.model.setProjectFileName(fileNameAbs);
      catch excp
        uiwait(errordlg(sprintf('There was a problem setting the project file: %s',excp.message), ...
                        'Error', ...
                        'modal'));
      end
      self.workingDirName=fileparts(fileNameAbs);
      self.view.update();
    end
    
    % ---------------------------------------------------------------------
    function chooseScoreFeatureJabFileButtonPressed(self,scoreFeatureMatFileName)
      [fileNameRel,pathName] = ...
        uigetfile({'*.jab','JAABA files (*.jab)'}, ...
                  'Choose JAABA File', ...
                  self.workingDirName);
      if fileNameRel == 0; return; end;
      fileNameAbs = fullfile(pathName,fileNameRel);
      self.model.setScoreFeatureJabFileName(scoreFeatureMatFileName,fileNameAbs);
      self.workingDirName=fileparts(fileNameAbs);
      self.view.update();
    end
    
    % ---------------------------------------------------------------------
    function chooseClassifierFileButtonPressed(self)
      [fileNameRel,pathName] = ...
        uigetfile({'*.mat','Matlab MAT files (*.mat)'}, ...
                  'Choose JAABA Classifier File', ...
                  self.workingDirName);
      if fileNameRel == 0; return; end;
      fileNameAbs = fullfile(pathName,fileNameRel);
      self.model.setClassifierFileName(fileNameAbs);
      self.workingDirName=fileparts(fileNameAbs);
      self.view.update();
    end
    
    % ---------------------------------------------------------------------
    function gtExpDirsAddButtonPressed(self)
      dirNameAbs = ...
        uigetdir(self.workingDirName, ...
                 'Add Ground-Truth Experiment Directory');
      if dirNameAbs == 0; return; end;
      self.model.addGTExpDirName(dirNameAbs);
      self.workingDirName=fileparts(dirNameAbs);
      self.view.update();
    end
    
    % ---------------------------------------------------------------------
    function gtExpDirsRemoveButtonPressed(self)
      self.model.removeCurrentGTExpDirName();
      self.view.update();      
    end
    
    % ---------------------------------------------------------------------
    function convertButtonPressed(self)
      [~,baseName]=fileparts(self.model.classifierFileName);
      suggestedFileNameAbs=fullfile(self.workingDirName,[baseName '.jab']);
      [fileNameRel, pathName] = ...
        uiputfile({'*.jab','JAABA files (*.jab)'}, ...
                  'Save JAABA File As...', ...
                  suggestedFileNameAbs);
      if fileNameRel == 0; return; end;  % user hit cancel
      fileNameAbs = fullfile(pathName,fileNameRel);
      %projectFileName=self.model.projectFileName;
      %classifierFileName=self.model.classifierFileName;
      %gtExpDirNames=self.model.gtExpDirNames;
      self.view.spin();  
      try
        self.model.convert(fileNameAbs)
%         everythingFileFromOldStyleProjectAndClassifierFiles(...
%           fileNameAbs, ...
%           projectFileName, ...
%           classifierFileName, ...
%           gtExpDirNames);
      catch excp
        self.view.unspin();  
        uiwait(errordlg(excp.message,'Error','modal'));
        return
      end  % try/catch
      h=helpdlg('Conversion was successful','Success!');
      set(h,'windowstyle','modal');
      self.view.unspin();
      uiwait(h);
      self.quit();
    end  % method

    % ---------------------------------------------------------------------
    function quit(self)
      self.view.close();
      self.view=[];
      self.model=[];
    end
  end
end
