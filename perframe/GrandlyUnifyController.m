classdef GrandlyUnifyController < handle
  properties
    model  % a ref to the model
    view  % a ref to the view
  end
  methods
    function self=GrandlyUnifyController()
      self.model=GrandlyUnifyModel();
      self.view=GrandlyUnifyView(self.model,self);
    end  % constructor method

    % ---------------------------------------------------------------------
    function chooseProjectFileButtonPressed(self)
      [fileNameRel,pathName] = ...
        uigetfile({'*.mat','Matlab MAT files (*.mat)'}, ...
                  'Choose JAABA Project File');
      if fileNameRel == 0; return; end;
      fileNameAbs = fullfile(pathName,fileNameRel);
      self.model.projectFileName=fileNameAbs;
      self.view.update();
    end
    
    % ---------------------------------------------------------------------
    function chooseClassifierFileButtonPressed(self)
      [fileNameRel,pathName] = ...
        uigetfile({'*.mat','Matlab MAT files (*.mat)'}, ...
                  'Choose JAABA Classifier File');
      if fileNameRel == 0; return; end;
      fileNameAbs = fullfile(pathName,fileNameRel);
      self.model.classifierFileName=fileNameAbs;
      self.view.update();
    end
    
    % ---------------------------------------------------------------------
    function gtExpDirsAddButtonPressed(self)
      dirNameAbs = ...
        uigetdir('', ...
                 'Add Ground-Truth Experiment Directory');
      if dirNameAbs == 0; return; end;
      self.model.addGTExpDirName(dirNameAbs);
      self.view.update();
    end
    
    % ---------------------------------------------------------------------
    function gtExpDirsRemoveButtonPressed(self)
      self.model.removeCurrentGTExpDirName();
      self.view.update();      
    end
    
    % ---------------------------------------------------------------------
    function convertButtonPressed(self)
      [fileNameRel, pathName] = ...
        uiputfile({'*.jab','JAABA files (*.jab)'}, ...
                  'Save JAABA File As...');
      if fileNameRel == 0; return; end;  % user hit cancel
      self.view.spin();  
      fileNameAbs = fullfile(pathName,fileNameRel);
      projectFileName=self.model.projectFileName;
      classifierFileName=self.model.classifierFileName;
      gtExpDirNames=self.model.gtExpDirNames;
      try
        everythingFileFromOldStyleProjectAndClassifierFiles(...
          fileNameAbs, ...
          projectFileName, ...
          classifierFileName, ...
          gtExpDirNames);
      catch excp
        uiwait(errordlg(excp.message,'Error','modal'));
        return
      end  % try/catch
      h=helpdlg('Conversion was successful','Success!');
      set(h,'windowstyle','modal');
      self.view.unspin();  
      uiwait(h);
    end  % method
  end
end
