classdef ChangeScoresAsFeaturesDialog < handle
  properties
    figureJLabel
    fileNameList={}  
      % The list of .jab files with classifiers whose scores will be used 
      % as per-frame features.
      % this will store _absolute_ file name
    timeStampList=[]
      % the time stamp of the classifier in fileNameList{i}
    scoreBaseNameList={}
      % the score file name (local, without the .mat) for the classifier in
      % fileNameList{i}
    iCurrent=[]
      % The index of the currently selected file in fileNameList.  Empty
      % iff fileNameList is empty.
    graphicsObject
    listBoxLabelText
    listBox
    addButton
    removeButton
    cancelButton
    doneButton
  end
  methods
    function self=ChangeScoresAsFeaturesDialog(fileNameList, ...
                                               timeStampList, ...
                                               scoreBaseNameList, ...
                                               figureJLabel)
      % need to keep this around to tell it when we're done
      self.figureJLabel=figureJLabel;
      self.fileNameList=fileNameList;
      self.timeStampList=timeStampList;
      self.scoreBaseNameList=scoreBaseNameList;
      self.iCurrent= ...
        fif(isempty(fileNameList),[],1);
      
      % Set figure layout parameters
      figureWidth=570;
      figureHeight=200;

      sideButtonWidth=80;  % add and remove are the side buttons
      sideButtonHeight=26;
      interSideButtonHeight=15;
      
      bottomButtonWidth=100;  % done and cancel are the bottom buttons
      buttonButtonHeight=30;
      interBottomButtonWidth=20;
      
      listBoxWidth=290;
      listBoxHeight=100;
      listBoxXOffset=158;
      listBoxYOffsetFromTop=20;
      listBoxYOffset=figureHeight-listBoxYOffsetFromTop-listBoxHeight;

      listBoxLabelTextXOffset=15;
      
      bottomButtonYOffset=18;
      bottomButtonRowXOffsetFromRight=20;
      
      % Figure the placement of the figure, want it centered relative to
      % figureJLabel
      jLabelPosition=get(figureJLabel,'position');
      figureXOffset=jLabelPosition(1)+(jLabelPosition(3)-figureWidth )/2;
      figureYOffset=jLabelPosition(2)+(jLabelPosition(4)-figureHeight)/2;
      
      % Create the figure
      self.graphicsObject=figure('tag','changeFeatureLexiconFigure', ...
                                 'position',[figureXOffset figureYOffset figureWidth figureHeight], ...
                                 'resize','off', ...
                                 'numbertitle','off', ...
                                 'menubar','none', ...
                                 'closeRequestFcn',@(g,e)(self.cancelButtonPressed()), ...
                                 'name','Change Target Type...', ...
                                 'windowstyle','normal');
      
      % make the popupmenu label, but make it invisible, so we can get its extent
      self.listBoxLabelText=uicontrol('parent',self.graphicsObject, ...
                                      'style','text', ...
                                      'tag','listBoxLabelText', ...
                                      'backgroundcolor',get(self.graphicsObject,'color'), ...
                                      'string',{'List of files', ...
                                                'whose output scores', ...
                                                'will be used as', ...
                                                'per-frame features'}', ...
                                      'HorizontalAlignment','right',...
                                      'visible','on', ...
                                      'fontsize',12);
      extent=get(self.listBoxLabelText,'extent');
      extentSize=extent(3:4);
      listBoxLabelTextSize=extentSize+4;  % pad a little
      
      % align the list box label with the list box in y
      listBoxLabelTextYOffset=listBoxYOffset+(listBoxHeight-listBoxLabelTextSize(2))/2;
      
      % done button is right-justified, with a margin of bottomButtonRowXOffsetFromRight
      doneButtonXOffset=figureWidth-bottomButtonRowXOffsetFromRight-bottomButtonWidth;
      
      % cancel button is to the left off the done button, set off by interBottomButtonWidth
      cancelButtonXOffset=doneButtonXOffset-interBottomButtonWidth-bottomButtonWidth;
            
      % add and remove buttons are aligned in X, and have space of
      % interSideButtonHeight between them in Y
      % the column of two buttons is aligned in Y with the listbox, and 
      % each button is centered in X in the space between the right side of
      % the listbox and the right edge of the figure
      allSideButtonsHeight=sideButtonHeight+interSideButtonHeight+sideButtonHeight;
      removeButtonYOffset=listBoxYOffset+(listBoxHeight-allSideButtonsHeight)/2;
      addButtonYOffset=removeButtonYOffset+sideButtonHeight+interSideButtonHeight;
      widthOfSpaceToTheRightOfListBox=figureWidth-(listBoxXOffset+listBoxWidth);
      removeButtonXOffset=listBoxXOffset+listBoxWidth+(widthOfSpaceToTheRightOfListBox-sideButtonWidth)/2;
      addButtonXOffset=removeButtonXOffset;
      
      % Create and/or position the controls
      set(self.listBoxLabelText, ...
          'position',[listBoxLabelTextXOffset listBoxLabelTextYOffset listBoxLabelTextSize], ...
          'visible','on');
      self.listBox=uicontrol('parent',self.graphicsObject, ...
                             'style','listbox', ...
                             'fontsize',12, ...
                             'position',[listBoxXOffset listBoxYOffset listBoxWidth listBoxHeight]);
      self.addButton=uicontrol('parent',self.graphicsObject, ...
                               'style','pushbutton', ...
                               'string','Add', ...
                               'fontsize',12, ...
                               'position',[addButtonXOffset addButtonYOffset sideButtonWidth sideButtonHeight], ...
                               'callback',@(g,e)(self.addButtonPressed()));
      self.removeButton=uicontrol('parent',self.graphicsObject, ...
                                'style','pushbutton', ...
                                'string','Remove', ...
                                'fontsize',12, ...
                                'position',[removeButtonXOffset removeButtonYOffset sideButtonWidth sideButtonHeight], ...
                                'callback',@(g,e)(self.removeButtonPressed()));
      self.cancelButton=uicontrol('parent',self.graphicsObject, ...
                                  'style','pushbutton', ...
                                  'string','Cancel', ...
                                  'fontsize',12, ...
                                  'position',[cancelButtonXOffset bottomButtonYOffset bottomButtonWidth buttonButtonHeight], ...
                                  'callback',@(g,e)(self.cancelButtonPressed()));
      self.doneButton=uicontrol('parent',self.graphicsObject, ...
                                'style','pushbutton', ...
                                'string','Done', ...
                                'fontsize',12, ...
                                'position',[doneButtonXOffset bottomButtonYOffset bottomButtonWidth buttonButtonHeight], ...
                                'callback',@(g,e)(self.doneButtonPressed()));
                              
      % sync the view with the model
      self.updateView();
    end

    function addButtonPressed(self)
      [fileNameRel,pathName] = ...
        uigetfile({'*.jab','JAABA Everything Files (*.jab)'}, ...
                  'Add .jab file containing classifier to be used as input');
      if fileNameRel == 0; return; end;
      fileNameAbs = fullfile(pathName,fileNameRel);
      everythingParams = load(fileNameAbs,'-mat');
      if isempty(everythingParams.classifier)
        uiwait(errordlg(sprintf('%s does not contain a classifier.',fileNameRel), ...
                        'Error', ...
                        'modal'));
        return
      end
      % Check that the classifier has a time stamp
      classifier=everythingParams.classifier;
      if isfield(classifier,'timeStamp');
        timeStamp = classifier.timeStamp;
      else
        uiwait(errordlg('The classifier in the selected file lacks a timestamp.', ...
                        'Error', ...
                        'modal'));
        return
      end
      % Add the name of the score file (without the .mat extension)
      if isfield(classifier,'file') && isfield(classifier.file,'scorefilename')
        scoreFileName = classifier.scorefilename;
        [~,scoreBaseName] = fileparts(scoreFileName);
      elseif isfield(everythingParams,'behaviors') && ...
             isfield(everythingParams.behaviors,'names') && ...
             ~isempty(everythingParams.behaviors.names)
        behaviorName=everythingParams.behaviors.names{1};
        scoreBaseName = sprintf('scores_%s',behaviorName);
      else
        uiwait(errordlg('Unable to determine score file name for classifier.', ...
                        'Error', ...
                        'modal'));
        return
      end
      self.fileNameList{end+1}=fileNameAbs;
      self.timeStampList(end+1)=timeStamp;
      self.scoreBaseNameList{end+1}=scoreBaseName;
      self.iCurrent = length(self.fileNameList);
      self.updateView();
    end
    
    function removeButtonPressed(self)
      i = self.iCurrent;
      if isempty(i), return; end
      nScoresAsInput=length(self.fileNameList);
      % delete the i'th entries
      self.fileNameList(i) = [];
      self.timeStampList(i) = [];
      self.scoreBaseNameList(i) = [];
      % update i if we just deleted the n'th element
      if (i==nScoresAsInput)
        i=fif(i>1,i-1,[]);
      end
      self.iCurrent=i;
      % Update the "view" to reflect the changed "model"
      self.updateView();      
    end
    
    function cancelButtonPressed(self)
      delete(self.graphicsObject);
    end

    function doneButtonPressed(self)
      % Delete the figure
      delete(self.graphicsObject);
      % Call the appropriate function to notify the JLabel "object" that 
      % we're done.
      JLabel('changeScoresAsFeaturesDone', ...
             self.figureJLabel, ...
             self.fileNameList, ...
             self.timeStampList, ...
             self.scoreBaseNameList);
    end
    
    function updateView(self)
      % Update the "view" to reflect the current "model"
      
      % Update the list of scores-an-inputs
      set(self.listBox,'String',self.fileNameList);
      set(self.listBox,'Value',self.iCurrent);

      % Disble the Remove button iff the list of scores-as-inputs is empty
      set(self.removeButton,'enable',offIff(isempty(self.fileNameList)));
    end  % method
    
  end
end
