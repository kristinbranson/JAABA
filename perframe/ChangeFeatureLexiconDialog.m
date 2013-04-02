classdef ChangeFeatureLexiconDialog < handle
  properties
    figureJLabel
    featureLexiconNameList
    graphicsObject
    popupMenuLabelText
    featureLexiconPopupMenu
    cancelButton
    doneButton
  end
  methods
    function self=ChangeFeatureLexiconDialog(featureLexiconName,figureJLabel)
      % need to keep this around to tell it when we're done
      self.figureJLabel=figureJLabel;
      
      % Set figure layout parameters
      figureWidth=350;
      figureHeight=130;
      popupMenuHeight=24;
      popupMenuWidth=120;
      interLabelPopupMenuWidth=13;
      buttonWidth=100;
      buttonHeight=30;
      interButtonWidth=40;
      popupMenuYOffset=82;
      buttonYBaseline=22;
      
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
                                 'windowstyle','modal');
      
      % make the popupmenu label, but make it invisible, so we can get its extent
      self.popupMenuLabelText=uicontrol('parent',self.graphicsObject, ...
                                        'style','text', ...
                                        'tag','popupMenuLabelText', ...
                                        'backgroundcolor',get(self.graphicsObject,'color'), ...
                                        'string','Target Type:', ...
                                        'HorizontalAlignment','right',...
                                        'visible','off', ...
                                        'fontsize',12);
      extent=get(self.popupMenuLabelText,'extent');
      extentSize=extent(3:4);
      popupMenuLabelTextSize=extentSize+2;  % pad a little
      
      % calculate placement of popupmenu and label
      popupMenuAndLabelWidth=popupMenuLabelTextSize(1)+interLabelPopupMenuWidth+popupMenuWidth;
      popupMenuLabelTextXOffset=(figureWidth-popupMenuAndLabelWidth)/2;
      popupMenuXOffset=popupMenuLabelTextXOffset+popupMenuLabelTextSize(1)+interLabelPopupMenuWidth;
      popupMenuLabelTextYOffset=popupMenuYOffset+(popupMenuHeight-popupMenuLabelTextSize(2))/2;
      
      % calculate placement of buttons
      allButtonsWidth=buttonWidth+interButtonWidth+buttonWidth;
      cancelButtonXOffset=(figureWidth-allButtonsWidth)/2;
      doneButtonXOffset=cancelButtonXOffset+buttonWidth+interButtonWidth;
      
      % Get the list of possible feature lexicon names, and the index of
      % the current one
      standardFeatureLexiconNameList = getFeatureLexiconListsFromXML();
      if isequal(featureLexiconName,'custom')
        self.featureLexiconNameList={'custom';standardFeatureLexiconNameList};
        iFeatureLexiconName=1;
      else
        self.featureLexiconNameList = standardFeatureLexiconNameList;
        iFeatureLexiconName=find(strcmp(featureLexiconName,self.featureLexiconNameList));
      end
      
      % Create and/or position the controls
      set(self.popupMenuLabelText, ...
          'position',[popupMenuLabelTextXOffset popupMenuLabelTextYOffset popupMenuLabelTextSize], ...
          'visible','on');
      self.featureLexiconPopupMenu=uicontrol('parent',self.graphicsObject, ...
                                             'style','popupmenu', ...
                                             'tag','featureLexiconPopupMenu', ...
                                             'String',self.featureLexiconNameList, ...
                                             'Value',iFeatureLexiconName, ...
                                             'fontsize',12, ...
                                             'position',[popupMenuXOffset popupMenuYOffset popupMenuWidth popupMenuHeight]);
      self.cancelButton=uicontrol('parent',self.graphicsObject, ...
                                  'style','pushbutton', ...
                                  'tag','cancelButton', ...
                                  'string','Cancel', ...
                                  'fontsize',12, ...
                                  'position',[cancelButtonXOffset buttonYBaseline buttonWidth buttonHeight], ...
                                  'callback',@(g,e)(self.cancelButtonPressed()));
      self.doneButton=uicontrol('parent',self.graphicsObject, ...
                                'style','pushbutton', ...
                                'tag','doneButton', ...
                                'string','Done', ...
                                'fontsize',12, ...
                                'position',[doneButtonXOffset buttonYBaseline buttonWidth buttonHeight], ...
                                'callback',@(g,e)(self.doneButtonPressed()));
    end

    function cancelButtonPressed(self)
      delete(self.graphicsObject);
    end

    function doneButtonPressed(self)
      % Get the info we need out of the handles
      figureJLabel=self.figureJLabel;
      featureLexiconNameList=self.featureLexiconNameList;
      % get the currently selected feature lexicon name
      featureLexiconName=featureLexiconNameList{get(self.featureLexiconPopupMenu,'value')};
      % Delete the figure
      delete(self.graphicsObject);
      % Call the appropriate function to notify the JLabel "object" that 
      % we're done.
      JLabel('changeFeatureLexiconDone', ...
             figureJLabel, ...
             featureLexiconName);
    end
  end
end
