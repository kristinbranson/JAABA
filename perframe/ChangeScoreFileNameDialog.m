classdef ChangeScoreFileNameDialog < handle
  properties
    figureJLabel
    scoreFileName
    editBoxEnterKeyPressed=false
    graphicsObject
    labelText
    editBox
    cancelButton
    doneButton
  end
  methods
    function self=ChangeScoreFileNameDialog(scoreFileName,figureJLabel)
      % need to keep this around to tell it when we're done
      self.scoreFileName=scoreFileName;
      self.figureJLabel=figureJLabel;
      
      % Set figure layout parameters
      figureWidth=350;
      figureHeight=130;
      editBoxHeight=24;
      editBoxWidth=120;
      interLabelEditBoxWidth=13;
      buttonWidth=100;
      buttonHeight=30;
      interButtonWidth=40;
      editBoxYOffset=82;
      buttonYBaseline=22;
      
      % Figure the placement of the figure, want it centered relative to
      % figureJLabel
      jLabelPosition=get(figureJLabel,'position');
      figureXOffset=jLabelPosition(1)+(jLabelPosition(3)-figureWidth )/2;
      figureYOffset=jLabelPosition(2)+(jLabelPosition(4)-figureHeight)/2;
      
      % Create the figure
      self.graphicsObject=figure('position',[figureXOffset figureYOffset figureWidth figureHeight], ...
                                 'resize','off', ...
                                 'numbertitle','off', ...
                                 'menubar','none', ...
                                 'closeRequestFcn',@(g,e)(self.cancelButtonPressed()), ...
                                 'name','Change Score File Name...', ...
                                 'windowstyle','modal');
      
      % make the editBox label, but make it invisible, so we can get its extent
      self.labelText=uicontrol('parent',self.graphicsObject, ...
                               'style','text', ...
                               'backgroundcolor',get(self.graphicsObject,'color'), ...
                               'string','Score File Name:', ...
                               'HorizontalAlignment','right',...
                               'visible','off', ...
                               'fontsize',12);
      extent=get(self.labelText,'extent');
      extentSize=extent(3:4);
      labelTextSize=extentSize+2;  % pad a little
      
      % calculate placement of editBox and label
      editBoxAndLabelWidth=labelTextSize(1)+interLabelEditBoxWidth+editBoxWidth;
      labelTextXOffset=(figureWidth-editBoxAndLabelWidth)/2;
      editBoxXOffset=labelTextXOffset+labelTextSize(1)+interLabelEditBoxWidth;
      labelTextYOffset=editBoxYOffset+(editBoxHeight-labelTextSize(2))/2;
      
      % calculate placement of buttons
      allButtonsWidth=buttonWidth+interButtonWidth+buttonWidth;
      cancelButtonXOffset=(figureWidth-allButtonsWidth)/2;
      doneButtonXOffset=cancelButtonXOffset+buttonWidth+interButtonWidth;
            
      % Create and/or position the controls
      set(self.labelText, ...
          'position',[labelTextXOffset labelTextYOffset labelTextSize], ...
          'visible','on');
      self.editBox=uicontrol('parent',self.graphicsObject, ...
                             'style','edit', ...
                             'String',scoreFileName, ...
                             'fontsize',12, ...
                             'callback',@(g,e)(self.editBoxEdited()), ...
                             'keypressfcn',@(g,e)(self.editBoxKeyPressed(e)), ...
                             'position',[editBoxXOffset editBoxYOffset editBoxWidth editBoxHeight]);
      self.cancelButton=uicontrol('parent',self.graphicsObject, ...
                                  'style','pushbutton', ...
                                  'string','Cancel', ...
                                  'fontsize',12, ...
                                  'position',[cancelButtonXOffset buttonYBaseline buttonWidth buttonHeight], ...
                                  'callback',@(g,e)(self.cancelButtonPressed()));
      self.doneButton=uicontrol('parent',self.graphicsObject, ...
                                'style','pushbutton', ...
                                'string','Done', ...
                                'fontsize',12, ...
                                'position',[doneButtonXOffset buttonYBaseline buttonWidth buttonHeight], ...
                                'callback',@(g,e)(self.doneButtonPressed()));
    end

    function editBoxKeyPressed(self,event)
      if isequal(event.Key,'return')
        self.editBoxEnterKeyPressed=true;
        % editBoxEdited() will get called next, which checks this and exits
        % if it was pressed and the new behavior name is valid.
      end
    end
    
    function editBoxEdited(self)
      % This gets called when user hits the enter key after clicking in the
      % edit box, or when they click another control having clicked in the
      % edit box, or if they click off the editBox having edited it.
      %fprintf('editBoxEdited() entry.\n');
      scoreFileNameNew=get(self.editBox,'string');
      if JLabelData.isValidScoreFileName(scoreFileNameNew), ...
        self.scoreFileName=scoreFileNameNew;
        if self.editBoxEnterKeyPressed,
          self.closeAndSendMessageToJLabel();
        end
      else
        set(self.editBox,'string',self.scoreFileName);
      end
      self.editBoxEnterKeyPressed=false;  % clear for future use
      %fprintf('editBoxEdited() exit.\n');
    end
    
    function cancelButtonPressed(self)
      delete(self.graphicsObject);
    end

    function doneButtonPressed(self)
      %fprintf('doneButtonPressed() entry.\n');
      % It should never happen that we get here with the editBox string
      % and self.scoreFileName out-of-sync, so there's not much to do here.
      self.closeAndSendMessageToJLabel();
      %fprintf('doneButtonPressed() exit.\n');
    end
    
    function closeAndSendMessageToJLabel(self)
      % Called when done button clicked or enter pressed in the editBox, 
      % to close this figure and signal back to JLabel.
      %fprintf('closeAndSendMessageToJLabel() entry.\n');
      % Delete the figure
      delete(self.graphicsObject);
      % Call the appropriate function to notify the JLabel "object" that 
      % we're done.
      JLabel('changeScoreFileNameDone', ...
             self.figureJLabel, ...
             self.scoreFileName);
      %fprintf('closeAndSendMessageToJLabel() exit.\n');
    end
    
  end
end
