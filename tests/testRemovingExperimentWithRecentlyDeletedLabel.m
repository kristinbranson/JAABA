function success=testRemovingExperimentWithRecentlyDeletedLabel()

jabFileName='/groups/branson/bransonlab/adam/may01ing.jab';
gtMode=false;
data=JLabelData('setstatusfn',@(str)(fprintf('%s\n',str)), ...
                'clearstatusfn',@()(nop()));
data.openJabFile(jabFileName,gtMode);

% Add a new exp dir
nameOfExpDirToAdd='/groups/branson/bransonlab/adam/GMR_71G01_AE_01_TrpA_Rig2Plate14BowlD_20110707T154929';
data.AddExpDir(nameOfExpDirToAdd);

% Change the current target to the first fly of the new experiment
iExp=whichstr(nameOfExpDirToAdd,data.expdirs);
if isempty(iExp)
  error('testRemovingExperimentWithRecentlyDeletedLabel:addedExperimentIsMissing', ...
        'Added experiment seems to be missing.');
elseif length(iExp)>1
  error('testRemovingExperimentWithRecentlyDeletedLabel:addedExperimentPresentMoreThanOnce', ...
        'Added experiment seems to be present more than once in the JLabelData object.');
end  
iFly=1;
data.setCurrentTarget(iExp,iFly);

% Add a label
ts=(5050:5100);  % frame indices to set
iBehavior=whichstr('may01ing',data.labelnames);
isImportant=true;
data.SetLabel(iExp,iFly,ts,iBehavior,isImportant);

% Delete the just-created label
ts=(5040:5110);  % frame indices to set
iBehavior=0;  % by convention, this means to clear the labels from those frames
isImportant=false;  % I don't think this is used for iBehavior==0
data.SetLabel(iExp,iFly,ts,iBehavior,isImportant);

% Remove experiments with no labels
data.removeExperimentsWithNoLabels();

% Clean up
data.closeJabFile();
data=[];  %#ok
success=true;

end
