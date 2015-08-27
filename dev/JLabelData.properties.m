%% PROJECT-LEVEL CONFIG
targettype
featureLexiconName % such as 'flies', or 'vivek_mice'
featureLexicon
isST % Logical scalar, true for ST mode
usePastOnly % whether to only use past information when predicting the current frame

landmark_params
perframe_params % computing per-frame properties

filetypes % constant: files per experiment directory, {'movie','trx','label','gt_label','perframedir','clipsdir','scores'};
moviefilename
trxfilename
scorefilename % cellstr !IDXCLS
perframedir
clipsdir
scores
stfeatures

windowdatachunk_radius % shared by all classifiers.
predictwindowdatachunk_radius % shared by all classifiers.

%% EXPERIMENTS
expdirs % cell array of input experiment directory paths for the current GT mode
nflies_per_exp % array of number of flies in each experiment
firstframes_per_exp % firstframes_per_exp{expi}(fly)
endframes_per_exp % endframes_per_exp{expi}(fly)
frac_sex_per_exp % sex per experiment, fly
sex_per_exp
hassex % scalar, whether sex is computed
hasperframesex % scalar, whether sex is computed on a per-frame basis

fileexists % matrix of size nexps x numel(file_types), where fileexists(expi,filei)
filetimestamps % same structure as fileexists
allfilesexist % whether all necessary files for all experiments exist
filesfixable % whether we can generate any missing files
perframeGenerate % whether user has given permission to generate the perframe files
perframeOverwrite % to overwrite or keep the perframe files.
arenawarn % warn about removing arena features.
hasarenaparams

%% TARGET
expi % currently selected experiment
flies % currently selected flies
trx % last-used trajectories (one experiment, all flies)
t0_curr % first frame that all flies currently selected are tracked
t1_curr % last frame that all flies currently selected are tracked

%% LABELS
labelnames %!IDXCLS
nbehaviors 
ntimelines %!IDXCLS
iLbl2iCls % 2*nclassifiers-by-1 array. iCls = iLbl2iCls(iLbl) where 
% iLbl/iCls reference .labelnames/.classifiernames resp. !IDXCLS
iCls2iLbl % nclassifiers-by-1 cell array, each el is 1-by-2 array.
% [iLblPos iLblNeg] = iCls2iLbl{iCls} where iLblPos, iLblNeg index
% .labelnames. !IDXCLS
labels % labels(expi).t0s, etc
labelidx %!IDXCLS
labelidx_off
labelstats % labelstats(expi).nflies_labeled

%% DATA/FEATURES
perframedata % last-used per-frame data (one fly)
perframeunits % units for perframedata
windowdata % windowdata(iCls).X,  computed and cached window features !IDXCLS
windowfeaturesparams % windowfeaturesparams{iCls} is: struct from curperframefns->[structs containing window feature parameters] !IDXCLS
windowfeaturescellparams % windowfeaturescellparams{iCls} is: % struct from curperframefns->[cell array of window feature parameters] !IDXCLS
allperframefns % Common to all classifiers
curperframefns % curperframefns{iCls} !IDXCLS
scoreFeatures  % A structure array with number of elements equal to the % number of score features
savewindowdata % 1-by-nclassifiers double !IDXCLS
loadwindowdata % 1-by-nclassifiers double !IDXCLS

%% PREDICTIONS
predictdata % predictdata{expi}{flyi}(timelinei) !IDXCLS
predictblocks % predictblocks(iCls).t0, .t1, .expi, .flies which frames, for which tracks, get predicted when training happens. !IDXCLS
fastPredict % See Predict.fastPredict fastPredict(iCls).<etc> !IDXCLS
predictedidx % cache of self.predictdata, like labelidx for labels !IDXCLS
scoresidx  % Is this really an index?  Isn't it just the score itself? !IDXCLS
scoresidx_old %!IDXCLS
scoreTS % timestamp for the above. ALTODO: does not appear to be used for anything !IDXCLS
erroridx % whether the predicted label matches the true label. 0/1/2 !IDXCLS
suggestedidx % predictedidx for unlabeled data, same representation as labelidx !IDXCLS

%% CLASSIFIER/PP
classifiertype % 1-by-nclassifiers cellstr !IDXCLS
classifier % 1-by-nclassifiers cell !IDXCLS
classifier_old % !IDXCLS
lastFullClassifierTrainingSize
classifierTS  % 1-by-nclassifiers vector of default time stamps. Number of days since Jan 0, 0000 (typically non-integer) !IDXCLS
trainstats % 1-by-nclassifiers cell array !IDXCLS
classifier_params % 1-by-nclassifiers cell !IDXCLS
postprocessparams % 1-by-nclassifiers cell !IDXCLS
confThresholds % nclassifiers-by-2 array !IDXCLS

%% UI
labelGraphicParams
trxGraphicParams
labelcolors %!IDXCLS
unknowncolor

%% SIMILARFRAMES
% These have not been updated for multiple classifiers.
% SimilarFrames/Bagging functionality should all be asserted to run
% only on single-classifier projects.
frameFig % Show similar frames
distMat % Show similar frames
bagModels % Show similar frames
fastPredictBag % Show similar frames, fastPredictBag.classifier

%% GT MODE
gtMode
randomGTSuggestions % randomGTSuggestions{iExp}(iFly).start, .end (scalars)
thresholdGTSuggestions
loadedGTSuggestions % same format/structure as randomGTSuggestions, except .start/.end can be vectors
balancedGTSuggestions
GTSuggestionMode
otherModeLabelsEtc  % A place to store things for the other mode.

%% MISC IMPL
everythingFileNameAbs

defaultpath % last-used path for loading jabfiles
expdefaultpath % last-used path for loading experiments
setstatusfn    % functions for writing text to a status bar
clearstatusfn

version
isInteractive
deterministic = false; % scalar logical or double, for testing purposes.

doUpdate % Retrain properly % AL: does not appear to be used
cacheSize % ALTODO comments suggests obsolete

thereIsAnOpenFile
userHasSpecifiedEverythingFileName
needsave  % true iff a file is open and there are unsaved changes.  False if no file is open.
