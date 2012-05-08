function importData(labelFilename,expDir, featureFile, varargin)

[analysis_protocol,settingsdir] = myparse(varargin,...
  'analysis_protocol','current',...
  'settingsdir','settings');

    function behaviors = safelyAddBehavior(behaviors, newBehavior)
        if ~ismember(newBehavior, behaviors)
            behaviors = [behaviors {newBehavior}];
        end
    end
    function [segstarts segends labels] = removeLabelSegment(id, segstarts, segends, labels)
        segstarts(id) = [];
        segends(id) = [];
        labels(id) = [];
    end
    function [t0 t1 segstarts segends labels] = getLabeledWindow(t0tolabelcurr, t1tolabelcurr, segstarts, segends, labels)
        t0 = t0tolabelcurr;
        t1 = t1tolabelcurr;
        if t0 == segstarts(1) % do not use labels that may have started right before...
            t0 = segstarts(2);
            [segstarts segends labels] = removeLabelSegment(1, segstarts, segends, labels);
        end
        if t1 == segends(end) % ... or ended right after the labeling window
            t1 = segends(end-1);
            [segstarts segends labels] = removeLabelSegment(length(segends), segstarts, segends, labels);
        end
    end
    function [NEWbehaviors NEWsegstarts NEWsegends NEWlabels] = convertLabelsCWAU(behaviors, t0, t1, segstarts, segends, labels)
        [~,labelIdx] = ismember(labels, behaviors);
        t1 = t1 - t0 + 1;
        segstarts = segstarts - t0 + 1;
        segends = segends - t0 + 1;
        seq = NaN(t1, 1);
        for si = 1 : length(segstarts)
%            seq(max(1, segstarts(si)) : min(t1, segends(si))) = labelIdx(si);
            seq(segstarts(si) : segends(si)) = labelIdx(si);
        end
        
        [~,idxUnknown] = ismember('unknown', behaviors); % unknown is last label in behavior set
%        NEWbehaviors = [behaviors(1:idxUnknown-1) {'none'}];
        NEWbehaviors = safelyAddBehavior(behaviors, 'none'); % already added by SVMBehaviorParams.txt
        
        allUnlabeled = isnan(seq);
        [~,idxNone] = ismember('none', behaviors);
        seq(allUnlabeled) = idxNone;
%        allUnknown = seq(seq == idxUnknown);
%        idxNone = length(NEWbehaviors);
%        seq(unlabeled) = idxNone;
        grad = [seq; 0] - [0; seq];
        pos = find(grad ~= 0);
        
        numSegments = length(pos) -1;
        NEWlabelIdx = seq(pos(1 : numSegments));
        
        NEWsegstarts = NaN(1, numSegments);
        NEWsegends = NaN(1, numSegments);
        NEWlabels = cell(1, numSegments);
        for si = 1 : numSegments
            NEWsegstarts(si) = pos(si);
            NEWsegends(si) = pos(si+1)-1;
            NEWlabels(si) = NEWbehaviors(NEWlabelIdx(si));
        end
        NEWsegstarts = NEWsegstarts + t0 - 1;
        NEWsegends = NEWsegends + t0 - 1;
        
        NEWbehaviors(idxUnknown) = [];
        del = NEWlabelIdx == idxUnknown;

        NEWsegstarts(del) = [];
        NEWsegends(del) = [];
        NEWlabels(del) = [];
    end


%% build filenames
    
    [path featureFilename ext] = fileparts(featureFile);
    outputDirPrefix = fileparts(path);
    
%    pathDirs = textscan(path, '%s', 'delimiter', filesep);
%     pathDirs = textscan(path, '%s', 'delimiter', '\\');
%     outputDirPrefix = [];
%     for i = 1 : length(pathDirs{1})-1
%         outputDirPrefix = [outputDirPrefix pathDirs{1}{i} filesep];
%     end
    outputDir = fullfile(outputDirPrefix,'TrainingData');
    
    
    featureName = featureFilename(1 : strfind(featureFilename,'SVMFeatureParams')-1);
    behaviorFilename = [featureName 'SVMBehaviorParams'];
    behaviorSet = loadBehaviorFile([path filesep behaviorFilename ext]);

%% read expected labels and features

    behaviors = {}; % set of valid 'to-be-expected' labels
    fns = fieldnames(behaviorSet);
    numBehaviors = length(fns);
    for fi = 1 : numBehaviors
        currBehavior = behaviorSet.(fns{fi});
        for bi = 1 : length(currBehavior)
            behaviors = [behaviors {currBehavior(bi).short}];
        end
    end
    
    
    
    featureFileHandle = fopen(featureFile, 'r');
    featureFileContent = textscan(featureFileHandle, '%s', 'delimiter', ':');
    featureNames = {};
    for i = 1 : 2 : length(featureFileContent{1})
        featureNames = [featureNames; featureFileContent{1}(i)];
    end
    


%% Default Values taken from ~/projects/fromAlice/code/LabelLocomotionBehaviorScript.m
%    behaviors = {'walk','stop','sharpturn','crabwalk','backup','jump','pause','wallwalk','groom','unknown'};

  %CSC 20110321: depricated 
    fnsignore = {
      'x'
    'y'
    'a'
    'b'
    'id'
    'moviename'
    'firstframe'
    'arena'
    'off'
    'nframes'
    'endframe'
    'matname'
    'sex'
    'pxpermm'
    'fps'
    'annname'
    'units'
  }';


    behaviors2 = safelyAddBehavior(behaviors, 'unknown'); % case sensitive!
%    behaviors2 = safelyAddBehavior(behaviors2, 'tracker failure');
    
%   CSC 20110321: depricated, filename is now set in <settingsDir>/<analysis_protocol>/dataloc_params.txt, line 'trxfilestr'
%
     %tokenpos = strfind(labelFilename, '_');
     %inputTrxPath = fileparts(labelFilename);
     %inputTrxDir = labelFilename(tokenpos(2)+1 : tokenpos(end-1)-1);    
%     inputTrxFilename = sprintf('%s/%s/fixed_ctrax_results.mat', inputTrxPath, inputTrxDir);       

%   CSC 20110321: depricated, values now genetated by FlyBowlComputePerFrameFeatures (rather than directlyloaded)
%
%     % these variables are loaded from a registered_trx.mat file
%     
%     timestamps = NaN;
%     trx = NaN;
%
%    load(inputTrxFilename); % sets trx values above

    % these variables are loaded from the labeled*.mat-file
    
    behaviorcolors = NaN;
    behaviorcurr = NaN;
    fly = NaN;
    labelname = NaN;
    labels = NaN;
    moviename = NaN;
    propsmatname = NaN;
    segends = NaN;
    segstarts = NaN;
    t0tolabelcurr = NaN;
    t1tolabelcurr = NaN;
    
    load(labelFilename); % sets label values above
    
    usedLabels = unique(labels);
    keep = ismember(usedLabels,behaviors2);
    if ~all(keep)
        error('Potential mismatch: Labeled behavior not in specified behavior set. Aborting...\n');
    end
    
    %expDir = [inputTrxPath filesep inputTrxDir];
%   UNCOMMENT THIS PARAGRAPH !!    
    FlyBowlRegisterTrx(expDir, ...
      'analysis_protocol',analysis_protocol, ...
      'settingsdir',settingsdir,...
      'dotemporalreg',false);

%                    
%     FlyBowlClassifySex(expDir, ...
%                        'analysis_protocol','AliceFixed2Steve', ...
%                        'settingsdir', '~/projects/FlyBowlAnalysis/settings');
                   
    
    trx = FlyBowlComputePerFrameFeatures(expDir, 'analysis_protocol', analysis_protocol, ...
                                                      'settingsDir', settingsdir,...
                                                      'perframefns', featureNames);
                                                  

    
    [behaviors3 segstarts segends labels] = convertLabelsCWAU(behaviors2, t0tolabelcurr, t1tolabelcurr, segstarts, segends, labels);
    [t0 t1 segstarts segends labels] = getLabeledWindow(t0tolabelcurr, t1tolabelcurr, segstarts, segends, labels);
                                                                                                   
    [output, filename ext] = fileparts(labelFilename);
    outputTrxFilename = [outputDir filesep, filename ext '.trx'];
%    outputTrxFilename = sprintf('%s.trx', labelFilename);
    outputTrx2(trx, fly, moviename, propsmatname, outputTrxFilename, t0, t1, featureNames);
%    outputTrx(trx, fly, moviename, propsmatname, outputTrxFilename, trx(fly).firstframe, trx(fly).endframe, fnsignore);     
%    outputLabels(segstarts,segends,labels,behaviors,idx,moviename,matname,trxname,savename,t0,t1,flies)

    
%    [~, filename ext] = fileparts(labelFilename);
    %outputLabelFilename = sprintf('%s.label', labelFilename);
    outputLabelFilename = [outputDir filesep filename ext '.label'];
    [~, outputTrxFile, outputTrxExt] = fileparts(outputTrxFilename);
%    outputLabels(segstarts, segends, labels, {behavior}, fly, moviename, propsmatname, [outputTrxFile outputTrxExt], outputLabelFilename)

%     t0 = 1;
%     endframes = trx.endframes;
%     t1 = endframes(fly);
    
    outputLabels(segstarts, segends, labels, behaviors3, fly, moviename, propsmatname, [outputTrxFile outputTrxExt], outputLabelFilename, t0, t1)
%   use this call instead in case data labeled as 'unknown' should remain unlabeled and therefore become 'Unknown' in Steve's code
%    outputLabels(segstarts, segends, labels, behaviors2, fly, moviename, propsmatname, [outputTrxFile outputTrxExt], outputLabelFilename, t0, t1)
end