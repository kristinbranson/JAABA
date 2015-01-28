function classifier = trainDetectorSTCore(expdirs,labels,behname,nobehname,usePastOnly)
% classifier = trainDetectorSTCore(expdirs,labels,lblnames,usePastOnly)
% expdirs: cellstr 
% labels: .labels field from JLabelData
% behname: char. behavior name of interest. 
% nobehname: char. no-behavior name corresponding to behname. 
%
% classifier: trained classifier 
%
% Note: behname and nobehname must match (case-sensitive) bout names in
% labels!

FEATUREWINDOWRADIUS = 2;

if usePastOnly
  fprintf('Training a classifier that only uses past information...\n');
end

X = [];
L = zeros(1,0);

% Accumulate features and labels into X and L
for expi = 1:numel(expdirs)
  ftrsfname = fullfile(expdirs{expi},'features.mat');
  assert(exist(ftrsfname,'file')==2,'Missing features.mat file ''%s''.',ftrsfname);
  S = load(ftrsfname);
  
  for fly = 1:numel(labels(expi).flies)
    for bnum = 1:numel(labels(expi).t0s{fly})
      
      name = labels(expi).names{fly}{bnum};
      if ~any(strcmp(name,{behname nobehname}))
        % Bout irrelevant to current classifier; ignore.
        % This looping is theoretically inefficient.
        continue;
      else
        t0 = labels(expi).t0s{fly}(bnum);
        t1 = labels(expi).t1s{fly}(bnum)-1;
        
        % only train on parts of the bout that have information
        if usePastOnly,
          if t1 <= FEATUREWINDOWRADIUS+2,
            continue;
          elseif t0 <= FEATUREWINDOWRADIUS+2,
            t0 = FEATUREWINDOWRADIUS+2;
          end
          t0past = t0-FEATUREWINDOWRADIUS-1;
          t1past = t1-FEATUREWINDOWRADIUS-1;
          curFeatures = S.curFeatures(t0past:t1past,:);
        else
          curFeatures = S.curFeatures(t0:t1,:);
        end
        X = [X; curFeatures]; %#ok<AGROW>
        oldCondition = strcmpi(name,'None') || strncmp(name,'No_',3);
        isNoBeh = strcmp(name,nobehname);
        assert(isNoBeh==oldCondition);
        if isNoBeh
          curL = 2*ones(1,t1-t0+1);
        else
          curL = ones(1,t1-t0+1);
        end
        L = [L curL]; %#ok<AGROW>
      end
    end
  end
end

classifier_params = ...
  struct('iter',100, ...
  'iter_updates',10, ...
  'numSample',2500, ...
  'numBins',30, ...
  'CVfolds',7, ...
  'baseClassifierTypes',{'Decision Stumps'}, ...
  'baseClassifierSelected',1);

binVals = findThresholds(X,classifier_params);
bins = findThresholdBins(X,binVals);

classifier = boostingWrapper(X, ...
  L',[],...
  binVals,...
  bins, ...
  classifier_params);

