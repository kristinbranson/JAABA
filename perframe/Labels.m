classdef Labels 
  % Label-related datastructures
  %
  % labelIdx 
  % A scalar struct containing labels for a single track across multiple 
  % timelines/classifiers. Key props are .vals, .imp, .timestamp which 
  % are numTimelines-by-nsamps. 
  % 
  % labelsShort
  % A scalar struct containing labels for a single track across multiple 
  % timelines. Key props are .names, .t0s, .t1s, .imp. labelsShort is 
  % simlar to the labels structure and is in "bout" format.
  % 
  % labels
  % A struct array indexed by experiment and fly, ie labels(iExp){iFly}. 
  % Key props are .names, .t0s, .t1s. Labels are stored by bout.
  
  methods (Static) % labelidx methods
 
    function labelidx = labelIdx(labelnames,T0,T1)
      % labelIdx 'Constructor' 
      %
      % labelidx.vals - nTL-by-nFrm. A value of 0 indicates 'unlabeled'. 
      % Otherwise, the allowed values are the label indices, ie the values
      % 1:labelidx.nbeh.
      
      [nTL,idxBeh2idxTL,tl2IdxBeh] = Labels.determineNumTimelines(labelnames);
      
      n = T1-T0+1;
      off = 1 - T0;
      
      labelidx = struct();      
      labelidx.T0 = T0;
      labelidx.T1 = T1;
      labelidx.n = n;
      labelidx.off = off;
      labelidx.labelnames = labelnames;
      labelidx.nbeh = numel(labelidx.labelnames);            
      labelidx.nTL = nTL;
      labelidx.idxBeh2idxTL = idxBeh2idxTL;
      labelidx.TL2idxBeh = tl2IdxBeh;
      labelidx.vals = zeros(nTL,n);
      labelidx.imp = zeros(nTL,n);
      labelidx.timestamp = zeros(nTL,n);
    end
    
    function labelidx = labelIdxInit(labelidx,labelsShort)
      % Initialize from labelsShort struct
      % initialized props: .vals, .timestamp, .imp
      
      tfMultiCls = labelidx.nTL>1;
      
      for iBeh = 1:labelidx.nbeh
        beh = labelidx.labelnames{iBeh};
        iTL = labelidx.idxBeh2idxTL(iBeh);
        iLbls = find(strcmp(labelsShort.names,beh));
        assert(isequal(numel(labelsShort.names),numel(labelsShort.t0s),numel(labelsShort.t1s)));
        for iLbl = iLbls(:)'
          t0 = labelsShort.t0s(iLbl);
          t1 = labelsShort.t1s(iLbl);
          if t0>labelidx.T1 || t1<labelidx.T0
            continue;
          end
          t0 = max(labelidx.T0,t0);
          t1 = min(labelidx.T1+1,t1);
          idx = t0+labelidx.off:t1-1+labelidx.off;
          labelidx.vals(iTL,idx) = iBeh;
          if tfMultiCls
            labelidx.imp(iTL,idx) = 1; % ALXXX EXTENDED: for multicls, all bouts important
          end
          labelidx.timestamp(iTL,idx) = labelsShort.timestamp(iLbl); 
        end
      end
      if ~tfMultiCls
        % ALXXX EXTENDED. Single classifier, legacy code for writing 
        % importance 
        assert(numel(labelsShort.imp_t0s)==numel(labelsShort.imp_t1s));
        for iLbl = 1:numel(labelsShort.imp_t0s)
          t0 = labelsShort.imp_t0s(iLbl);
          t1 = labelsShort.imp_t1s(iLbl);
          if t0>labelidx.T1 || t1<labelidx.T0
            continue;
          end
          t0 = max(labelidx.T0,t0);
          t1 = min(labelidx.T1+1,t1);
          idx = t0+labelidx.off:t1-1+labelidx.off;

          % iTL = 1:labelidx.nTL;
          % warningNoTrace('Labels:labelidx:imp','labels.imp_t0s, .imp_t1s currently always map to ALL TLs.');
          labelidx.imp(1,idx) = 1;
        end
      end
    end
    
    function label012 = labelVec2label012(v,vPos,vNeg)
      % Convert a labeling vector from 0/vPos/vNeg to 0/1/2, ie remap
      % positive labels to 1 and negative values to 2

      assert(vPos~=vNeg && vPos~=0 && vNeg~=0);
      assert(isvector(v));
      assert(all(v==vPos | v==vNeg | v==0));
      
      label012 = zeros(size(v));
      label012(v==vPos) = 1;
      label012(v==vNeg) = 2;      
    end
    
  end
  
  methods (Static) % labels/labelShort methods
    
    % See JLabelData, 'labels' property, for description
    
    % AL 20140911 TODO
    % Current defn of .timelinetimestamp is awkward here in that only
    % "real" behaviors are included. Would be more natural and easier for
    % Labels.timelinetimestamp to include all labels, including
    % No-behaviors. Now currently in JLabel these timestamps may not
    % necessarily be 100% faithfully set, depending on implementation (eg
    % due to the "erasure problem" if a bout is erased it is unknown what
    % was erased), so in practice the behavior- and No-behavior timestamps
    % may always be the same. Nonetheless I think this would be an
    % improvement here.
    
%     ALXXX EXTENDED, defer until working on GT mode, when can confirm
%     usage of importance info. Until then:
%     - for single classifiers, importance captured in .imp_t0 and .imp_t1,
%     as today
%     - for multi-classiifers, all bouts assumed to be important.
% 
%     - do the above note, .timelinetimestamp -> .labeltimestamp
%     - add a prop, .allNames or similar, that contains all names/labels
%     that are allowed in .names (set of classifiers)
%     - add .imp, int vec same size as t0s/t1s giving importance (or whatever) level
%       of each bout. Only gotcha here might be a bout that has multiple
%       importance levels, eg first part important second half not; this
%       would need to get split up into two bouts which maybe isn't
%       convenient.
%     - .timestamp still present but not used and essentially deprecated
%     - .imp_t0s, .imp_t1s still present for legacy labels, but incorporated into .imp, 
%        and plan to deprecate moving forward. definitely, do not use .imp_t0s or .imp_t1s 
%        for multiclass.
    
    function labels = labels(nexp)
      % labels constructor
      
      % Data structure
      % .t0s, .t1s, .names, .timestamp, .imp_t0s, .imp_t1s are bout labels:
      % labels(iExp).t0s{iFly}.t0s(iBout)
      % - where iExp runs from 1:numel(labels) 
      % - where iFly runs from 1:numel(labels(iExp).flies)
      % - where iBout runs from 1:numel(labels(iExp).t0s{iFly}.t0s). 
      %
      % .flies is a column vector of length numflies-for-exp.
      % .off is a row vector of length numflies-for-exp.
      % 
      % .labels(iExp).timelinetimestamp{iFly} is a structure with
      % fieldnames as behaviors and values as
      % timestamp-that-the-labeling-timeline-for-that-behavior/target-was-last-edited.
      % This structure is guaranteed to contain fields for all 'real'
      % behaviors. Note, a target may not be present in labels (eg a fly in
      % some experiment has no labels), in which case there may be no
      % timelinetimestamp available. On the other hand, if a bout is
      % created and then erased, a timelinetimestamp may be present.
      labels = struct(...
        't0s',repmat({cell(1,0)},1,nexp),...
        't1s',{cell(1,0)},...
        'names',{cell(1,0)},...
        'flies',zeros(0,1),... 
        'off',zeros(1,0),...
        'timestamp',{cell(1,0)},...
        'timelinetimestamp',{cell(1,0)},...
        'imp_t0s',{cell(1,0)},...
        'imp_t1s',{cell(1,0)});      
    end
    
    function labels = initTimelineTimestamps(labels,realbehnames)
      % Initialize timelineTimestamps to 'now' for all exps/flies/behaviors 
      Nrealbeh = numel(realbehnames);
      nowvec = repmat({now},1,Nrealbeh);
      timelineTSVal = cell2struct(nowvec,realbehnames(:)',2);
      
      Nexp = numel(labels);
      for iExp = 1:Nexp
        Nfly = numel(labels(iExp).flies);
        for iFly = 1:Nfly
          labels(iExp).timelinetimestamp{iFly} = timelineTSVal;
        end
      end
    end 
    
    function [labels,tfModified] = modernizeLabels(labels,realbehnames)
      % Bring labels data structure up-to-date.
      % realbehnames: cellstr of real behavior names.
      % tfModified: true if a change was made
      % 
      % Changes to labels data structure:
      % 20140825: added .timelinetimestamp in order to optimize
      % multiclassifier training

      if isempty(labels) && ~isequal(labels,Labels.labels(0))
        % AL 20150709: some legacy .gtLabels do not have standard/expected 
        % fields        
        labels = Labels.labels(0);
        tfModified = true;
        return;
      end
      
      tfModified = false;
      
      Nexp = numel(labels);
      
      if ~isfield(labels,'timelinetimestamp')
        assert(~isempty(labels));
        for iExp = 1:Nexp
          Nfly = numel(labels(iExp).flies);
          labels(iExp).timelinetimestamp = cell(1,Nfly);
        end
        tmp = Labels.labels(0); % Just to get fieldnames in current/proper order
        labels = orderfields(labels,tmp);

        labels = Labels.initTimelineTimestamps(labels,realbehnames);
        
        tfModified = true;
      end
      
      % AL 20140902, some jabs got saved with labels.timestamp==0. This is
      % probably the issue addressed by SHA 73a8f. The zero timestamp(s) 
      % trigger an assert in JLabelData.
      for iExp = 1:Nexp
        Nfly = numel(labels(iExp).flies);
        for iFly = 1:Nfly
          ts = labels(iExp).timestamp{iFly};
          if any(ts==0)
            warningNoTrace('Labels:zeroTimestamp',...
              'Experiment %d, Fly %d: found timestamp=0. Setting zero timestamps to now (%s) and resetting timelinetimestamps.',...
              iExp,iFly,datestr(now));
                        
            % reset .timestamp for this experiment, fly. In theory all
            % bouts could be affected, not just the one that is zero.
            nowval = now;
            ts(:) = nowval;
            labels(iExp).timestamp{iFly} = ts;
            
            % reset .timelinetimestamp for this experiment, fly.
            behnames = fieldnames(labels(iExp).timelinetimestamp{iFly});
            assert(isequal(sort(behnames(:)),sort(realbehnames(:))));            
            for beh = behnames(:)',beh=beh{1}; %#ok<FXSET>
              labels(iExp).timelinetimestamp{iFly}.(beh) = nowval;
            end
            
            tfModified = true;
          end
        end
      end
    end
    
    function verifyLabels(labels)
      % Confirm labels invariants
      for i = 1:numel(labels)
        l = labels(i);
        
        nseq = numel(l.flies);
        if nseq>0
          assert(isequal(size(l.t0s),[1 nseq]));
          assert(isequal(size(l.t1s),[1 nseq]));
          assert(isequal(size(l.names),[1 nseq]));
          assert(isequal(size(l.flies),[nseq 1]));
          assert(isequal(size(l.off),[1 nseq]));
          assert(isequal(size(l.timestamp),[1 nseq]));
          assert(isequal(size(l.timelinetimestamp),[1 nseq]));
          assert(isequal(size(l.imp_t0s),[1 nseq]));
          assert(isequal(size(l.imp_t1s),[1 nseq]));
        else
          % Legacy labels won't have right-sized empties
          % TODO: do something here maybe just check empty
        end          
      end
    end
    
    function labels = removeClassifier(labels,behname,nobehname)
      % Remove all bouts for behavior behname and no-behavior nobehname;
      % remove timelinetimestamps.
      
      Labels.verifyLabels(labels);
      
      behnobeh = {behname nobehname};
      for iExp = 1:numel(labels)
        for iFly = 1:numel(labels(iExp).flies)
          tfKeep = ~ismember(labels(iExp).names{iFly},behnobeh);
          
          labels(iExp).t0s{iFly} = labels(iExp).t0s{iFly}(tfKeep);
          labels(iExp).t1s{iFly} = labels(iExp).t1s{iFly}(tfKeep);
          labels(iExp).names{iFly} = labels(iExp).names{iFly}(tfKeep);
          labels(iExp).timestamp{iFly} = labels(iExp).timestamp{iFly}(tfKeep);
          labels(iExp).timelinetimestamp{iFly} = rmfield(...
            labels(iExp).timelinetimestamp{iFly},behname);
        end
      end
    end
    
    function labels = addClassifier(labels,behname,nobehname) %#ok<INUSD>
      Labels.verifyLabels(labels);
      for iExp = 1:numel(labels)
        for iFly = 1:numel(labels(iExp).flies)
          assert(~isfield(labels(iExp).timelinetimestamp{iFly},behname));
          labels(iExp).timelinetimestamp{iFly}.(behname) = now;
        end
      end
    end
      
    function labels = clearLabels(labels,behname,realbehname)
      % Clears all bouts of behavior <behname> corresponding to
      % realbehavior <realbehname>. 
      % Timelinetimestamp is updated
      
      % AL: doesn't clear imp_t0s, imp_t1s
      
      Labels.verifyLabels(labels);
      
       for iExp = 1:numel(labels)
        for iFly = 1:numel(labels(iExp).flies)
          tfKeep = ~strcmp(behname,labels(iExp).names{iFly});
          
          labels(iExp).t0s{iFly} = labels(iExp).t0s{iFly}(tfKeep);
          labels(iExp).t1s{iFly} = labels(iExp).t1s{iFly}(tfKeep);
          labels(iExp).names{iFly} = labels(iExp).names{iFly}(tfKeep);
          labels(iExp).timestamp{iFly} = labels(iExp).timestamp{iFly}(tfKeep);          
          labels(iExp).timelinetimestamp{iFly}.(realbehname) = now;
        end
       end
    end
    
    function labels = renameBehavior(labels,behOld,behNew,nobehOld,nobehNew)
      % behOld: old behaviorname
      % behNew: new "
      % nobehOld: old no-behaviorname
      % nobehNew: new "
      %
      % This renames behOld->behNew and nobehOld->nobehNew throughout
      % labels. It also updates the fieldnames of timelinetimestamps. The
      % values of timelinetimestamps (eg the timestamps themselves) are not
      % modified.
      %
      % This does not do error checking, behOld should be present in
      % timelinetimestamps etc.
      
      for iExp = 1:numel(labels)
        for iFly = 1:numel(labels(iExp).flies)
          % rename bout lbls
          tf1 = strcmp(labels(iExp).names{iFly},behOld);
          tf2 = strcmp(labels(iExp).names{iFly},nobehOld);
          labels(iExp).names{iFly}(tf1) = {behNew};
          labels(iExp).names{iFly}(tf2) = {nobehNew};
          
          % rename timelinetimestamp
          if strcmp(behOld,behNew)
            % no operation necessary
          else
            timelineTS = labels(iExp).timelinetimestamp{iFly};
            assert(isfield(timelineTS,behOld));
            assert(~isfield(timelineTS,behNew));
            val = timelineTS.(behOld);
            timelineTS.(behNew) = val;
            labels(iExp).timelinetimestamp{iFly} = rmfield(timelineTS,behOld);
          end
        end
      end      
    end
      
    function labels = renameBehaviorRaw(labels,behOld,behNew)
      % behOld: old label
      % behNew: new label
      %
      % Replace label behOld with behNew throughout labels. behOld/behNew
      % may represent either real behaviors or no-behaviors.
      % WARNING: timelinetimestamps are unmodified. See
      % Lables.renameBehavior for a safer rename method.
      assert(ischar(behOld));
      assert(ischar(behNew));

      for iExp = 1:numel(labels)
        for iFly = 1:numel(labels(iExp).flies)
          tf = strcmp(labels(iExp).names{iFly},behOld);
          labels(iExp).names{iFly}(tf) = {behNew};
        end
      end
    end
    
    function tts = mostRecentTimelineTimestamps(labels)
      % tts = mostRecentTimelineTimestamps(labels)
      % Get most recent timeline timestamps for behaviors/classifiers,
      % across all experiments/flies.
      %
      % tts: scalar struct. Fieldnames are "real" behaviors. Values are
      % most recent timestamps for labeling marks made relevant to that
      % behavior/classifier.
      %
      % Note, not all behaviors will appear in tts. If a behavior does not
      % appear in tts, it means that its most recent timestamp is
      % indeterminate/unknown.

      tts = struct();
      
      Nexp = numel(labels);
      for iExp = 1:Nexp
        Nfly = numel(labels(iExp).flies);
        for iFly = 1:Nfly
          timelineTS = labels(iExp).timelinetimestamp{iFly};
          flds = fieldnames(timelineTS);
          for f = flds(:)',f=f{1}; %#ok<FXSET>
            val = timelineTS.(f);
            if ~isfield(tts,f)
              tts.(f) = val;
            else
              tts.(f) = max(tts.(f),val);
            end
          end
        end
      end
    end
    
    function tf = labelsSeen(labels,names)
      % names: cellstr
      % tf: logical array, same size as names. tf(i) is true if names{i} is
      % in at least one bout in labels.
      
      % MERGESTUPDATED
      
      assert(iscellstr(names));
      tf = false(size(names));      
      for i = 1:numel(labels)
        labelsCell = labels(i).names;
        assert(iscell(labelsCell));
        namesExp = [labelsCell{:}]; % concatenate bouts/names for all flies in exp
        if isempty(namesExp)
          namesExp = cell(1,0);
        end
        namesExp = unique(namesExp);
        tf = tf | ismember(names,namesExp);
        if all(tf)
          break;
        end
      end
    end
    
    function labelsComb = compileLabels(combExpDirNames,labels,expDirNames,labelnames)
      % Compile/combine multiple label struct arrays.
      % 
      % combExpDirNames: cellstr, desired/target master list of experiments
      % labels: struct array of Labels. 
      % expDirNames: cellstr of expdirnames for labels.
      % labelnames: cellstr of desired/expected "total set" of 
      %   labelnames (in standard order, <beh> <nobeh>). Used to init
      %   timelinetimestamps etc.
      % 
      % labelsComb: struct array of labels structures labeled by
      % combExpDirNames that contains all bouts in labels.
      %
      % This method works by iterating over labels and applying its bouts
      % consecutively onto an initially "blank slate".
      % IMPORTANT: The bouts in the resulting labels struct may overlap.
      
      assert(iscellstr(combExpDirNames));
      Labels.verifyLabels(labels);
      
      assert(iscellstr(expDirNames) && numel(expDirNames)==numel(labels));
      [tf,loc] = ismember(expDirNames,combExpDirNames);
      assert(all(tf),'Master list of experiments does not cover all labeled experiments.');
      
      % Initialize combLabels
      NCombExps = numel(combExpDirNames);
      labelsComb = Labels.labels(NCombExps);
      
      % Loop over input labels and add bouts to combLabels
      Nlbl = numel(labels);
      realBehsSeen = cell(0,1);
      for i = 1:Nlbl
        lbls = labels(i);
        iExpComb = loc(i);
        
        for iFly = 1:numel(lbls.flies)
          fly = lbls.flies(iFly);
          off = lbls.off(iFly);
          Nbout = numel(lbls.t0s{iFly});
          for iBout = 1:Nbout
            labelsComb = Labels.addBout(labelsComb,iExpComb,fly,...
              lbls.t0s{iFly}(iBout),...
              lbls.t1s{iFly}(iBout),...
              lbls.names{iFly}{iBout},...
              off,...
              lbls.timestamp{iFly}(iBout));
          end
          
          timelineTS = lbls.timelinetimestamp{iFly};
          realBehsSeen = union(realBehsSeen,fieldnames(timelineTS));
          % fieldnames for .timelinetimestamp is probably the same for 
          % all flies in lbls but to be safe etc
        end
      end
      
      clsNames = Labels.verifyBehaviorNames(labelnames);
      assert(all(ismember(realBehsSeen,clsNames)),...
        'Behavior names encountered that are not included in supplied ''labelnames''.');      
      labelsComb = Labels.initTimelineTimestamps(labelsComb,clsNames);      
      labelsComb = Labels.removeOverlappingBoutsRaw(labelsComb,labelnames);
    end
    
    function labels = removeOverlappingBoutsRaw(labels,labelnames)
      % Remove overlapping bouts by transforming
      % labels->labelsShort->labelIdx->labelsShort->labels, ie write bouts
      % to timelines, and then back.
      %
      % labelnames: cellstr of labelnames in standard ordering (<behs> then
      % <nobehs>). This is just for implementation convenience and is
      % theoretically unnecessary.
      %
      % IMPORTANT: 
      % * .imp_t0s, .imp_t1s may not be set correctly for multiclassifier labels.
      % * timelinetimestamp is not substantively updated.
      
      Nlbl = numel(labels);
      for i = 1:Nlbl
        nFly = numel(labels(i).flies);
        for iFly = 1:nFly
          % get labelsShort for this fly
          fly = labels(i).flies(iFly);
          lblShort = Labels.labelsShort();
          [lblShort,tffly] = Labels.labelsShortInit(lblShort,labels(i),fly);     
          assert(tffly);
          
          % create/init a "blank" labelIdx for this fly
          off = labels(i).off(iFly);
          T0 = 1-off;
          T1 = max([lblShort.t0s(:); lblShort.t1s(:); ... 
                    lblShort.imp_t0s(:); lblShort.imp_t1s(:); T0]) + 1; % make T1 as large as necessary to hold all bouts
          lblIdx = Labels.labelIdx(labelnames,T0,T1);

          lblIdx = Labels.labelIdxInit(lblIdx,lblShort);
          lblShort = Labels.labelsShortFromLabelIdx(lblIdx);
          labels(i) = Labels.assignFlyLabelsRaw(labels(i),lblShort,fly);
        end
      end
    end
    
    function m = labelMatrix(labels,T0,T1,behNames)
      % labels: labels array
      % T0,T1: scalar doubles
      % behNames: cellstr
      %
      % m: numel(labels) x nsamp x numel(behNames), nsamp=T0-T1+1. 1/-1 for
      % beh and no-beh, resp.
      %
      % Currently, bouts for all flies/targets within an experiment are
      % superimposed (only because the expected usage is in ST mode).
      
      % AL: in JAABA, could go Labels->LabelsShort->LabelsIdx to implement
      % this transformation.
      
      Labels.verifyLabels(labels);
      
      nExp = numel(labels);
      nBeh = numel(behNames);
      nFrm = T1-T0+1;
      off = 1-T0;
      m = zeros(nExp,nFrm,nBeh);
      
      behnobeh = Labels.behnames2labelnames(behNames);
      behnobeh = reshape(behnobeh,nBeh,2); 
      %behnobeh = [behNames(:) nobehNames(:)];
      
      for iExp = 1:nExp
        L = labels(iExp);
        for iFly = 1:numel(L.flies)          
          tts = L.timelinetimestamp{iFly};
          fns = fieldnames(tts);
          tf = ismember(behNames,fns);
          if ~all(tf)
            warning('Labels:behNotPresent','Exp %d, behaviors not present in labels timelinetimestamp: %s',...
              iExp,civilizedStringFromCellArrayOfStrings(behNames(~tf)));
          end
          
          for iBout = 1:numel(L.names{iFly})
            t0 = L.t0s{iFly}(iBout);
            t1 = L.t1s{iFly}(iBout);
            name = L.names{iFly}{iBout};
            if t0>T1 || t1<T0
              continue;
            end
            tf = strcmp(name,behnobeh);
            if nnz(tf)==0
              continue;
            end
            assert(nnz(tf)==1,'Label ''%s'' matches multiple reference labels.',name);
            [iTL,behOrNo] = find(tf);
            
            t0 = max(T0,t0);
            t1 = min(T1+1,t1);
            idx = t0+off:t1-1+off;
            if behOrNo==1 % positive behavior
              if any(m(iExp,idx,iTL)<0)
                warning('Labels:conflictingLabels','Conflicting labels across flies for exp %d, fly %d, bout %d (%s).',...
                  iExp,iFly,iBout,name);
              end
              m(iExp,idx,iTL) = 1;
            else % no-beh
              if any(m(iExp,idx,iTL)>0)
                warning('Labels:conflictingLabels','Conflicting labels across flies for exp %d, fly %d, bout %d (%s).',...
                  iExp,iFly,iBout,name);
              end
              m(iExp,idx,iTL) = -1;
            end
          end
        end
      end
    end
    
    function bouts = boutList(labels,expNames)
      % Return a struct array where each element represents a distinct bout
      % in labels.
      
      assert(numel(labels)==numel(expNames));
     
      bouts = struct('id',cell(0,1),'t0',[],'t1',[],'name',[],'exp',[],'fly',[],'flyoff',[],'timestamp',[]);      
      for iExp = 1:numel(labels)
        exp = expNames{iExp};
        lblsExp = labels(iExp);
        for iFly = 1:numel(lblsExp.flies)
          fly = lblsExp.flies(iFly);
          flyoff = lblsExp.off(iFly);          
          for iBout = 1:numel(lblsExp.t0s{iFly})
            t0 = lblsExp.t0s{iFly}(iBout);
            t1 = lblsExp.t1s{iFly}(iBout);
            nm = lblsExp.names{iFly}{iBout};
            
            id = sprintf('%s|%d|%s|%d|%d',exp,fly,nm,t0,t1);

            bouts(end+1,1).id = id; %#ok<AGROW>
            bouts(end).exp = exp; 
            bouts(end).fly = fly;
            bouts(end).flyoff = flyoff;
            bouts(end).t0 = t0;
            bouts(end).t1 = t1;
            bouts(end).name = nm;
            bouts(end).timestamp = lblsExp.timestamp{iFly}(iBout);
          end
        end
      end
    end
    
    function labels = assignFlyLabelsRaw(labels,labelsShort,fly)
      % Assign bouts for single target (labelsShort) into labels
      %
      % labels: scalar labels struct
      % labelShort: labelsShort for a single target
      % fly: target ID (value, not target index)
      %
      % IMPORTANT: labels.timelinetimestamp is expanded as necessary to
      % include target if it is new, but is NOT SUBSTANTIVELY UPDATED
      
      assert(isscalar(labels));
      assert(isscalar(fly)&&isnumeric(fly));
      
      if isempty(labels.flies)
        % AL: special-case branch may be unnecessary; labels.flies should 
        % be col vec even when empty
        tfFlyExists = false; 
      else
        [tfFlyExists,j] = ismember(fly,labels.flies,'rows');
      end
      if ~tfFlyExists
        j = size(labels.flies,1)+1;
      end

      labels.t0s{j} = labelsShort.t0s;
      labels.t1s{j} = labelsShort.t1s;
      labels.names{j} = labelsShort.names;
      labels.flies(j,:) = fly;
      labels.off(j) = labelsShort.off;
      labels.timestamp{j} = labelsShort.timestamp;
      NtimelineTS = numel(labels.timelinetimestamp);
      if NtimelineTS<j
        labels.timelinetimestamp(NtimelineTS+1:j) = {struct()};
      end
      labels.imp_t0s{j} = labelsShort.imp_t0s;
      labels.imp_t1s{j} = labelsShort.imp_t1s;
    end
    
    function labelsShort = labelsShort()
      % labelsShort constructor
      labelsShort = struct('t0s',[],'t1s',[],'names',{{}},'timestamp',[],...
        'off',0,'imp_t0s',[],'imp_t1s',[]);      
    end
    
    function labelsShort = labelsShortFromLabelIdx(labelidx)
      % Construct a labelsShort from a labelidx

      labelsShort = Labels.labelsShort();
      labelsShort.off = labelidx.off;
      
      for iTL = 1:labelidx.nTL
        iBehs = labelidx.TL2idxBeh{iTL};
        assert(isrow(iBehs));
        for iB = iBehs
          [i0s,i1s] = get_interval_ends(labelidx.vals(iTL,:)==iB);
          if ~isempty(i0s)
            n = numel(i0s);
            labelsShort.t0s(end+1:end+n) = i0s - labelidx.off;
            labelsShort.t1s(end+1:end+n) = i1s - labelidx.off;
            labelsShort.names(end+1:end+n) = repmat(labelidx.labelnames(iB),[1,n]);
            labelsShort.timestamp(end+1:end+n) = labelidx.timestamp(iTL,i0s); % first frames of bouts
            assert(all(labelidx.timestamp(iTL,i0s)>0),'Label with missing timestamp.');
          end
        end
      end
      % write importance
      if labelidx.nTL==1
        [i0s,i1s] = get_interval_ends(labelidx.imp);
        if ~isempty(i0s)
          labelsShort.imp_t0s = i0s - labelidx.off;
          labelsShort.imp_t1s = i1s - labelidx.off;
        end
      else
        % ALXXX EXTENDED
        % Multiclassifier importance for GT
      end
    end

    function [labelsShort,tffly] = labelsShortInit(labelsShort,labels,fly)
      % Init labelShort from SCALAR labels and fly specification
      %
      % labels: SCALAR labels, for experiment of interest
      % fly: scalar fly id
      %
      % labelShort: scalar labelsShort structure, for given experiment/fly
      % tffly: true if fly is present in labels and initialization is
      % nontrivial; false otherwise
      
      assert(isscalar(labels));
      assert(isscalar(fly));
          
      [tffly,ifly] = ismember(fly,labels.flies,'rows','legacy');
      if tffly
        labelsShort.t0s = labels.t0s{ifly};
        labelsShort.t1s = labels.t1s{ifly};
        labelsShort.names = labels.names{ifly};
        labelsShort.off = labels.off(ifly);
        if isfield(labels,'imp_t0s')
          labelsShort.imp_t0s = labels.imp_t0s{ifly};
          labelsShort.imp_t1s = labels.imp_t1s{ifly};
        end
        labelsShort.timestamp = labels.timestamp{ifly};
      end
    end
          
  end
  
  methods (Static,Access=private)
    function labels = addBout(labels,iExp,fly,t0,t1,name,off,timestamp)
      % Add a bout to a labels struct array
      % iExp: experiment index (index into labels)
      % fly: fly ID, NOT an index into labels.flies but the absolute ID
      % t0,t1: bout start/stop (scalar doubles)
      % name: bout behavior label (char)
      % off: offset for this fly. Can be [] if fly already exists in labels I guess...
      % timestamp: timestamp for bout
      %
      % WARNING: Does not set labels(iExp).timelinetimestamp(iFly). 
      %
      % For now, all bouts added are added to imp_t0s/imp_t1s as well
      %
      
      assert(iExp<=numel(labels));
      validateattributes(fly,{'numeric'},{'scalar' 'positive' 'integer'});
      
      [labels,iFly] = Labels.addFlyIfNec(labels,iExp,fly,off);
      
      labels(iExp).t0s{iFly}(1,end+1) = t0;
      labels(iExp).t1s{iFly}(1,end+1) = t1;
      labels(iExp).names{iFly}{1,end+1} = name;
      labels(iExp).timestamp{iFly}(1,end+1) = timestamp;
      %labels(iExp).timelinetimestamp{iFly} = struct();
      labels(iExp).imp_t0s{iFly}(1,end+1) = t0;
      labels(iExp).imp_t1s{iFly}(1,end+1) = t1;
    end
    function [labels,iFly] = addFlyIfNec(labels,iExp,fly,off)
      % if fly is not on the list labels(iExp).flies, add it
      %
      % iFly: fly index for iExp/fly, that is labels(iExp).flies(iFly)==fly
      %
      % WARNING: does not initialize timelinetimestamp.
      
      assert(iExp<=numel(labels));
      validateattributes(fly,{'numeric'},{'scalar' 'positive' 'integer'});
      
      tfFly = fly==labels(iExp).flies;
      if ~any(tfFly)
        labels(iExp).flies(end+1,1) = fly;
        labels(iExp).t0s{1,end+1} = zeros(1,0);
        labels(iExp).t1s{1,end+1} = zeros(1,0);
        labels(iExp).names{1,end+1} = cell(1,0);
        labels(iExp).off(1,end+1) = off;
        labels(iExp).timestamp{1,end+1} = zeros(1,0);
        labels(iExp).timelinetimestamp{1,end+1} = struct();
        labels(iExp).imp_t0s{1,end+1} = zeros(1,0);
        labels(iExp).imp_t1s{1,end+1} = zeros(1,0);
        
        iFly = numel(labels(iExp).flies);
      else
        assert(nnz(tfFly)==1,'labels.flies property contains duplicate flies');
        iFly = find(tfFly);

        if ~isempty(off)
          assert(labels(iExp).off(iFly)==off);
        end
      end
    end  
  end
  
  methods (Static) % behavior names, colors
    
    function [nTL,idxBeh2idxTL,tl2IdxBeh] = determineNumTimelines(labelnames)
      % Determine number of timelines to use from labelnames, assuming
      % behavior/label naming convention.
      
      [behs,~,idxBeh2idxTL,tl2IdxBeh] = Labels.verifyBehaviorNames(labelnames);
      nTL = numel(behs);
    end
    
    function [behs,nobehs,idxNames2Beh,tl2IdxNames] = verifyBehaviorNames(names)
      % behs, nobehs: guaranteed to be a (case-sensitive) partitioning of names
      n = numel(names);
      if n==2
        assert(strcmpi(names{2},'none'));
        % 'Classic': names{1} = <beh>, names{2} = 'None';

        behs = names(1);
        nobehs = names(2);
        idxNames2Beh = [1;1];
        tl2IdxNames = {[1 2]};
      elseif n>2 && mod(n,2)==0
        % multibehavior: {'beh1' 'beh2' 'No_beh1' 'No_beh2'};
        nbeh = n/2;        
        behs = cell(1,nbeh);
        nobehs = cell(1,nbeh);
        idxNames2Beh = zeros(2*nbeh,1); 
        tl2IdxNames = cell(nbeh,1);
        for i = 1:nbeh
          iNone = i+nbeh;
          assert(strcmp(Labels.noBehaviorName(names{i}),names{iNone}));
          
          behs{i} = names{i};
          nobehs{i} = names{iNone};
          idxNames2Beh([i iNone]) = i; % ith behavior and its None
          tl2IdxNames{i} = [i iNone];
        end
      else
        assert(false,'Unexpected number of behaviors/labels');
      end
    end
    
    function nobeh = noBehaviorName(beh)
      assert(~strncmp(beh,'No_',3),'Behavior starts with ''No_''.');
      nobeh = ['No_' beh];
    end
    
    function n = noneOrNoBehaviorName(beh,nCls)
      if nCls==1
        n = 'None';
      else
        n = Labels.noBehaviorName(beh);
      end
    end
    
    function lblnames = behnames2labelnames(behnames)
      % Augment behavior names with no-behavior names to get full set of
      % labelnames
      
      lblnames = behnames(:)';
      nbeh = numel(behnames);
      switch nbeh
        case 1 % classic/singleBehavior
          lblnames{1,2} = 'None';
        otherwise % Multiple behaviors
          lblnames = [lblnames cellfun(@Labels.noBehaviorName,lblnames,'uni',0)];
      end
    end
    
    function clr = noneColor(behColor)
      % return standard none color corresponding to behavior color
      % AL: See ShiftColor.m
      assert(numel(behColor)==3);
      clr = behColor/1.5;
    end
    
    function labelcolors = cropOrAugmentLabelColors(labelcolors,nlabels,...
        noBehColorStyle)
      % labelcolors (input): given labelcolors, may be too small or large
      % for nlabels
      % nlabels: number of labels (including no-behavior labels)
      % noBehColorStyle: either 'lines' or 'darkened'. 
      %
      % labelcolors(output): cropped or augmented to have nlabels*3
      % elements. Colors are assumed to correspond to
      % beh1,beh2,...,behN,nobeh1,...nobehN

      assert(mod(nlabels,2)==0,'Number of labels must be even (no-behavior labels should be included)');
      nclrs = numel(labelcolors);
      assert(mod(nclrs,3)==0,'Number of existing labelcolors must be a multiple of 3 (RGB)');

      tfarray = ~isvector(labelcolors);
      if tfarray
        % assumed shape is Nbehx3, cols are R-G-B
        tmp = labelcolors';
        labelcolors = tmp(:);
      end
      
      if nclrs>=3*nlabels
        % truncate existing colors
        labelcolors = labelcolors(1:3*nlabels);
      else
        % augment existing colors
        if nlabels==2
          % Single-classifier
          DEFAULTCOLORS_CLASSIC = [JLabelGUIData.COLOR_BEH_DEFAULT JLabelGUIData.COLOR_NOBEH_DEFAULT];
          labelcolors(nclrs+1:6) = DEFAULTCOLORS_CLASSIC(nclrs+1:6);
          % Don't worry about possible "color collisions" (same color for
          % beh and None)
        else
          % Multi-classifier
          switch noBehColorStyle
            case 'lines'
              COLORS = lines(nlabels)';
              COLORS = COLORS(:);
              labelcolors(nclrs+1:3*nlabels) = COLORS(nclrs+1:3*nlabels);              
            case 'darkened'
              % no-beh colors are darkened version of beh colors
              nrealbeh = nlabels/2;
              COLORS_BEH = lines(nrealbeh)';
              COLORS_BEH = COLORS_BEH(:);
              labelcolors(nclrs+1:nrealbeh*3) = COLORS_BEH(nclrs+1:nrealbeh*3);
              for i = 1:nrealbeh
                idxRealBeh = 3*(i-1)+1:3*(i-1)+3;
                idxNoBeh = idxRealBeh+3*nrealbeh;
                realClr = labelcolors(idxRealBeh);
                noneClr = Labels.noneColor(realClr);
                labelcolors(idxNoBeh) = noneClr;
              end
          end
        end
      end
      
      if tfarray
        labelcolors = reshape(labelcolors,3,[])';
      end
    end
    
    function labelcolors = addNoBehColors(labelcolors)
      assert(isrow(labelcolors));
      assert(mod(numel(labelcolors),3)==0);
      
      nrealbeh = numel(labelcolors)/3;
      for i = 1:nrealbeh
        labelcolors = [labelcolors Labels.noneColor(labelcolors(3*(i-1)+1:3*(i-1)+3))]; %#ok<AGROW>
      end      
    end

    function labelcolors = augmentColors(labelcolors,nlabels,noBehColorStyle)
      % Smarter color-augmentation than cropOrAugment, this does not
      % disturb assignment of existing colors
      
      assert(size(labelcolors,2)==3,'Expected RGB color array.');
      nOrig = size(labelcolors,1);
      assert(mod(nOrig,2)==0,'Expected an even number of colors.');
            
      labelcolors = Labels.cropOrAugmentLabelColors(labelcolors,nlabels,noBehColorStyle);
      labelcolors = [labelcolors(1:nOrig/2,:); labelcolors(end-1,:); ...
                     labelcolors(nOrig/2+1:nOrig,:); labelcolors(end,:)];
    end
      
  end
  
end