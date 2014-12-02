classdef Labels 
  % - Methods for the labelidx data structure (see JLabelData.labelIdx)
  % - Methods for the labels data structure (see JLabelData.labels)
  % - Some stuff related to behavior names
  
  % Label-related datastructures
  % - labelIdx is a scalar struct containing labels for a single track 
  % across multiple timelines. Key props are .vals, .imp, .timestamp which 
  % are numTimelines-by-nsamps.
  % - labelsShort is a scalar struct containing labels for a single track
  % across multiple timelines. Key props are .names, .t0s, .t1s, .imp. It
  % is simlar to the labels structure.
  % - labels is a struct array indexed by experiment and fly, ie
  % labels(iExp){iFly}. Key props are .names, .t0s, .t1s.
  
  methods (Static) % labelidx methods
 
    function labelidx = labelIdx(labelnames,T0,T1)
      % labelIdx 'Constructor' 
      
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
        % ALXXX EXTENDED. Single classifier, leegacy code for writing 
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
    
  end
  
  methods (Static) % labels methods
    
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

      tfModified = false;
      
      Nexp = numel(labels);
      
      if ~isfield(labels,'timelinetimestamp')
        if isempty(labels)
          % trick to add field to empty structure
          [labels(:).timelinetimestamp] = deal([]);
        else
          for iExp = 1:Nexp
            Nfly = numel(labels(iExp).flies);
            labels(iExp).timelinetimestamp = cell(1,Nfly);
          end        
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
    
    function labels = clearLabels(labels,behname,realbehname)
      % Clears all bouts of behavior <behname> corresponding to
      % realbehavior <realbehname>. 
      % Timelinetimestamp is updated
      
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
    
    function labels = renameBehavior(labels,behOld,behNew,realbehOld,realbehNew)
      % behOld: old label
      % realbehOld: "real behavior" corresponding to behOld (for real
      % behavior labels, will be same as behOld)
      % behNew: new label
      % realbehNew: etc
      
      % AL 20140910: This works in the one current callsite, where None is
      % renamed to No_<beh>, but this method is broken for other use cases 
      % b/c you typically need to rename the behavior and its No-behavior 
      % labels simultaneously. (Doing one then the other currently fires 
      % the assert about the timelinetimestamp etc)
            
      for iExp = 1:numel(labels)
        for iFly = 1:numel(labels(iExp).flies)
          % rename bout lbls
          labels(iExp).names{iFly} = regexprep(labels(iExp).names{iFly},behOld,behNew);
          
          % rename timelinetimestamp
          timelineTS = labels(iExp).timelinetimestamp{iFly};
          if strcmp(realbehOld,realbehNew)
            % no operation necessary
          else
            assert(isfield(timelineTS,realbehOld));
            assert(~isfield(timelineTS,realbehNew));
            val = timelineTS.(realbehOld);
            timelineTS.(realbehNew) = val;
            labels(iExp).timelinetimestamp{iFly} = rmfield(timelineTS,realbehOld);
          end
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
    
    function labelsComb = compileLabels(combExpDirNames,labelsCell,expDirNamesCell)
      % Compile/combine multiple label struct arrays.
      % 
      % combExpDirNames: cellstr, desired/target master list of experiments
      % labelsCell: cell array containing labels struct arrays. Typically
      % each element of labelsCell represents labels for a certain
      % behavior.
      % expDirNames: cell array containing cellstrs of expdirnames for
      % labelCell. expDirNames{i} are the expdirnames for labelsCell{i}. 
      % 
      % labelsCell and expDirNamesCell must have the same number of
      % elements; labelsCell{i} and expDirNames{i} must have the same 
      % number of elements.
      %
      % labelsComb: struct array of labels structures labeled by
      % combExpDirNames that contains all bouts in labelsCell.
      
      assert(iscellstr(combExpDirNames));
      assert(iscell(labelsCell));
      assert(iscell(expDirNamesCell));
      assert(numel(labelsCell)==numel(expDirNamesCell));
      
      % Initialize combLabels
      NCombExps = numel(combExpDirNames);
      labelsComb = Labels.labels(NCombExps);
      
      % Loop over input labels and add bouts to combLabels
      N = numel(labelsCell);
      realbehnames = cell(N,1); % realbehnames{i} holds realbehnames for labelsCell{i}
      for i = 1:N
        lbls = labelsCell{i};
        edirs = expDirNamesCell{i};
        assert(numel(lbls)==numel(edirs));
        
        [tf,loc] = ismember(edirs,combExpDirNames);
        assert(all(tf),'Master list of experiments does not cover all labeled experiments.');
        for iExp = 1:numel(lbls)
          iExpComb = loc(iExp);
          for iFly = 1:numel(lbls(iExp).flies)
            fly = lbls(iExp).flies(iFly);
            off = lbls(iExp).off(iFly);
            Nbout = numel(lbls(iExp).t0s{iFly});
            for iBout = 1:Nbout
              labelsComb = Labels.addBout(labelsComb,iExpComb,fly,...
                lbls(iExp).t0s{iFly}(iBout),...
                lbls(iExp).t1s{iFly}(iBout),...
                lbls(iExp).names{iFly}{iBout},...
                off,...
                lbls(iExp).timestamp{iFly}(iBout));
            end            
            
            if iFly==1
              timelineTS = lbls(iExp).timelinetimestamp{iFly};
              realbehnames{i} = fieldnames(timelineTS);
              % fieldnames for .timelinetimestamp should be same for all
              % exps/flies of lbls
            end
          end
        end
      end
      
      realbehnames = cat(1,realbehnames{:});
      assert(numel(realbehnames)==numel(unique(realbehnames)),...
        'Labels to be merged contain duplicate behavior names.');
      
      labelsComb = Labels.initTimelineTimestamps(labelsComb,realbehnames);
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
    
%     function beh = labelName2ClassifierName(name)
%       % eg, return 'Lift' for either 'Lift' or 'No_lift', the point being
%       % that both labels pertain to the Lift classifier.
%       
%       assert(false,'This is not safe, what if name==''None'');
%       if strncmp(name,'No_',3)
%         beh = name(4:end);
%       else
%         beh = name;
%       end
%     end
%     
    function clr = noneColor(behColor)
      % return standard none color corresponding to behavior color
      assert(numel(behColor)==3);
      clr = behColor/1.5;
    end
    
    function labelcolors = cropOrAugmentLabelColors(labelcolors,nRealBeh)
      % labelcolors (input): given labelcolors, may be too small or large for nRealBeh
      % nRealBeh: number of desired real behaviors (not counting None or No_beh)
      %
      % labelcolors(output): cropped or augmented to have 3*nRealBeh
      % elements

      assert(isvector(labelcolors));
      nlblClrs = numel(labelcolors);
      
      if nlblClrs>=3*nRealBeh
        % truncate existing colors
        labelcolors = labelcolors(1:3*nRealBeh);
      else
        % add new colors
        nbehOld = floor(nlblClrs)/3; % nlblClrs should be multiple of 3
        lncolors = lines(nRealBeh);
        newcolors = lncolors(nbehOld+1:nRealBeh,:)';
        newcolors = newcolors(:);
        labelcolors(3*nbehOld+1:3*nRealBeh) = newcolors;
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
    
  end
  
end