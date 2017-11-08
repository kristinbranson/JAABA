classdef RepeatPlotBase < handle
  
  properties
    repeatCacheHandleObj;
  end
  
  methods
    
    function obj = RepeatPlotBase(hobj)
      obj.repeatCacheHandleObj = hobj;
    end
    
    function [ok,Nrpt,tfKeep,datesAllC,gRptC,gCtlExpWashC] = collectData(obj,data)
      % Collect grouping vars: dates, repeats, ctl-exp-wash.
      %
      % tfKeep: numel(data)-by-1 logical, true if trial is included in
      % selected repeats
      % 
      % calls getRepeatsAndSetCache
      
      
      datesAll = {data.date}';
      datesUn = unique(datesAll);
      
      rpts = obj.getRepeatsAndSetCache(datesUn);
      if isempty(rpts)
        ok = false;
        Nrpt = [];
        tfKeep = [];
        datesAllC = [];
        gRptC = [];
        gCtlExpWashC = [];
        return;
      end
      
      Nrpt = numel(rpts);
      datesRpt = [{rpts.ctl}' {rpts.exp}' {rpts.wsh}']; % row i contains repeat i
      tfEmpty = cellfun(@isempty,datesRpt);
      UNMATCHABLE_STR = '__UNMATCHABLE__!@#$%^&*';
      datesRpt(tfEmpty) = {UNMATCHABLE_STR}; % rpts.wsh is optional and can be []. In this case, fill in with random string.
      DATESRPTCOL2LBL = {'ctl' 'exp' 'wsh'};
      
      % restrict to dates in repeats, auxvar/success subset
      tfKeep = ismember(datesAll,datesRpt(:));
      assert(isequal([numel(data) 1],size(datesAll),size(tfKeep)));
      %bstat = bstat(tfKeep,:);
      datesAll = datesAll(tfKeep);
      
      % Form rpt grouping var, ctl/exp grouping var
      gRpt = nan(size(datesAll));
      gCtlExpWash = cell(size(datesAll));
      for iDate = 1:numel(datesAll)
        tf = strcmp(datesAll{iDate},datesRpt);
        assert(nnz(tf)==1);
        [i,j] = find(tf);
        gRpt(iDate) = i;
        gCtlExpWash{iDate} = DATESRPTCOL2LBL{j};
      end
      
      datesTripLong = datesRpt';
      datesTripLong = datesTripLong(:); % c1 e1 w1 c2 e2 w2 ... 
      datesTripLong(strcmp(datesTripLong,UNMATCHABLE_STR),:) = [];
      datesAllC = categorical(datesAll,datesTripLong); 
      gRptC = categorical(gRpt,1:Nrpt);
      gCtlExpWashC = categorical(gCtlExpWash,{'ctl' 'exp' 'wsh'});
      
      ok = true;
    end
    
    function repeats = getRepeatsAndSetCache(obj,datesUn)
      % repeats: if [] or empty, failed to get repeats. Otherwise, nonempty
      % struct array. See Repeat.SelectRepeats.

      % Try to reuse existing repeats; if successful, early return
      if ~isempty(obj.repeatCacheHandleObj.data)
        tmpRpts = obj.repeatCacheHandleObj.data;
        tf = Repeat.uiVerify(tmpRpts,'prestr','Reuse last selected repeats?');
        if tf
          repeats = tmpRpts;
          return;
        end
      end
      
      % lastRepeats was not set, or user does not want to use them
      [repeats,ok] = Repeat.SelectRepeats(datesUn);
      if ok
        obj.repeatCacheHandleObj.data = repeats;
      else
        repeats = [];
      end
    end
  end
  
  methods (Static)
    
    function descStr = FailedRepeatSelection(ax)
      xl = xlim(ax);
      yl = ylim(ax);
      hTxt = text(mean(xl),mean(yl),'Experimental/Control dates unselected. Push Refresh to retry.');
      set(hTxt,'horizontalalignment','center');
      
      descStr = {'Experiment/Control dates unselected'};
    end
        
  end
  
end
