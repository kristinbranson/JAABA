classdef JScores < handles
  
  properties (Access = public)
    scores = {};
    flies = [];
    t = [];
    timestamp = [];
    classifierfilenames = [];
  end
  
  methods ( Access = public)
    
    function scores = GetScores(obj,flies)
      scores = obj.scores{flies};
    end
    
    function SetScores(obj,flies,scores)
      obj.scores{flies} = scores;
    end
    
    function scores = NormalizeScores(obj,scores)
      
      if isempty(obj.windowdata.scoreNorm) || isnan(obj.windowdata.scoreNorm)
        isLabeled = obj.windowdata.labelidx_cur~=0;
        wScores = obj.windowdata.scores(isLabeled);
        scoreNorm = prctile(abs(wScores),80);
        obj.windowdata.scoreNorm = scoreNorm;
      end
      
      scoreNorm = obj.windowdata.scoreNorm;
      scores(scores<-scoreNorm) = -scoreNorm;
      scores(scores>scoreNorm) = scoreNorm;
      scores = scores/scoreNorm;
    end
    
    function SaveScores(obj,allScores,expi)
      % Save prediction scores for the whole experiment.
      % The scores are stored as a cell array.
      sfn = obj.GetFile('scores',expi,true);
      obj.SetStatus('Saving scores for experiment %s to %s',obj.expnames{expi},sfn);
      
      didbak = false;
      if exist(sfn,'file'),
        [didbak,msg] = copyfile(sfn,[sfn,'~']);
        if ~didbak,
          warning('Could not create backup of %s: %s',sfn,msg);
        end
      end
      timestamp = obj.classifierTS;
      save(sfn,'allScores','timestamp');
    end
    
    function LoadScores(obj,expi,sfn)
      
      obj.SetStatus('Loading scores for experiment %s from %s',obj.expnames{expi},sfn);
      if ~exist(sfn,'file')
        warndlg('Score file %s does not exist. Not loading scores',sfn);
        return;
      end
      load(sfn,'allScores','timestamp');
      if ~isempty(whos('-file',sfn,'classifierfilename'))
        S = load(sfn,'classifierfilename');
        obj.scoredata.classifierfilenames{expi} = S.classifierfilename;
      else
        obj.scoredata.classifierfilenames{expi} = '';
      end
      
      for ndx = 1:numel(allScores.scores)
        if obj.scoredata.exp
          idxcurr = obj.scoredata.exp == expi & all(bsxfun(@eq,obj.scoredata.flies,ndx),2);
        else
          idxcurr = [];
        end
        if any(idxcurr),
          obj.scoredata.scores(idxcurr) = [];
          obj.scoredata.predicted(idxcurr) = [];
          obj.scoredata.exp(idxcurr,:) = [];
          obj.scoredata.flies(idxcurr,:) = [];
          obj.scoredata.t(idxcurr) = [];
          obj.scoredata.timestamp(idxcurr) = [];
        end
        tStart = allScores.tStart(ndx);
        tEnd = allScores.tEnd(ndx);
        sz = tEnd-tStart+1;
        curScores = allScores.scores{ndx}(tStart:tEnd);
        obj.scoredata.scores(end+1:end+sz) = curScores;
        obj.scoredata.predicted(end+1:end+sz) = -sign(curScores)*0.5+1.5;
        obj.scoredata.exp(end+1:end+sz,1) = expi;
        obj.scoredata.flies(end+1:end+sz,1) = ndx;
        obj.scoredata.t(end+1:end+sz) = tStart:tEnd;
        obj.scoredata.timestamp(end+1:end+sz) = timestamp;
      end
      obj.UpdatePredictedIdx();
      
      if isempty(obj.windowdata.scoreNorm) || isnan(obj.windowdata.scoreNorm)
        if ~isempty(obj.scoredata.scores)
          scoreNorm = prctile(abs(obj.scoredata.scores),80);
          obj.windowdata.scoreNorm = scoreNorm;
        end
      end
      
      obj.ClearStatus();
      
    end
    
  end
  
end