classdef WindowData
  
  methods (Static)
    
    function wd = windowdata(nCls)
      % Windowdata constructor
      % nCls: number of classifiers
      % wd: struct array, indexed by classifier
      %
      % 
      % .labelidx_cur
      % .labelidx_new
      % .labelidx_old: These three arrays take the same values as
      % Labels.labelidx.vals. A 0 value indicates 'unlabeled'; Otherwise
      % the values are label indices 1:nlabels.
      
      % ALTODO see below, NUA='not used anymore',
      % ALTODO: correct-sized empties (pass windowdataVerify)
      
      wd = struct(...
        'X',single([]), ...
        'exp',cell(nCls,1), ...
        'flies',[], ...
        't',[], ...
        'labelidx_cur',[], ...
        'labelidx_new',[], ...
        'labelidx_old',[], ...
        'labelidx_imp',[], ...
        'featurenames',{{}}, ...
        'predicted',[], ... % NUA
        'predicted_probs',[], ... % NUA
        'isvalidprediction',[], ... % NUA
        'distNdx',[], ... % NUA?
        'scores',[], ... % NUA
        'scoreNorm',[], ... 
        'binVals',[], ... 
        'scores_old',[], ... % NUA
        'scores_validated',[], ... % NUA?
        'postprocessed',[]); % NUA?
    end
    
    function wd = windowdataVerify(wd)
      for iCls = 1:numel(wd)        
        w = wd(iCls);        
        [n,nftr] = size(w.X);
        
        assert(isequal([n 1],size(w.exp),size(w.flies),size(w.t),...
          size(w.labelidx_cur),size(w.labelidx_new),size(w.labelidx_old),...
          size(w.labelidx_imp)) || ...
          isempty(w.exp) && isempty(w.flies) && isempty(w.t) && ... % AL: alternate condition here for legacy compat
          isempty(w.labelidx_cur) && isempty(w.labelidx_new) && ...
          isempty(w.labelidx_old) && isempty(w.labelidx_imp));
        
        assert(numel(w.featurenames)==nftr);

        % scoreNorm?
        % binVals?
      end
    end
    
    function wd = windowdataClear(wd)
      for iCls = 1:numel(wd)
        wd(iCls).X = single([]);
        wd(iCls).exp = [];
        wd(iCls).flies = [];
        wd(iCls).t = [];
        wd(iCls).labelidx_cur = [];
        wd(iCls).labelidx_new = [];
        wd(iCls).labelidx_imp = [];
        wd(iCls).labelidx_old = [];      
        wd(iCls).featurenames = {{}};

        wd(iCls).predicted = [];
        wd(iCls).predicted_probs = [];
        wd(iCls).isvalidprediction = [];
        wd(iCls).distNdx = [];
        wd(iCls).scores = [];
        wd(iCls).scoreNorm = [];
        wd(iCls).binVals = [];
        wd(iCls).scores_old = [];
        wd(iCls).scores_validated = [];
        wd(iCls).postprocessed  = [];
      end
    end
    
    function wd = windowdataTrim(wd,predFcn)
      assert(isa(predFcn,'function_handle'));
      
%       sizeLimit = obj.cacheSize*1e6;
%       classSize = 4;
%       ratioLimit = 0.2;
%       
%       numUnlabeled = nnz(obj.windowdata.labelidx_new==0);
%       numLabeled = nnz(obj.windowdata.labelidx_new);
%       
%       if (nargin < 2 || ~doforce) && (numel(obj.windowdata.X)*classSize < sizeLimit || numUnlabeled/numLabeled<ratioLimit);
%         return;
%       end
%       
%       idx2remove = obj.windowdata.labelidx_new==0 & ...
%         ~obj.FlyNdx(obj.expi,obj.flies);
      
      for i = 1:numel(wd)
        n = size(wd(i).X,1);
        % MAYANK_JAN16_2015: Don't trim if there is no windowdata.
        if n == 0, continue ; end
        tfRm = predFcn(wd(i));
        assert(islogical(tfRm)&&isvector(tfRm)&&numel(tfRm)==n);
        
        wd(i).X(tfRm,:) = [];
        wd(i).exp(tfRm,:) = [];
        wd(i).flies(tfRm,:) = [];
        wd(i).t(tfRm,:) = [];
        wd(i).labelidx_cur(tfRm,:) = [];
        wd(i).labelidx_new(tfRm,:) = [];
        wd(i).labelidx_imp(tfRm,:) = [];
        wd(i).labelidx_old(tfRm,:) = [];
        wd(i).binVals = [];

%         wd(i).distNdx = [];
%         wd(i).isvalidprediction(idxcurr(1:numel(wd(i).isvalidprediction))) = [];
%         wd(i).predicted(idxcurr(1:numel(wd(i).predicted))) = [];
%         wd(i).predicted_probs(idxcurr(1:numel(wd(i).predicted_probs))) = [];
%         wd(i).scores(idxcurr(1:numel(wd(i).scores))) = [];
%         wd(i).scores_old(idxcurr(1:numel(wd(i).scores_old))) = [];
%         wd(i).scores_validated(idxcurr(1:numel(wd(i).scores_validated)),:) = [];
%         wd(i).postprocessed(idxcurr(1:numel(wd(i).postprocessed)),:) = [];
      end
    end
    
    function wd = windowdataSetFeaturenames(wd,featurenames)
      assert(iscell(featurenames)&&numel(featurenames)==numel(wd));
      for i = 1:numel(wd)
        nFtr = numel(featurenames{i});
        if isequal(wd(i).X,[])
          wd(i).X = zeros(0,nFtr);          
        else
          nCol = size(wd(i).X,2);
          assert(nCol==nFtr,...
            'Number of featurenames (%d) does not match number of cols in feature matrix (%d).',...
            nFtr,nCol);
        end
        wd(i).featurenames = featurenames{i};
      end
    end
        
  end
  
end

