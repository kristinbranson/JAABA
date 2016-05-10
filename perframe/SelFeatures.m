classdef SelFeatures
  
  methods (Static)

    function sf = createEmpty()
      sf = struct('f',[],...
                  'ts',0,... time stamp
                  'count',0,...
                  'npos',0,... number of positive
                  'nneg',0,...
                  'extra1',[],... 
                  'extra2',0,...
                  'extra3',[],...
                  'do',false,... do feature selection at next training
                  'use',false... use feature selection
                  );
    end
    
    function newcls = convertClassifier(cls,sf)
      newcls = cls;
      for ndx = 1:length(cls)
        dim = cls(ndx).dim;
        newdim = sf.f(dim);
        newcls(ndx).dim = newdim;
      end
    end
    
    function newsf = select(sf,data,labels,obj,binVals,bins,params,str)
      if nargin<7,
        str = '';
      end
      
      newsf = sf;
      
      nftrs = params.nselfeatures;
      obj.SetStatus('%s selecting window features to use',str);

      assert(isvector(labels) && numel(labels)==size(data,1));
      assert(all(labels==1 | labels==2));
      assert(nftrs>0);

      boostIterations = nftrs;
      % Learn classifier with all the data.
      posEx = labels == 1;
      negEx = ~posEx;
      numPos = sum(posEx);
      numNeg = sum(negEx);

      if numPos<1 || numNeg<1,
        return;
      end

      modLabels = sign( (labels==1)-0.5);
      wt = getWeights(modLabels);

      [~, allDataModel] = loglossboostLearnRandomFeatures(data,modLabels,boostIterations,wt,binVals,bins,params,[],str);
      newsf.f = unique([allDataModel.dim]);
      newsf.npos = numPos;
      newsf.nneg = numNeg;
      newsf.ts = now();
      newsf.do = false;
      
      
    end
    
    function [newsf,didwarn] = checkForOpt(sf,labels)
      newsf = sf;
      modLabels = sign( (labels==1)-0.5);
      nneg = nnz(modLabels<0);
      npos = nnz(modLabels>0);
      didwarn = false;
      wstr = {'Parameters optimized for training classifiers may',
              'be outdated. Please run optimization again'};
      if abs(sf.npos-npos)/(sf.npos+1) > 0.1 || ...
         abs(sf.nneg-nneg)/(sf.nneg+1) > 0.1 || ...
         sf.count > 20
        didwarn = true;
        warndlg(wstr,'Optimize again..');
      end
      newsf.count = newsf.count+1;
    end
    
    function [sf,mod] = modernize(sf)
      mod = false;
    end
    
  end
  
  
  
end