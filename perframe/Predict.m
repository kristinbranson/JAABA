classdef Predict
  
  methods (Static) % fastPredict
    
    function fastPredict = fastPredict(nClassifier)
      % fastPredict constructor
      % 
      % fastPredict(ibeh).classifier % struct array, reduced/indexed classifier
      % fastPredict(ibeh).windowfeaturescellparams % struct, pffnames -> cell  
      % fastPredict(ibeh).wfs % window features
      % fastPredict(ibeh).pffs % all distinct pffs for wfs; fieldnames of windowfeaturescellparams; 
      % fastPredict(ibeh).ts % classifier TS for this classifier
      % fastPredict(ibeh).wfidx % indices into 'external'/supplied featurename list corresponding to .wfs
      % fastPredict(ibeh).wfidx_valid = false;           
      
      fastPredict = struct(...
        'classifier',cell(nClassifier,1),...
        'windowfeaturescellparams',[],...
        'wfs',{cell(1,0)},...
        'pffs',{cell(1,0)}, ...
        'ts',[], ...
        'wfidx',[], ...
        'wfidx_valid',false);
    end
    
    function fastPredict = fastPredictInit(fastPredict,classifiers,classifierTS,feature_names) %#ok<INUSL>
      % classifiers: cell array of classifiers
      % classifierTS: classifier timestamps
      % feature_names: cellstr
      
        nClassifier = numel(classifiers);
        assert(isequal(nClassifier,numel(classifierTS),numel(feature_names)));
        
        fastPredict = Predict.fastPredict(nClassifier); % No part of fastPredict is reused 
        
        for iCls = 1:nClassifier
          
          if isempty(classifiers{iCls})
            continue;
          end
          
          % which features are actually used
          dims = [classifiers{iCls}.dim];
          feature_names_cls = feature_names{iCls}(dims);
          
          % put these in with the rest of the classifiers' window features
          wfs = {};
          for j = 1:numel(feature_names_cls),
            wfidxcurr = find(WindowFeatureNameCompare(feature_names_cls{j},wfs),1);
            if isempty(wfidxcurr),
              wfidxcurr = numel(wfs)+1;
              wfs{wfidxcurr} = feature_names_cls{j}; %#ok<AGROW>
            end
          end
          
          wf2pff = cellfun(@(x)x{1},wfs,'UniformOutput',false);
          [pffs,~,wf2pffidx] = unique(wf2pff);
          
          windowfeaturescellparams = struct;
          for pfi = 1:numel(pffs),
            pf = pffs{pfi};
            wfidx_cur = wf2pffidx==pfi;
            windowfeaturescellparams.(pf) = WindowFeatureName2Params(wfs(wfidx_cur));
          end
          
          classifiers_indexed = classifiers{iCls};
          for j = 1:numel(classifiers_indexed),
            classifiers_indexed(j).dim = j;
          end
          
          fastPredict(iCls).classifier = classifiers_indexed;
          fastPredict(iCls).windowfeaturescellparams = windowfeaturescellparams;
          fastPredict(iCls).wfs = feature_names_cls;
          fastPredict(iCls).pffs = pffs;
          fastPredict(iCls).ts = classifierTS(iCls);
          fastPredict(iCls).wfidx = [];
          fastPredict(iCls).wfidx_valid = false;
        end
    end
    
    function fastPredict = fastPredictFindWfidx(fastPredict,feature_names)
      feature_names_char = cellfun(@cell2str,feature_names,'UniformOutput',false);
      for i = 1:numel(fastPredict)
        wfs = fastPredict(i).wfs;

        wfs_char = cellfun(@cell2str,wfs,'UniformOutput',false);
        [ism,wfidx] = ismember(wfs_char,feature_names_char);
        if ~all(ism),
          error('Error matching wfs for classifier with window features computed');
        end

  %       wfidx = nan(1,numel(wfs));
  %       for j = 1:numel(wfs),
  %         idxcurr = find(WindowFeatureNameCompare(wfs{j},feature_names));
  %         if numel(idxcurr) ~= 1,
  %           error('Error matching wfs for classifier with window features computed');
  %         end
  %         wfidx(j) = idxcurr;
  %       end

        fastPredict(i).wfidx = wfidx;
        fastPredict(i).wfidx_valid = true;
      end
    end
        
  end
  
  methods (Static) % predictblocks
    
    function pb = predictblocks(nCls)
      pb = struct('t0',cell(nCls,1),'t1',[],'expi',[],'flies',[]);
    end
    
  end
    
end