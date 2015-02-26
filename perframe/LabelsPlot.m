classdef LabelsPlot 
  % timeline images/plots for UI
  
  methods (Static)
    
    function labels_plot = labelsPlot(t0,t1,nTL,nBeh,beh2TL,TL2beh,nflies)
      % labels_plot constructor
      
      n = t1-t0+1;      

      labels_plot = struct();      
      labels_plot.n = n;
      labels_plot.off = 1-t0;
      labels_plot.nflies = nflies;
      labels_plot.nbeh = nBeh;
      labels_plot.nTL = nTL;
      labels_plot.beh2TL = beh2TL;
      labels_plot.TL2beh = TL2beh;

      labels_plot.im = zeros(nTL,n,3);
      labels_plot.predicted_im = []; % image for predicted timeline
      labels_plot.suggest_xs = nan;
      labels_plot.error_xs = nan;
      labels_plot.suggest_gt = nan;
      %labels_plot.suggested_im = zeros([1,n,3]);
      %labels_plot.error_im = zeros([1,n,3]);
      labels_plot.x = nan(2,n,nBeh,nflies);
      labels_plot.y = nan(2,n,nBeh,nflies);
      labels_plot.predx = nan(2,n,nBeh,nflies);
      labels_plot.predy = nan(2,n,nBeh,nflies);      
    end
    
    function labels_plot = labelsPlotInitXY(labels_plot,ifly,x,y,...
        labelidx,predictedidx,scores,confThreshs)
      % ifly: scalar fly index (into eg labels_plot.x(:,:,:,ifly)
      % x, y: vectors, track positions
      % labelidx: .vals field from labelIdx structure. Values are in
      %   1:labels_plot.nbeh      % 
      % predictedidx: optional, can be [] if no prediction. Values are in
      %   {0,1,2,1.5}, see JLabelData.GetPredictedIdx.
      % scores: optional, can be [] if no prediction. See GetPredictedIdx.
      % confThreshs: vector of length labels_plot.nbeh, or [] if no pred
      
      %MERGESTUPDATED
      
      assert(isscalar(ifly) && any(ifly==1:labels_plot.nflies));
      assert(isequal(numel(x),numel(y),labels_plot.n));
      assert(isequal([labels_plot.nTL labels_plot.n],size(labelidx)));
      tfPred = ~isempty(predictedidx);
      if tfPred
        assert(isequal([labels_plot.nTL labels_plot.n],size(predictedidx),size(scores)));
        assert(isvector(confThreshs)&&numel(confThreshs)==labels_plot.nbeh);
      end
      
      for ibeh = 1:labels_plot.nbeh
        iTL = labels_plot.beh2TL(ibeh);
        tfBeh = labelidx(iTL,:)==ibeh; 
        
        idx = find(tfBeh);
        idx1 = min(idx+1,labels_plot.n);
        labels_plot.x(1,idx,ibeh,ifly) = x(idx);
        labels_plot.x(2,idx,ibeh,ifly) = x(idx1);
        labels_plot.y(1,idx,ibeh,ifly) = y(idx);
        labels_plot.y(2,idx,ibeh,ifly) = y(idx1);
      end
      if tfPred
        for iTL = 1:labels_plot.nTL
          iBehs = labels_plot.TL2beh{iTL};
          assert(numel(iBehs)==2);
          iBehPos = iBehs(1);
          iBehNeg = iBehs(2);
          
          idxPos = predictedidx(iTL,:)==1 & ...
                    abs(scores(iTL,:))>confThreshs(iBehPos);
          idxNeg = predictedidx(iTL,:)==2 & ...
                    abs(scores(iTL,:))>confThreshs(iBehNeg);            
          idxPos = find(idxPos);
          idxNeg = find(idxNeg);
          idxPos1 = min(idxPos+1,labels_plot.n);
          idxNeg1 = min(idxNeg+1,labels_plot.n);
          
          labels_plot.predx(1,idxPos,iBehPos,ifly) = x(idxPos);
          labels_plot.predx(2,idxPos,iBehPos,ifly) = x(idxPos1);
          labels_plot.predx(1,idxNeg,iBehNeg,ifly) = x(idxNeg);
          labels_plot.predx(2,idxNeg,iBehNeg,ifly) = x(idxNeg1);
          labels_plot.predy(1,idxPos,iBehPos,ifly) = x(idxPos);
          labels_plot.predy(2,idxPos,iBehPos,ifly) = x(idxPos1);
          labels_plot.predy(1,idxNeg,iBehNeg,ifly) = x(idxNeg);
          labels_plot.predy(2,idxNeg,iBehNeg,ifly) = x(idxNeg1);          
        end
      end
    end
      
    function labels_plot = labelsPlotInitIm(labels_plot,labelidx,labelcolors)
      assert(labels_plot.nbeh==labelidx.nbeh);
      assert(isequal(size(labelcolors),[labels_plot.nbeh 3]));
      
      labels_plot.im(:) = 0;
      
      % Permute image matrix to make it easier to assign behavior colors
      im = permute(labels_plot.im,[3 1 2]); % im is now 3xnTLxn
      
      for behaviori = 1:labels_plot.nbeh
        tf = (labelidx.vals == behaviori) & labelidx.imp;
        curColor = labelcolors(behaviori,:);
        for channel = 1:3,
          im(channel,tf) = curColor(channel);
        end
        
        tf = (labelidx.vals == behaviori) & ~labelidx.imp;
        curColor = ShiftColor.decreaseIntensity(labelcolors(behaviori,:));
        for channel = 1:3,
          im(channel,tf) = curColor(channel);
        end
      end
      
      labels_plot.im = permute(im,[2 3 1]);

      % JLabel/SetLabelPlot does direct access to labels_plot.im
    end
      
    function labels_plot = labelsPlotSetPredIm(labels_plot,...
        scores,pred,confThresh,bottomType,bottomDist,...
        scores_bottom,pred_bottom,labelcolors,scorecolor)
      % Set labels_plot.predicted_im, single classifier version
      % - Top: binary prediction
      % - Middle: raw scores
      % - Bottom: configureable
      %
      % scores: 1 x n double, normalized scores
      % pred: 1 x n. 1==beh, 2==no-beh (also allowed, 0 or 1.5. See JLabelData.GetPredictedIdx)
      % confThresh: 1x2 array for beh/noBeh resp
      % bottomType: char enum
      % bottomDist: n-long vector, only used if bottomType=='Distance'
      % scores_bottom: 1 x n double
      % pred_bottom: 1 x n. 1==beh, 2==no-beh (not used if bottomType=='Distance')
      % labelcolors: nbeh x 3 (nbeh==2)
      % scorecolors: 63x3x3
      
      %MERGESTUPDATED
      
      assert(isequal(size(scores),size(pred),size(scores_bottom),size(pred_bottom)));
      assert(numel(confThresh)==2);
      %assert(size(labels_plot.predicted_im,1)==6);
      
      labels_plot.predicted_im = zeros(6,labels_plot.n,3);

      idxBottomScores = ~isnan(scores_bottom);
      bottomScoreNdx = ceil(scores_bottom(idxBottomScores)*31)+32;      

      im = labels_plot.predicted_im;
      im(:) = 0;
      for iBeh = 1:2
        idxScores = pred==iBeh;
        idxPredict = idxScores & abs(scores)>confThresh(iBeh);
        scoreNdx = ceil(scores(idxScores)*31)+32; % scoreNdx in [1,63]
        for channel = 1:3
          im(1,idxPredict,channel) = labelcolors(iBeh,channel);
          im(2,idxPredict,channel) = labelcolors(iBeh,channel);
          im(3,idxScores,channel) = scorecolor(scoreNdx,channel,1);
          im(4,idxScores,channel) = scorecolor(scoreNdx,channel,1);
          
          % Bottom row (configureable)
          if strcmp(bottomType,'Distance')
            im(5:6,:,channel) = repmat(1-bottomDist(:)',[2 1 1]);
            im(5:6,isnan(bottomDist(:)'),channel) = 0;
          else
            im(5,idxBottomScores,channel) = scorecolor(bottomScoreNdx,channel,1);
            im(6,pred_bottom==iBeh,channel) = labelcolors(iBeh,channel);
          end
        end        
      end
      
      labels_plot.predicted_im = im;      
    end
    
    function labels_plot = labelsPlotSetPredImMultiClsAnalog(labels_plot,scores,scorecolors)
      % scores: 1 x n double, normalized scores
      % scorecolors: cell vector with nclassifiers elements. Each element 
      %   scorecolors{iCls} is 63x3x3
      
      nTL = labels_plot.nTL;
      assert(isequal(size(scores),[nTL labels_plot.n]));
      assert(numel(scorecolors)==nTL);

      labels_plot.predicted_im = zeros(nTL,labels_plot.n,3);
      im = labels_plot.predicted_im;
      im(:) = 0;
      for iTL = 1:nTL
        scoreNdx = ceil(scores(iTL,:)*31)+32; % scoreNdx in [1,63]
        scoreclr = scorecolors{iTL};
        for channel = 1:3
          im(iTL,:,channel) = scoreclr(scoreNdx,channel,1);
        end
      end
      
      labels_plot.predicted_im = im;
    end

    function labels_plot = labelsPlotSetPredImMultiClsBinary(labels_plot,...
        predIdx,scores,confThreshs,labelcolors)
      % Set labels_plot.predicted_im, multiclass version
      % At the moment, we just show the predicted (binary) bouts, as the 
      % regular tripartite timeline is likely too cluttered for multiple
      % classifiers.
      %
      % predIdx: ntimelines x nSamp array, values in {0,1,2,1.5}, see JLabelData.GetPredictedIdx
      % scores: ntimelines x nSamp
      % confThreshs: ntimelines x 2. confThreshs(i,:) is [lo hi] thresholds
      %   for ith classifier/timeline 
      % labelcolors: nbehaviors x 3 color array
      
      %MERGESTUPDATED
      
      nTL = labels_plot.nTL;
      nSamp = labels_plot.n;
      assert(isequal([nTL nSamp],size(predIdx),size(scores)));
      assert(isequal(size(confThreshs),[nTL 2]));
      assert(isequal(size(labelcolors),[labels_plot.nbeh 3]));
      
      labels_plot.predicted_im = zeros(nTL,nSamp,3);
      %assert(isequal(size(labels_plot.predicted_im),[nTL labels_plot.n 3]));
      
      im = labels_plot.predicted_im;
      im(:) = 0;
      for iTL = 1:nTL
        iLbls = labels_plot.TL2beh{iTL};
        for i12 = [1 2]
          iBeh = iLbls(i12);
          ct = confThreshs(iBeh);
          assert(ct==confThreshs(iTL,i12));
          
          tfPred = predIdx(iTL,:)==i12 & abs(scores(iTL,:))>ct;
          
          curColor = labelcolors(iBeh,:);
          for channel = 1:3,
            im(iTL,tfPred,channel) = curColor(channel);
          end
        end
      end
      
      labels_plot.predicted_im = im;
    end
    
  end
  
end
