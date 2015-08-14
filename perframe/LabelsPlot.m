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
    
    function labels_plot = labelsPlotClearXYLabels(labels_plot,ifly,iBeh,t0,t1)
      % Clear .x, .y for given flies/behaviors/time interval.
      % ifly: vector of fly indices. if -1, use 1:.nflies.
      % iBeh: vector of behavior indices. if -1, use 1:.nbeh.
      % t0, t1: scalar time indices
      
      if isequal(ifly,-1)
        ifly = 1:labels_plot.nflies;
      end
      if isequal(iBeh,-1)
        iBeh = 1:labels_plot.nbeh;
      end
      assert(t0<=t1);
      idx = (t0:t1)+labels_plot.off;
      
      labels_plot.x(:,idx,iBeh,ifly) = nan;
      labels_plot.y(:,idx,iBeh,ifly) = nan;      
    end
    
    function labels_plot = labelsPlotWriteXYLabelsInterval(labels_plot,iFly,...
        iBeh,t0,t1,x,y)
      % Write a track (x,y) to .x, .y over a certain time interval.
      %
      % iFly: scalar fly index
      % iBeh: scalar behavior index
      % t0, t1: scalar time indices
      % x: vector of length t1-t0+2. NOTE: this is ONE LONGER than
      %   numel(t0:t1)
      % y: vector of length t1-t0+2. ONE LONGER than numel(t0:t1).
      
      assert(isscalar(iFly));
      assert(isscalar(iBeh));
      assert(t0<=t1);
      idx = (t0:t1)+labels_plot.off;
      n = t1-t0+1;
      assert(numel(x)==n+1);
      assert(numel(y)==n+1);
      x = x(:)';
      y = y(:)';
      
      labels_plot.x(:,idx,iBeh,iFly) = [x(1:end-1);x(2:end)];
      labels_plot.y(:,idx,iBeh,iFly) = [y(1:end-1);y(2:end)];
    end
      
    function labels_plot = labelsPlotWriteXYFull(labels_plot,ifly,x,y,...
        labelidx,predictedidx,scores,confThreshs)
      % Update labels_plot.x, .y, .predx, .predy from labels/predictions.
      %
      % ifly: scalar fly index (into eg labels_plot.x(:,:,:,ifly)
      % x, y: vectors, track positions
      % labelidx: .vals field from labelIdx structure. Values are in
      %   1:labels_plot.nbeh
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
          labels_plot.predy(1,idxPos,iBehPos,ifly) = y(idxPos);
          labels_plot.predy(2,idxPos,iBehPos,ifly) = y(idxPos1);
          labels_plot.predy(1,idxNeg,iBehNeg,ifly) = y(idxNeg);
          labels_plot.predy(2,idxNeg,iBehNeg,ifly) = y(idxNeg1);          
        end
      end
    end
      
    % Note: labels_plot.x and .y are directly modified in
    % JLabel/SetLabelPlot
    
%     function labels_plot = labelsPlotWriteIm(labels_plot,t0,t1,labelidx,...
%         labelcolors,unkcolor)
%       % Write labels to labels_plot.im in time interval [t0,t1]. Image data
%       % outside that range is unaffected.
%       % 
%       % labelidx: .vals, .imp, should be constrained to time interval [t0 t1]
%       
%       assert(t0<=t1);
%       n = t1-t0+1;
%       assert(labelidx.nbeh==labels_plot.nbeh);
%       assert(isequal(size(labelidx.vals),size(labelidx.imp),...
%         [labels_plot.nTL n]));
%       assert(isequal(size(labelcolors),[labels_plot.nbeh 3]));
%       
%       idx = (t0:t1)+labels_plot.off;
%       imslice = labels_plot.im(:,idx,:);
%       for chan = 1:3
%         imslice(:,:,chan) = unkcolor(chan);
%       end
%       
%       % Permute image matrix to make it easier to assign behavior colors
%       imslice = permute(imslice,[3 1 2]); % imslice is now 3xnTLxn
%       
%       for behaviori = 1:labels_plot.nbeh
%         tf = (labelidx.vals==behaviori) & labelidx.imp;
%         curColor = labelcolors(behaviori,:);
%         for channel = 1:3,
%           imslice(channel,tf) = curColor(channel);
%         end
%         
%         tf = (labelidx.vals==behaviori) & ~labelidx.imp;
%         curColor = ShiftColor.decreaseIntensity(labelcolors(behaviori,:));
%         for channel = 1:3,
%           imslice(channel,tf) = curColor(channel);
%         end
%       end
%       
%       labels_plot.im(:,idx,:) = permute(imslice,[2 3 1]);
% 
%       % JLabel/SetLabelPlot does direct access to labels_plot.im
%     end
    
    function labels_plot = labelsPlotWriteImInterval(labels_plot,iTL,t0,t1,color)
      % Write a single color to a time interval
      % 
      % iTL: vector of timeline indices. If equal to -1, use 1:nTL.
      % t0, t1: time indices
      % color: 3-element RGB color
      
      if isequal(iTL,-1)
        iTL = 1:labels_plot.nTL;
      end
      assert(t0<=t1);
      idx = (t0:t1)+labels_plot.off;
      for chan = 1:3
        labels_plot.im(iTL,idx,chan) = color(chan);
      end      
    end
    
    function labels_plot = labelsPlotWriteIm(labels_plot,labelidx,labelcolors,unkcolor)
      % Write all labels from labelidx to labels_plot.
      
      assert(labels_plot.nbeh==labelidx.nbeh);
      assert(isequal(size(labelcolors),[labels_plot.nbeh 3]));
      
      for chan = 1:3
        labels_plot.im(:,:,chan) = unkcolor(chan);
      end
      
      % Permute image matrix to make it easier to assign behavior colors
      % MK-- permuting maybe slow, removing it.
%       im = permute(labels_plot.im,[3 1 2]); % im is now 3xnTLxn
      
      for behaviori = 1:labels_plot.nbeh
        curTL = labels_plot.beh2TL(behaviori);
        tf = (labelidx.vals(curTL,:) == behaviori) & labelidx.imp(curTL,:);
        curColor = labelcolors(behaviori,:);
        for channel = 1:3,
%           im(channel,tf) = curColor(channel); % if using the permute above.
          labels_plot.im(curTL,tf,channel) = curColor(channel);
        end
        
        tf = (labelidx.vals(curTL,:) == behaviori) & ~labelidx.imp(curTL,:);
        curColor = ShiftColor.decreaseIntensity(labelcolors(behaviori,:));
        for channel = 1:3,
%           im(channel,tf) = curColor(channel);
          labels_plot.im(curTL,tf,channel) = curColor(channel);
        end
      end
      
%       labels_plot.im = permute(im,[2 3 1]);

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
%       im(:) = 0;
      for iBeh = 1:2
        idxScores = pred==iBeh;
        idxPredict = idxScores & abs(scores)>confThresh(iBeh);
        idxBottom = pred_bottom==iBeh;
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
            im(6,idxBottom,channel) = labelcolors(iBeh,channel);
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
      
      scores(isnan(scores)) = 0;

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
