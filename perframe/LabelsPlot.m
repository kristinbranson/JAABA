classdef LabelsPlot 
  
  methods (Static)
    
    function labels_plot = labelsPlot(t0,t1,nTL,nBeh,nflies)      
      n = t1-t0+1;      

      labels_plot = struct();      
      labels_plot.n = n;
      labels_plot.off = 1-t0;
      labels_plot.nbeh = nBeh;
      labels_plot.nflies = nflies;

      labels_plot.im = zeros(nTL,n,3);
      labels_plot.predicted_im = zeros(nTL,n,3); % TODO: multiclass shortcut: just show raw scores, plot gets busy
      labels_plot.suggest_xs = nan;
      labels_plot.error_xs = nan;
      %labels_plot.suggested_im = zeros([1,n,3]);
      %labels_plot.error_im = zeros([1,n,3]);
      labels_plot.x = nan(2,n,nBeh,nflies);
      labels_plot.y = nan(2,n,nBeh,nflies);
      labels_plot.predx = nan(2,n,nBeh,nflies);
      labels_plot.predy = nan(2,n,nBeh,nflies);      
    end
    
    function labels_plot = labelsPlotInitXY(labels_plot,labelidx,ifly,x,y)
      assert(isequal(numel(x),numel(y),labels_plot.n));
      assert(ismember(ifly,1:labels_plot.nflies));
      
      for ibeh = 1:labels_plot.nbeh
        tfBeh = labelidx.vals==ibeh; 
        tfBeh = any(tfBeh,1); % Doesn't matter which timeline(s) the behavior is labeled in
        
        idx = find(tfBeh);
        idx1 = min(idx+1,labels_plot.n);
        labels_plot.x(1,idx,ibeh,ifly) = x(idx);
        labels_plot.x(2,idx,ibeh,ifly) = x(idx1);
        labels_plot.y(1,idx,ibeh,ifly) = y(idx);
        labels_plot.y(2,idx,ibeh,ifly) = y(idx1);
        
        classifierPresent = false;
        if classifierPresent,
          assert(false,'XXXAL');
          % idx = find(predictedidx == behaviori);
          idx = find((predictedidx == ibeh) & ...
            (abs(scores)>data.GetConfidenceThreshold(ibeh)));
          idx1 = min(idx+1,numel(x));
          labels_plot.predx(1,idx,ibeh,ifly) = x(idx);
          labels_plot.predx(2,idx,ibeh,ifly) = x(idx1);
          labels_plot.predy(1,idx,ibeh,ifly) = y(idx);
          labels_plot.predy(2,idx,ibeh,ifly) = y(idx1);
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
    end
      
    function labels_plot = labelsPlotInitPredIm(labels_plot,predidx,labelcolors)
      % labels_plot = labelsPlotInitPredIm(labels_plot,predidx,cmap)
      % predidx: ntimelines x n logical array, true means behavior is predicted.
      % labelcolors: nbehx3 color array
      
      nTL = size(labels_plot.predicted_im,1);
      assert(isequal(size(labels_plot.predicted_im),[nTL labels_plot.n 3]));
      assert(isequal(size(predidx),[nTL labels_plot.n]));
      
      im = labels_plot.predicted_im;
      im(:) = 0;
      for iTL = 1:nTL
        tf = predidx(iTL,:);
        % rows of labelcolors are indexed by behaviors, not timelines,
        % but currently the first nTL behaviors correspond to the
        % timelines.
        curColor = labelcolors(iTL,:);
        for channel = 1:3,
          im(iTL,tf,channel) = curColor(channel);
        end
      end
      
      labels_plot.predicted_im = im;
    end
    
    % JLabel/SetLabelPlot does direct access
    
  end  
  
end