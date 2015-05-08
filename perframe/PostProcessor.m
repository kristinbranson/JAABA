classdef PostProcessor
  
  methods (Static)    
    
    function posts = PostProcess(curs,postprocessparams,scoreNorm)      
      if isempty(postprocessparams);
        posts = curs;
        return;
      end;
      
      if strcmpi(postprocessparams.method,'hysteresis')
        posts = PostProcessor.ApplyHysteresis(curs,postprocessparams,scoreNorm);
      else
        posts = PostProcessor.ApplyFiltering(curs,postprocessparams);
      end
      
      posts = PostProcessor.RemoveSmallBouts(posts,postprocessparams);      
    end     
    
    % ---------------------------------------------------------------------
    function posts = RemoveSmallBouts(posts,postprocessparams)      
      if postprocessparams.blen > 1 && numel(posts)>0,
        if numel(posts)<= postprocessparams.blen
          if nnz(posts>0) > numel(posts)/2,
            posts(:) = 1;
          else
            posts(:) = -1;
          end
          return;
        end
        
        while true,
          tposts = [posts(:)' 1-posts(end)];
          ends = find(tposts(1:end-1)~=tposts(2:end));
          ends = [1 ends+1];
          blens = ends(2:end)-ends(1:end-1);
          [minblen,smallbout] = min(blens);
          if minblen>=postprocessparams.blen, break; end
          posts(ends(smallbout):ends(smallbout+1)-1) = ...
            1 - posts(ends(smallbout):ends(smallbout+1)-1);
        end
      end
    end
    
    
    % ---------------------------------------------------------------------
    function posts = ApplyHysteresis(curs,params,scoreNorm)
      % Use imfill to find the regions.
      if isempty(curs), posts = curs; return; end
      
      
      % Select pos bouts that have at least one frame about the high
      % threshold.
      hthresh = curs > params.hystopts(1).value*scoreNorm;
      lthresh = curs > 0;
      if( nnz(hthresh)>0)
        pos = imfill(~lthresh,find(hthresh(:))) & lthresh;
        computeNeg = true;
      else
        pos = false(size(curs));
        computeNeg = false;
      end
      % Select neg bouts that have at least one frame below the low
      % threshold.
      hthresh = curs < params.hystopts(2).value*scoreNorm;
      lthresh = curs < params.hystopts(1).value*scoreNorm;
      if nnz(hthresh)>0 && computeNeg,
        neg = imfill(~lthresh,find(hthresh(:))) & lthresh;
      else
        neg = true(size(curs));
      end
      
      posts = pos | ~neg;
      
    end
    
    
    % ---------------------------------------------------------------------
    function posts = ApplyFiltering(curs,params)
      % Use filt to find the regions.
      if isempty(curs), posts = curs; return; end
      filts = conv(curs,ones(1,params.filtopts(1).value),'same');
      posts = filts>0;
    end
    
  end
  
end
  
  
