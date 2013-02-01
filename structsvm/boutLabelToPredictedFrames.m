function boutLabelToPredictedFrames(obj, y, expi, flies)
    %ndx = obj.FlyNdxPredict(expi,flies);
        
    for i=1:length(y.bouts)
        b = 1-strcmpi(y.bouts(i).behavior,'none')*2; 
%         curndx = ndx & ...
%             obj.predictdata.t>=y.bouts(i).start_frame & ...
%             obj.predictdata.t<=(y.bouts(i).end_frame-1);
%         obj.predictdata.cur(curndx) = b;
%         obj.predictdata.cur_valid(curndx) = true;
 
        curndx = obj.predictdata{expi}{flies}.t>=y.bouts(i).start_frame & ...
          obj.predictdata{expi}{flies}.t<=(y.bouts(i).end_frame-1);
        obj.predictdata{expi}{flies}.cur(curndx) = b;
        obj.predictdata{expi}{flies}.cur_valid(curndx) = true;
        
    end
end

