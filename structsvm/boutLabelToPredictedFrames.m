function predictdata = boutLabelToPredictedFrames(y, expi, fly_id)
    predictdata = {};
    predictdata.t = [];
    predictdata.cur = [];
    for i=1:length(y.bouts)
        b = 1-strcmpi(y.bouts(i).behavior,'none')*2; 
        predictdata.t = [predictdata.t y.bouts(i).start_frame:(y.bouts(i).end_frame-1)];
        predictdata.cur = [predictdata.cur, b*ones(1,y.bouts(i).end_frame-y.bouts(i).start_frame)];
    end
    predictdata.cur_valid = ones(1,length(predictdata.t));
    predictdata.exp = expi*ones(1,length(predictdata.t));
    predictdata.flies = fly_id*ones(1,length(predictdata.t));
end

