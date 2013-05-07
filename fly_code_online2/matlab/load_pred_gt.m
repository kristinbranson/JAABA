function [pred_frames, gt_frames, max_dur, pred_bouts, gt_bouts, gt_movie_names, gt_fly_ids, gt_start_frames]=load_pred_gt(pred_file)
    strs = textread(pred_file, '%s', 'bufsize', 100000);
    max_dur = 0;
    pred_frames = {};
    gt_frames = {};
    num_tracks = length(strs);
    gt_movie_names = cell(1,num_tracks);
    gt_start_frames = zeros(1,num_tracks);
    gt_fly_ids = zeros(1,num_tracks);
    gt_bouts = {};
    pred_bouts = {};
    
    for ii=1:length(strs),
        o = loadjson(strs{ii});
        
        t0=o.ground_truth.y.firstframe;
        gt_fly_ids(ii) = o.ground_truth.y.fly_id;
        gt_start_frames(ii) = t0;
        gt_movie_names{ii} = o.ground_truth.y.moviename;
        gt_frames{ii} = [];
        pred_frames{ii} = [];
        gt_bouts{ii} = zeros(length(o.ground_truth.y.bouts),6);
        pred_bouts{ii} = zeros(length(o.predicted.y.bouts),6);
        for jj=1:length(o.ground_truth.y.bouts),
            s = o.ground_truth.y.bouts(jj).start_frame+1-t0;
            e = o.ground_truth.y.bouts(jj).end_frame-t0;
            b = o.ground_truth.y.bouts(jj).behavior;
            gt_frames{ii}(s:e) = b;
            gt_bouts{ii}(jj,:) = [b ii s e 0 gt_fly_ids(ii)];
        end
        for jj=1:length(o.predicted.y.bouts),
            s = o.predicted.y.bouts(jj).start_frame+1-t0;
            e = o.predicted.y.bouts(jj).end_frame-t0;
            b = o.predicted.y.bouts(jj).behavior;
            pred_frames{ii}(s:e) = b;
            pred_bouts{ii}(jj,:) = [b ii s e 0 gt_fly_ids(ii)];
        end
        pred_frames{ii} = pred_frames{ii}(1:length(gt_frames{ii}));
        if length(gt_frames{ii}) > max_dur
            max_dur = length(gt_frames{ii});
        end
    end
end
