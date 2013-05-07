function export_predictions_to_eyruns_format(pred_file_in, behs2, actions_file_pred, actions_file_gt)
    [~, ~, ~, pred_bouts, gt_bouts, gt_movie_names, gt_fly_ids, gt_start_frames]=load_pred_gt(pred_file_in);

    % assumes behs2{1}=='none'
    behs = cell(length(behs2)-1,1);
    for kk=2:length(behs2)
        behs{kk-1} = behs2{kk};
    end
    
    movies = unique(gt_movie_names);
    nflies = max(gt_fly_ids);
    for ii = 1:length(movies),
        bouts_empty = cell(nflies,length(behs));
        for jj=1:nflies,
            for kk=1:length(behs)
                bouts_empty{jj,kk} = zeros(0,3);
            end
        end
        [~,f]=fileparts(movies{ii});
        
        bouts = bouts_empty;
        n = zeros(nflies, length(behs));
        for jj=1:length(pred_bouts)
            id = gt_fly_ids(jj);
            t0 = gt_start_frames(jj);
            for kk=1:size(pred_bouts{jj},1)
                b = pred_bouts{jj}(kk,1);
                if b > 0 && strcmp(gt_movie_names{jj}, movies{ii}),
                   n(id,b) = n(id,b)+1;
                   bouts{id,b}(n(id,b),1) = pred_bouts{jj}(kk,3)+t0;
                   bouts{id,b}(n(id,b),2) = pred_bouts{jj}(kk,4)+t0;
                end
            end
        end
        if ~isempty(actions_file_pred),
            save([actions_file_pred,'_',f,'_actions.mat'], 'behs', 'bouts');
        end
        
        bouts = bouts_empty;
        n = zeros(nflies, length(behs));
        for jj=1:length(gt_bouts)
            id = gt_fly_ids(jj);
            t0 = gt_start_frames(jj);
            for kk=1:size(gt_bouts{jj},1)
                b = gt_bouts{jj}(kk,1);
                if b > 0 && strcmp(gt_movie_names{jj}, movies{ii}),
                   n(id,b) = n(id,b)+1;
                   bouts{id,b}(n(id,b),1) = gt_bouts{jj}(kk,3)+t0;
                   bouts{id,b}(n(id,b),2) = gt_bouts{jj}(kk,4)+t0;
                end
            end
        end
        if ~isempty(actions_file_gt),
            save([actions_file_gt,'_',f,'_actions.mat'], 'behs', 'bouts');
        end
    end
end

