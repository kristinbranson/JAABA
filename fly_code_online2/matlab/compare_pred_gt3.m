%% Load data
function OUT = compare_pred_gt3(behs, pred_file)%compare_pred_gt3(model_file, pred_file)

if nargin < 1, behs = {'none','lunge'}; end
if nargin < 2, pred_file = 'data_ssvm_fragments/pred_lunge.txt'; end

%repfile = 'reportsUpdated/model5iter17_new';
repfile = 'dummy';
repfile2 = 'dummy2';

%[params,behs]=load_parameters(model_file); 

[pred_frames, gt_frames, max_dur]=load_pred_gt(pred_file);

plot_pred_gt(pred_frames, gt_frames, max_dur, behs, repfile);

OUT=plot_conf_mat(pred_frames, gt_frames, behs, repfile2);






%% Plot GT vs Prediction
function plot_pred_gt(pred_frames, gt_frames, max_dur, behs, repfile)
    num_tracks = length(gt_frames);
    num_behs = numel(behs);

    h1=figure; clf
    colors = jet(num_behs);
    handles = [];
    for b=1:num_behs
        plot(-10,-10,'.','Color', colors(b,:))
        hold on
    end
    legend(behs)
    step = 0.4;
    width = 0.15;
    for i=1:num_tracks
        % plot baselines that span track's duration
        %plot([0 numel(gt_frames{i})],[i i],'-k')
        %plot([0 numel(pred_frames{i})],[i+step i+step],'-','Color',[.5,.5,.5])
        % color each frame according to its predicted behavior
        for b=1:num_behs
            % ground truth
            frames = find(gt_frames{i} == b-1);
            temp = zeros(size(gt_frames{i}));
            temp(frames) = 1;
            cc = bwconncomp(temp);
            for c=1:cc.NumObjects
                frms = cc.PixelIdxList{c};
                if numel(frms) > 0
                    fill([frms(1) frms(1) frms(end)+1 frms(end)+1], ...
                        [i+width i-width i-width i+width ], colors(b,:), 'EdgeColor','none')
                end
            end
            hold on
            % prediction
            frames = find(pred_frames{i} == b-1);
            temp = zeros(size(pred_frames{i}));
            temp(frames) = 1;
            cc = bwconncomp(temp);
            for c=1:cc.NumObjects
                frms = cc.PixelIdxList{c};
                if numel(frms) > 0
                    fill([frms(1) frms(1) frms(end)+1 frms(end)+1], ...
                        [i+step+width i+step-width i+step-width i+step+width], colors(b,:), 'EdgeColor','none')
                end
            end
        end
    end
    %legend(handles,behs)
    %title(['Ground truth vs predicted segmentation'])
    axis([1 max_dur 0.7 i+0.7])
    %xlabel('frame')
    %ylabel('sequence id')
    % set(gca,'xtick',[])
    % set(gca,'ytick',[])
    % set(gca,'xticklabel',[])
    % set(gca,'yticklabel',[])

    orient landscape
    print('-dpsc',repfile);
end


%% Confusion matrix
function OUT=plot_conf_mat(pred_frames, gt_frames, behs, repfile2)
    num_tracks = length(gt_frames);
    num_behs = numel(behs);
    % percent frames correctly labeled
    total_frames = 0;
    true_frames = 0;
    gt_preds = cell(1,num_behs);
    pred_gts = cell(1,num_behs);
    for i=1:num_tracks
        diff = gt_frames{i}-pred_frames{i};
        true_frames = true_frames + numel(find(diff==0));
        total_frames = total_frames + numel(gt_frames{i});
        for b=1:num_behs
            gt_preds{b} = [gt_preds{b} pred_frames{i}(gt_frames{i}==(b-1))];
            pred_gts{b} = [pred_gts{b} gt_frames{i}(pred_frames{i}==(b-1))];
        end
    end

    % confusion matrix for automatic annotator
    confusion_mat = zeros(num_behs,num_behs);
    for i=1:num_behs
        for j=1:num_behs
            confusion_mat(i,j) = numel(find(gt_preds{i}==(j-1)));
        end
        if sum(confusion_mat(i,:)) == 0
            confusion_mat(i,:) = 0;
        else
            confusion_mat(i,:) = confusion_mat(i,:) / sum(confusion_mat(i,:));
        end
    end

    %valid = 1:10;%[1:5 8 10:13];
    diagonal = diag(confusion_mat)
    %diagonal = temp(valid)
    mean(diagonal)
    save confusion_mat_mice confusion_mat


    confusion_mat_precision = zeros(num_behs,num_behs);
    for i=1:num_behs
        for j=1:num_behs
            confusion_mat_precision(i,j) = numel(find(pred_gts{i}==(j-1)));
        end
        confusion_mat_precision(i,:) = confusion_mat_precision(i,:) / sum(confusion_mat_precision(i,:));
    end

    diagonal = diag(confusion_mat_precision)
    mean(diagonal)

    OUT.confusion_mat_recall = confusion_mat;
    OUT.confusion_mat_precision = confusion_mat_precision;

    %h3=figure(2); clf;
    figure
    subplot(2,2,1)
    imagesc(confusion_mat);
    colorbar
    set(gca,'XTick',1:num_behs,'XTickLabel',behs,'YTick',1:num_behs,'YTickLabel',behs)
    title(['Confusion matrix for ground truth (frame based)'])


    %figure; clf;
    subplot(2,2,2)
    imagesc(confusion_mat_precision)
    colorbar
    set(gca,'XTick',1:num_behs,'XTickLabel',behs,'YTick',1:num_behs,'YTickLabel',behs)
    title(['Confusion matrix for prediction (frame based)'])
    drawnow
    
    %return
    %% Precision / Recall
    percent_correct_frames = true_frames/total_frames;
    diagon = diag(confusion_mat);
    valid = find(~isnan(diagon));
    frame_pred_accuracy = mean(diagon(valid));



    overlap_THRESH = .1;
    beh_results = cell(1,num_behs);
    for i=1:num_behs
        beh_count = 0;
        beh_pred_count = 0;
        false_neg = 0;
        false_pos = 0;
        for t=1:num_tracks
            % search for false negatives
            temp = zeros(1,numel(gt_frames{t}));
            temp(gt_frames{t}==(i-1)) = 1;
            cc_beh = bwconncomp(temp);
            beh_count = beh_count + cc_beh.NumObjects;
            % loop through all ground truth bouts
            for r=1:cc_beh.NumObjects
                % see how many frames intersect
                s = cc_beh.PixelIdxList{r}(1);
                e = cc_beh.PixelIdxList{r}(end);
                overlap = numel(find(pred_frames{t}(s:e)==(i-1))) / (e-s+1);
                if overlap < overlap_THRESH
                    false_neg = false_neg + 1;
                end
            end
            % search for false positives
            temp = zeros(1,numel(pred_frames{t}));
            temp(pred_frames{t}==(i-1)) = 1;
            cc_beh = bwconncomp(temp);        
            beh_pred_count = beh_pred_count + cc_beh.NumObjects;
            % loop through all pred bouts
            for p=1:cc_beh.NumObjects
                % see how many frames intersect
                s = cc_beh.PixelIdxList{p}(1);
                e = cc_beh.PixelIdxList{p}(end);
                overlap = numel(find(gt_frames{t}(s:e)==(i-1))) / (e-s+1);
                if overlap < overlap_THRESH
                    false_pos = false_pos + 1;
                end
            end
        end
        result.name = behs{i};
        result.beh_count = beh_count;
        result.beh_pred_count = beh_pred_count;
        result.false_neg = false_neg;
        result.false_pos = false_pos;
        beh_results{i} = result;
    end
    OUT.beh_results = beh_results;
    %return


    %h4=figure; clf
    subplot(2,2,3:4)
    axis off
    axis ij
    incr = 1.5;
    text(0,0, ['Percent correct frames: ' num2str(percent_correct_frames)]);
    hold on
    text(0,incr, ['Frame pred accuracy: ' num2str(frame_pred_accuracy)]);

    text(0,2.5*incr, 'Precision recall:')
    for i=1:numel(beh_results)
        tp = beh_results{i}.beh_count - beh_results{i}.false_neg;
        pred = tp + beh_results{i}.false_pos;
        text(incr,i*incr + 2.5*incr, [beh_results{i}.name ... 
            ':   beh count = ' num2str(beh_results{i}.beh_count) ...
             ', recall = ' num2str(tp/beh_results{i}.beh_count) ...
             ', precision = ' num2str(tp/pred)]);
    %         num2str(beh_results{i}.beh_count) ', beh pred count = ' ...
    %         num2str(beh_results{i}.beh_pred_count) ', false neg = ' ...
    %         num2str(beh_results{i}.false_neg) ', false pos = ' ...
    %         num2str(beh_results{i}.false_pos)]);
    end
    axis([0 10 0 10])

    orient landscape
    print('-dpsc','-append',repfile2);
end

end
