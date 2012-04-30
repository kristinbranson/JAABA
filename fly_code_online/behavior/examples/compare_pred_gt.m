
DIR = 'Data/synth_flies/TrainingData/';
%file containing the list of file that were tested
gt_filelist = 'testdata.txt';   
%file containing behaviors used for this experiment
beh_file = 'Data/synth_flies/Params/LungeWalkSVMBehaviorParams.txt';

%% Read behaviors
lines = textread(beh_file,'%s','delimiter','\n');
behs = {};
for i=1:numel(lines)
    if numel(str2num(lines{i}(1))) || strcmp(lines{i}(1),'*')
        continue
    end
    behs{end+1} = strtok(lines{i});
end
num_behs = numel(behs);

%% Read ground truth bouts
gt_files = textread([DIR gt_filelist],'%s','delimiter','\n');
num_tracks = numel(gt_files);

gt_frames = cell(1,num_tracks);   % frameswise labels per track
for i=1:num_tracks
    filename = [strtok(gt_files{i},'.'),'.label']; % in case file ends with .trx rather than .label
    label = readLabels([DIR filename]);
    gt_frames{i} = zeros(1,label.t1);
    for j=1:numel(label.segends)
        for b=1:num_behs
            if strcmp(label.labels{j},behs{b})
                gt_frames{i}(((label.segstarts(j)+1):label.segends(j))) = b-1;
                break
            end
        end
    end       
end

%% Read predicted bouts
pred_frames = cell(1,num_tracks);   % frameswise labels per track
for i=1:num_tracks
    filename = [strtok(gt_files{i},'.'),'.label.pred'];
    label = readLabels([DIR filename]);
    pred_frames{i} = zeros(1,label.t1);
    for j=1:numel(label.segends)
        for b=1:num_behs
            if strcmp(label.labels{j},behs{b})
                pred_frames{i}(((label.segstarts(j)+1):label.segends(j))) = b-1;
                break
            end
        end
    end       
end

%% Plot
figure
colors = {'y','g','r','b','m','c'};
handles = [];
for i=1:num_tracks
    % color each frame according to its predicted behavior
    for b=1:num_behs
%         if strcmp(behs{b},'none') || strcmp(behs{b},'other')
%             continue
%         end
        % ground truth
        frames = find(gt_frames{i} == b-1);
        handles(end+1) = plot(frames,ones(1,numel(frames))*i,['.' colors{b}]);
        hold on
        % prediction
        frames = find(pred_frames{i} == b-1);
        plot(frames,ones(1,numel(frames))*i+0.2,['.' colors{b}]);        
    end
    % plot baselines that span track's duration
    plot([0 numel(gt_frames{i})],[i i],'-k')
    plot([0 numel(pred_frames{i})],[i+0.2 i+0.2],'-','Color',[.5,.5,.5])
end
legend(handles,behs)
title('Ground Truth vs SVM prediction')
axis([0 numel(gt_frames{1}) 0 i+1.2])
xlabel('frame number')
ylabel('GT(lower) and predictions(upper)')

%% Confusion matrix
% percent frames correctly labeled
total_frames = 0;
true_frames = 0;
gt_preds = cell(1,num_behs);
for i=1:num_tracks
    diff = gt_frames{i}-pred_frames{i};
    true_frames = true_frames + numel(find(diff==0));
    total_frames = total_frames + numel(gt_frames{i});
    for b=1:num_behs
        gt_preds{b} = [gt_preds{b} pred_frames{i}(gt_frames{i}==(b-1))];
    end
end

% confusion matrix for automatic annotator
confusion_mat = zeros(num_behs,num_behs);
for i=1:num_behs
    for j=1:num_behs
        confusion_mat(i,j) = numel(find(gt_preds{i}==(j-1)));
    end
    confusion_mat(i,:) = confusion_mat(i,:) / sum(confusion_mat(i,:));
end
figure
imagesc(confusion_mat);
colorbar
set(gca,'XTick',1:num_behs,'XTickLabel',behs,'YTick',1:num_behs,'YTickLabel',behs)
title('Confusion matrix (frame based)')

percent_correct_frames = true_frames/total_frames
frame_pred_accuracy = mean(diag(confusion_mat))

%% Precision / Recall
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
beh_results{:}