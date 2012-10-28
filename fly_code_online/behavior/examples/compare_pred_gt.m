%% Load data
function OUT = compare_pred_gt(gt_filelist, predSuffix)
%DIR = 'Data/synth_flies/TrainingData/';
%DIR = 'Data/midres_flies/TrainingDataAllMoviesUpdated/';
DIR = 'Data/midres_flies/CleanedMovies/';
%DIR = 'Data/mice/XaviData10behsUseLog/';
%DIR = 'Data/honeybee/rot_and_noise/';
%DIR = 'Data/eric_flies/many_features_and_behs/';
%file containing the list of file that were tested
%gt_filelist = 'traindata_5.txt';
%gt_filelist = 'testdata.txt';
%gt_filelist = 'testdata_special.txt';
%file containing behaviors used for this experiment
beh_file = 'Data/midres_flies/Params/BoyMeetsBoySVMBehaviorParams.txt';
%beh_file = 'Data/mice/Params/MouseBehaviorParams10.txt';
%beh_file = 'Data/honeybee/Params/BehaviorParams.txt';
%beh_file = 'Data/eric_flies/Params/BehaviorParams_more.txt';
%beh_file = 'Data/synth_flies/Params/LungeWalkSVMBehaviorParams.txt';

%repfile = 'reportsUpdated/model5iter17_new';
repfile = 'dummy';
repfile2 = 'dummy2';

%predSuffix = '.label.pred.model53.iter16';
if nargin < 2, gt_filelist = 'traindata_5.txt'; end
if nargin < 2, predSuffix = '.label.pred.model0.iter0'; end


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
max_dur = 0;
%num_tracks = 3;

gt_frames = cell(1,num_tracks);   % frameswise labels per track
for i=1:num_tracks
    filename = [strtok(gt_files{i},'.'),'.label']; % in case file ends with .trx rather than .label
    label = readLabels([DIR filename]);
    duration = label.t1-label.t0+1;
    if duration > max_dur
        max_dur = duration;
    end
    gt_frames{i} = zeros(1,duration);
    for j=1:numel(label.segends)
        for b=1:num_behs
            if strcmp(label.labels{j},behs{b})
                gt_frames{i}(((label.segstarts(j)+1):label.segends(j))) = b-1;
                break
            end
        end
    end 
    %gt_frames{i} = mod(gt_frames{i}-1,10);
end

%% Read predicted boutsc
pred_frames = cell(1,num_tracks);   % frameswise labels per track
invalid = [];
for i=1:num_tracks
    filename = [strtok(gt_files{i},'.'),predSuffix]; % .label.pred.model5.iter17
    try
        label = readLabels([DIR filename]);
    catch
        invalid = [invalid i];
        filename
        continue
    end
    pred_frames{i} = zeros(1,numel(gt_frames{i}));
    for j=1:numel(label.segends)
        for b=1:num_behs
            if strcmp(label.labels{j},behs{b})
                e = min(label.segends(j),numel(gt_frames{i}));
                pred_frames{i}(((label.segstarts(j)+1):e)) = b-1;
                break
            end
        end
    end  
    %pred_frames{i} = mod(pred_frames{i}-1,10);
end

valid = setdiff(1:num_tracks,invalid);
gt_frames = gt_frames(valid);
pred_frames = pred_frames(valid);
num_tracks = numel(valid);
num_invalid = numel(invalid)

%% Plot GT vs Prediction
h1=figure; clf
%colors = {'y','g','r','b','m','c','k'};
colors = jet(num_behs);
% colors(1,:) = [.5 .5 .5];
% colors = [1 0 0;
%           0 0 1;
%           0 1 0];
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
%         if strcmp(behs{b},'none') || strcmp(behs{b},'other')
%             continue
%         end
        % ground truth
        frames = find(gt_frames{i} == b-1);
%         if numel(frames) > 0
%             handles(end+1) = plot(frames,ones(1,numel(frames))*i,['.' colors{b}]);
%         end
        temp = zeros(size(gt_frames{i}));
        temp(frames) = 1;
        cc = bwconncomp(temp);
        for c=1:cc.NumObjects
            frms = cc.PixelIdxList{c};
            if numel(frms) > 0
%                 if b == 1
%                 fill([frms(1) frms(1) frms(end) frms(end)], ...
%                     [i+width i-width i-width i+width ], 'k','FaceColor',[.5 .5 .5], 'EdgeColor','none')
%                 else
                fill([frms(1) frms(1) frms(end)+1 frms(end)+1], ...
                    [i+width i-width i-width i+width ], colors(b,:), 'EdgeColor','none')
%                 end
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
%                 if b == 1
%                 fill([frms(1) frms(1) frms(end) frms(end)], ...
%                     [i+step+width i+step-width i+step-width i+step+width], 'k','FaceColor',[.5 .5 .5], 'EdgeColor','none')   
%                 else
                fill([frms(1) frms(1) frms(end)+1 frms(end)+1], ...
                    [i+step+width i+step-width i+step-width i+step+width], colors(b,:), 'EdgeColor','none')
%                 end
            end
        end
%         plot(frames,ones(1,numel(frames))*i+0.2,['.' colors{b}]);    
        %if numel(frames) > 0
        %    fill([frames(1) frames(1) frames(end) frames(end)], ...
        %        [i+.2+.07 i+.2-.07 i+.2-.07 i+.2+.07 ], colors{b}, 'EdgeColor','none')
        %end
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

%% Plot GT vs Prediction separating ``other" and interesting behaviors
% h2=figure(2); clf
% colors = {'y','g','r','b','m','y','y'};
% handles = [];
% for i=1:num_tracks
%     % Find frame where ``other" behaviors occur    
%     frames_still = find(gt_frames{i} == 0);
%     frames_walk = find(gt_frames{i} == 5);
%     frames_fly = find(gt_frames{i} == 6);
%     frames_other = [frames_still frames_walk frames_fly];
%     frames_other = sort(frames_other);
%     
%     % Find frames where ``interesting" behaviors occur
%     frames_touch = find(gt_frames{i} == 1);
%     frames_lunge = find(gt_frames{i} == 2);
%     frames_wing_th = find(gt_frames{i} == 3);
%     frames_wing_ext = find(gt_frames{i} == 4);
%     frames_nonother = [frames_touch frames_lunge frames_wing_th frames_wing_ext];
%     frames_nonother = sort(frames_nonother);
%     
%     % New frame order
%     frames_reorder = [frames_other frames_nonother];
%     [~, idx] = sort(frames_reorder);
%     
%     for b=1:num_behs
% %         if strcmp(behs{b},'none') || strcmp(behs{b},'other')
% %             continue
% %         end
%         % ground truth
%         frames = find(gt_frames{i} == b-1);
%         if numel(frames) > 0
%             handles(end+1) = plot(idx(frames),ones(1,numel(frames))*i,['.' colors{b}]);
%         end
%         hold on
%         % prediction
%         frames = find(pred_frames{i} == b-1);
%         plot(idx(frames),ones(1,numel(frames))*i+0.2,['.' colors{b}]);        
%     end
%     % plot baselines that span track's duration
%     plot([-2 numel(gt_frames{i})],[i i],'-k')
%     plot([0 numel(pred_frames{i})],[i+0.2 i+0.2],'-','Color',[.5,.5,.5])
% end
% %legend(handles,behs)
% title(['Ground Truth vs SVM prediction'])
% axis([0 numel(gt_frames{1}) 0 i+1.2])
% xlabel('frame number')
% ylabel('GT(lower) and predictions(upper)')
% 
% orient landscape

%% Confusion matrix
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
%% Confusion matrix only main behaviors
% % percent frames correctly labeled
% total_frames = 0;
% true_frames = 0;
% gt_preds = cell(1,num_behs);
% for i=1:num_tracks
%     for j=1:numel(gt_frames{i})
%         if (gt_frames{i}(j) > 4)
%             gt_frames{i}(j)=0; %Set all ``other behaviors to only 1 ``other" behavior
%         end
%         if (pred_frames{i}(j) > 4)
%             pred_frames{i}(j) = 0;
%         end
%     end
% end
%     
% for i=1:num_tracks
%     diff = gt_frames{i}-pred_frames{i};
%     true_frames = true_frames + numel(find(diff==0));
%     total_frames = total_frames + numel(gt_frames{i});
%     for b=1:num_behs-2
%         gt_preds{b} = [gt_preds{b} pred_frames{i}(gt_frames{i}==(b-1))];
%     end
% end
% 
% % confusion matrix for automatic annotator
% confusion_mat = zeros(num_behs-2,num_behs-2);
% for i=1:num_behs-2
%     for j=1:num_behs-2
%         confusion_mat(i,j) = numel(find(gt_preds{i}==(j-1)));
%     end
%     confusion_mat(i,:) = confusion_mat(i,:) / sum(confusion_mat(i,:));
% end
% h5=figure(5); clf;
% imagesc(confusion_mat);
% behs{1}='none';
% behs{4}='wing thrt';
% behs{5}='wing ext';
% colorbar
% set(gca,'XTick',1:num_behs,'XTickLabel',behs,'YTick',1:num_behs,'YTickLabel',behs)
% title(['Confusion matrix (frame based)'])
% 
% orient landscape
% %print('-dpsc','-append',plot_name);
% 
% %% Precision / Recall inrteresting behaviors
% percent_correct_frames = true_frames/total_frames;
% frame_pred_accuracy = mean(diag(confusion_mat));
% 
% overlap_THRESH = .1;
% beh_results = cell(1,num_behs-2);
% for i=1:num_behs-2
%     beh_count = 0;
%     beh_pred_count = 0;
%     false_neg = 0;
%     false_pos = 0;
%     for t=1:num_tracks
%         % search for false negatives
%         temp = zeros(1,numel(gt_frames{t}));
%         temp(gt_frames{t}==(i-1)) = 1;
%         cc_beh = bwconncomp(temp);
%         beh_count = beh_count + cc_beh.NumObjects;
%         % loop through all ground truth bouts
%         for r=1:cc_beh.NumObjects
%             % see how many frames intersect
%             s = cc_beh.PixelIdxList{r}(1);
%             e = cc_beh.PixelIdxList{r}(end);
%             overlap = numel(find(pred_frames{t}(s:e)==(i-1))) / (e-s+1);
%             if overlap < overlap_THRESH
%                 false_neg = false_neg + 1;
%             end
%         end
%         % search for false positives
%         temp = zeros(1,numel(pred_frames{t}));
%         temp(pred_frames{t}==(i-1)) = 1;
%         cc_beh = bwconncomp(temp);        
%         beh_pred_count = beh_pred_count + cc_beh.NumObjects;
%         % loop through all pred bouts
%         for p=1:cc_beh.NumObjects
%             % see how many frames intersect
%             s = cc_beh.PixelIdxList{p}(1);
%             e = cc_beh.PixelIdxList{p}(end);
%             overlap = numel(find(gt_frames{t}(s:e)==(i-1))) / (e-s+1);
%             if overlap < overlap_THRESH
%                 false_pos = false_pos + 1;
%             end
%         end
%     end
%     result.name = behs{i};
%     result.beh_count = beh_count;
%     result.beh_pred_count = beh_pred_count;
%     result.false_neg = false_neg;
%     result.false_pos = false_pos;
%     beh_results{i} = result;
% end
% 
% h6=figure(6); clf
% axis off
% axis ij
% incr = 0.5;
% text(0,0, ['Percent correct frames: ' num2str(percent_correct_frames)]);
% hold on
% text(0,incr, ['Frame pred accuracy: ' num2str(frame_pred_accuracy)]);
% 
% text(0,2.5*incr, 'Precision recall:')
% for i=1:numel(beh_results)
%     text(incr,i*incr + 2.5*incr, [beh_results{i}.name ':   beh count = ' ...
%         num2str(beh_results{i}.beh_count) ', beh pred count = ' ...
%         num2str(beh_results{i}.beh_pred_count) ', false neg = ' ...
%         num2str(beh_results{i}.false_neg) ', false pos = ' ...
%         num2str(beh_results{i}.false_pos)]);
% end
% axis([0 10 0 10])
% %print('-dpsc','-append',plot_name);
% orient landscape
