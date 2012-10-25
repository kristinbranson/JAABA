function visualizer()
    use_trx = 1;
    use_tracking = 1;
    modelId = 2;
    iterId = 17;

    %% Load data
    %DIR = 'Data/midres_flies/TrainingDataAllMoviesUpdated/';
    DIR = 'Data/midres_flies/CleanedMovies/';
    %DIR = 'Data/eric_flies/many_features_and_behs/';
    %DIR = 'Data/honeybee/honeybee_data_mine_new/';
    % file containing the list of file that were tested
    %gt_filelist = ['traindata_' num2str(2) '.txt'];
    %gt_filelist = ['traindata.txt'];
    %gt_filelist = 'train_all.txt';
    gt_filelist = 'traindata_2.txt';
    % file containing behaviors used for this experiment
    beh_file = 'Data/midres_flies/Params/BoyMeetsBoySVMBehaviorParams.txt';
    %beh_file = 'Data/eric_flies/Params/BehaviorParams_more.txt';
    %beh_file = 'Data/honeybee/Params/BehaviorParams.txt';
    
    %% Read behaviors
    lines = textread(beh_file,'%s','delimiter','\n');
    behs = {};
    for i=1:numel(lines)
        if numel(str2num(lines{i}(1))) || strcmp(lines{i}(1),'*')
            continue
        end
        behs{end+1} = strtok(lines{i});
    end
%     behs = {'none','wingext'};
    num_behs = numel(behs);

    %% Read ground truth bouts
    gt_files = textread([DIR gt_filelist],'%s','delimiter','\n');
    max_frame = 0;
    num_tracks = numel(gt_files);
    gt_movie_names = cell(1,num_tracks);
    gt_start_frames = zeros(1,num_tracks);
    gt_fly_ids = zeros(1,num_tracks);
    gt_frames = cell(1,num_tracks);   % frameswise labels per track
    gt_bouts = cell(1,num_tracks);
    for i=1:num_tracks
        filename = [strtok(gt_files{i},'.'),'.label']; % in case file ends with .trx rather than .label
        label = readLabels([DIR filename]);
        gt_movie_names{i} = label.moviename;%(1:20);
        gt_start_frames(i) = label.t0;
        gt_fly_ids(i) = label.flies;        
        % frame-wise ground truth
        gt_frames{i} = zeros(1,label.t1-label.t0+1);
        if label.t1-label.t0+1 > max_frame
            max_frame = label.t1-label.t0+1;
        end
        for b=1:num_behs
            inds = find(strcmp(label.labels,behs{b}));
            for ind=inds
                gt_frames{i}(label.segstarts(ind)+1:label.segends(ind)) = b-1;
            end
        end
        % bout-wise ground truth
        gt_bouts{i} = zeros(numel(label.labels),6);
        for ind=1:numel(label.labels)
            b = find(strcmp(behs,label.labels{ind})) - 1;
            gt_bouts{i}(ind,:) = [b i label.segstarts(ind)+1 label.segends(ind) 0 label.flies];
        end   
    end

    %% Read predicted bouts
    pred_frames = cell(1,num_tracks);   % frameswise labels per track
    pred_bouts = cell(1,num_tracks);
    for i=1:num_tracks
        %filename = [strtok(gt_files{i},'.') '.label.pred.model' num2str(modelId) '.iter' num2str(iterId)]; %.label.interesting  .pred
        filename = [strtok(gt_files{i},'.') '.label.pred.model53.iter18'];%.pred.model1.iter5'
        label = readLabels([DIR filename]);
        % frame-wise ground truth
        pred_frames{i} = zeros(1,label.t1-label.t0+1);
        for b=1:num_behs
            inds = find(strcmp(label.labels,behs{b}));
            for ind=inds
                pred_frames{i}(label.segstarts(ind)+1:label.segends(ind)) = b-1;
            end
        end
        % bout-wise ground truth
        pred_bouts{i} = zeros(numel(label.labels),6);
        for ind=1:numel(label.labels)
            b = find(strcmp(behs,label.labels{ind})) - 1;
            pred_bouts{i}(ind,:) = [b i label.segstarts(ind)+1 label.segends(ind) 0 gt_bouts{i}(1,6)];
        end   
    end
    
    %% Find false negative and false positive bouts
    THRESH_overlap = .2;

    for i=1:num_tracks
        % FALSE NEGATIVES
        % loop through all gt bouts
        for j=1:size(gt_bouts{i},1)
            b_id = gt_bouts{i}(j,1);
            gt_fr = gt_bouts{i}(j,3):gt_bouts{i}(j,4);
            pred_fr = find(pred_frames{i}(gt_fr) == b_id);
            overlap = numel(pred_fr)/numel(gt_fr);
            if overlap > THRESH_overlap
                gt_bouts{i}(j,5) = 1;
            end
        end
        
        % FALSE NEGATIVES
        % loop through all pred bouts
        for j=1:size(pred_bouts{i},1)
            b_id = pred_bouts{i}(j,1);
            pred_fr = pred_bouts{i}(j,3):pred_bouts{i}(j,4);
            gt_fr = find(gt_frames{i}(pred_fr) == b_id);
            overlap = numel(gt_fr)/numel(pred_fr);
            if overlap > THRESH_overlap
                pred_bouts{i}(j,5) = 1;
            end
        end
    end
    
    %% Find gt bouts that are segmented wrong, display movie and params for that
    movies = unique(gt_movie_names);
    vinfos = cell(1,numel(movies));
    features = cell(1,numel(movies));
    if use_tracking
        trks = cell(1,numel(movies));
    end
    for i=1:numel(movies)
        % load video files
        %vinfos{i} = video_open(['/scratch/eeyjolfs/seqs/movies/' movies{i} '.seq']);
        
        vid_dir = '/scratch/eeyjolfs/TRACKER/movies/';
        %vid_dir = '/scratch/eeyjolfs/TRACKER/Eric_scored_movies/';
        tmp = dir([vid_dir movies{i} '/*.seq']); vid_name = tmp.name;
        vinfos{i} = video_open([vid_dir movies{i} '/' vid_name]);
        
        %vid_dir = '/scratch/eeyjolfs/honeybee_movies/';
        %tmp = dir([vid_dir movies{i} '.avi']); vid_name = tmp.name;
        %vinfos{i} = video_open([vid_dir vid_name]);
        
        if use_tracking
            % load trk files
            %trks{i} = load(['/scratch/eeyjolfs/trks_new/' movies{i} '/' movies{i} '-track.mat']);
            trks{i} = load([vid_dir movies{i} '/' vid_name(1:end-4) '-track_1_new.mat']);
        end
        % load feature files
        if use_trx
            features{i}.features = cell(1,2);
%             features{i}.features{1} = zeros(numel(gt_frames{i}),38);
%             features{i}.features{2} = zeros(numel(gt_frames{i}),38);
            if use_tracking
                features{i}.features{1} = zeros(numel(trks{i}.trk.frame_ids),38);
                features{i}.features{2} = zeros(numel(trks{i}.trk.frame_ids),38);
            end
        else
            %features{i} = load(['/scratch/eeyjolfs/seqs/tracks/' movies{i} '_run_features.mat']); 
            features{i} = load([vid_dir movies{i} '/' vid_name(1:end-4) '_features_old.mat']);
        end
    end
    
    if use_trx
%        temp = load('keypoint_fcn'); keypoint_fcn = temp.keypoint_fcn;
        for i=1:num_tracks
            filename = [strtok(gt_files{i},'.'),'.trx']; %.label.interesting
            trx = readTrx([DIR filename]);
            %mov_id = find(strcmp(movies,trx.moviename(1:20)));
            mov_id = find(strcmp(movies,trx.moviename));
            frames = (trx.t0+1):(trx.t1+1);
            if use_tracking
                features{mov_id}.features{trx.flies}(frames,1:size(trx.fields,2)) = trx.fields; 
            else
                features{mov_id}.features{trx.flies} = trx.fields;
            end
%            features{mov_id}.features{trx.flies}(frames,38) = keypoint_fcn{i}';
            feat_names = trx.fieldnames;
%            feat_names{end+1} = 'keypoint_fcn';
        end

        if use_tracking
            for i=1:numel(movies)
                fly1_numframes = numel(trks{i}.trk.sequences{1}.pos.x);
                fly2_numframes = numel(trks{i}.trk.sequences{2}.pos.x);
                features{i}.features{1} = features{i}.features{1}(end-fly1_numframes+1:end,:);
                features{i}.features{2} = features{i}.features{2}(end-fly2_numframes+1:end,:);
            end
        end
    end
    
    %% UI variables
    f_radius = 60;
    if use_trx == 0
        feat_names = {'vel','accel','ang_vel','ang_accel','wing_ang_diff','wing_mean_ang','wing_length',...
                      'dist_food','dist_edge','axis_ratio','axis_ratio_vel','dist_to_other','vel_to_other',...
                      'accel_to_other','ang_between','vel_ang_between','accel_ang_between','facing_angle',...
                      'leg_dist','body_wing_dist'};
%         core_feats = {'vel','ang_vel','min_wing_ang','max_wing_ang','wing_length',...%'dist_food','dist_edge',...
%                   'axis_ratio','dist_to_other','angle_between','facing_angle','leg_dist',...%'body_wing_dist',...
%                   'fg_body_ratio','contrast'};          
%         feat_names = core_feats;
%         for s=1:numel(core_feats)
%             feat_names{end+1} = [core_feats{s} '_gaus1x1'];
% %             feat_names{end+1} = [core_feats{s} '_gaus1x2'];
% %             feat_names{end+1} = [core_feats{s} '_gaus1x4'];
%             feat_names{end+1} = [core_feats{s} '_gaus2x1'];
% %             feat_names{end+1} = [core_feats{s} '_gaus2x2'];
% %             feat_names{end+1} = [core_feats{s} '_gaus2x4'];
%         end
    end
%     colors = [0 0 1;
%               0 1 1;
%               0 1 0;
%               1 1 0;
%               1 0 1;
%               0.7 0.5 1;
%               1 0 0]; 
    if use_tracking
        colors = jet(num_behs);
    else
        colors = [1 0 0;
                  0 0 1;
                  0 1 0];
    end

    curr_beh = 0;  % all behaviors     
    curr_type = 1; % gt intervals
    curr_bout = 1; % current bout under examination    
    do_play = 0;
    feature_id = 1;
    feat_name = feat_names{feature_id};
    bouts = [];
    
    feat_frames = [];
    seg_frames = [];
    trk_id = 1;
    fly_id = 1;
    idx = 1;    
    f = 1;
    f_limit = 1;
    frame_offset = 0;
    start_frame = 1;
    
    %% UI
    % MAIN WINDOW
    scrsz = get(0,'ScreenSize');
    fig_h = figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)]);
    set(fig_h,'Color',[.7 .7 .7])
    %set(fig_h,'WindowStyle','modal')
    %set(fig_h,'MenuBar','none')
    set(fig_h,'Resize','off')
    
    % FEATURE SELECTION LIST
    feat_str = feat_names{1};
    for i=2:numel(feat_names)
        feat_str = [feat_str '|' feat_names{i}];
    end
    uicontrol('Style', 'popup',...
      'String', feat_str,...
      'Position', [scrsz(3)/2-56 148 130 50],...
      'Callback', @set_feature ...
    );       
        
    marginy = -110;
    marginx = scrsz(3)/2 + 70;

    % BEHAVIOR SELECTION LIST
    beh_str = 'all';
    for i=1:num_behs
        beh_str = [beh_str '|' behs{i}];
    end
    uicontrol('Style', 'text', 'String', 'Behavior:', ...
      'Position',[marginx+80 marginy+315 100 23], ...
      'HorizontalAlignment', 'left', ...
      'BackgroundColor', get(fig_h,'Color') ...
    );
    uicontrol('Style', 'popup',...
      'String', beh_str,...
      'Position', [marginx+180 marginy+290 130 50],...
      'Callback', @select_behavior ...
    ); 
    % TYPE SELECTION LIST
    uicontrol('Style', 'text', 'String', 'Bout type:', ...
      'Position',[marginx+80 marginy+280 100 23], ...
      'HorizontalAlignment', 'left', ...
      'BackgroundColor', get(fig_h,'Color') ...
    );    
    type_str = 'ground truth|prediction|false negatives|false positives';
    uicontrol('Style', 'popup',...
      'String', type_str,...
      'Position', [marginx+180 marginy+255 130 50],...
      'Callback', @select_interval_type ...
    ); 
    
    % NUM BOUTS INFORMATION
    num_bouts_h = uicontrol('Style', 'text', ...
      'String', ['Total number of bouts = ' num2str(size(bouts,1))], ...
      'HorizontalAlignment', 'left', ...
      'Position',[marginx+80 marginy+235 250 23], ...
      'BackgroundColor', get(fig_h,'Color') ...
    );
    
    % CLIP NAVIGATION
    uicontrol('Style', 'pushbutton', ...
      'String', '<', ...
      'Callback', @ui_down, ...
      'Position', [marginx+80 marginy+205 20 20] ...
    );   
    curr_bout_h = uicontrol('Style', 'edit', ...
      'String', '1', ... 
      'Callback', @ui_edit_text_box, ...
      'Position', [marginx+105 marginy+203 35 25], ...
      'BackgroundColor', get(fig_h,'Color') ...
    );
    uicontrol('Style', 'pushbutton', ...
      'String', '>', ...
      'Callback', @ui_up, ...
      'Position', [marginx+145 marginy+205 20 20] ...
    );
    
    % PLAY BUTTON     
    butt_h = uicontrol('Style', 'pushbutton', 'String', 'Play',...
      'Position', [marginx+182 marginy+200 123 28],...
      'Callback', @play ...
    );         
          
    %% Run interface
    update_bouts();
    plot_comparison();
    

    %% HELP FUNCTIONS
    function plot_comparison()
        sub_h = subplot('Position',[0.6 0.35 0.3 0.58]);
        hold off
        % plot dots for legend
        for k=1:num_behs
            plot(-10,-10,'.','Color',colors(k,:),'MarkerSize',20)
            hold on
        end
        legend(behs)%,'Location','NorthWestOutside')
        
%         % highlight the current bout being run
%         buff = 20;
%         plot([bouts(curr_bout,3)-buff bouts(curr_bout,4)+buff],[bouts(curr_bout,2) bouts(curr_bout,2)]+0.5,'-k','LineWidth',1);
%         plot([bouts(curr_bout,3)-buff bouts(curr_bout,4)+buff],[bouts(curr_bout,2) bouts(curr_bout,2)]-0.2,'-k','LineWidth',1);
%         plot([bouts(curr_bout,3) bouts(curr_bout,3)]-buff,[bouts(curr_bout,2)-0.2 bouts(curr_bout,2)+0.5],'-k','LineWidth',1);
%         plot([bouts(curr_bout,4) bouts(curr_bout,4)]+buff,[bouts(curr_bout,2)-0.2 bouts(curr_bout,2)+0.5],'-k','LineWidth',1);
        
        % plot the segmentations
        if curr_beh == 0
            plot_behs = 1:num_behs;
        else
            plot_behs = curr_beh;
        end
        for k=1:num_tracks
            % color each frame according to its predicted behavior
            for bb=plot_behs
                indices = find(gt_bouts{k}(:,1) == bb-1);
                for jj=1:numel(indices)
                    ind = indices(jj);
                    plot([gt_bouts{k}(ind,3) gt_bouts{k}(ind,4)],[1 1]*k,'-','Color',colors(bb,:),'LineWidth',4);
                    hold on
                end
                % prediction
                indices = find(pred_bouts{k}(:,1) == bb-1);
                for jj=1:numel(indices)
                    ind = indices(jj);
                    plot([pred_bouts{k}(ind,3) pred_bouts{k}(ind,4)],[1 1]*(k+0.3),'-','Color',colors(bb,:),'LineWidth',4);
                    hold on
                end
            end
            % plot baselines that span track's duration
            plot([0 numel(gt_frames{k})],[k k],'-k')
            plot([0 numel(pred_frames{k})],[k+0.3 k+0.3],'-','Color',[.5,.5,.5])
        end
        title('Ground Truth vs SVM prediction')
        axis([0.5 max_frame 0 k+1.3])
        ylabel('gt  |  prediction')          
        set(sub_h,'YTick',(1:num_tracks)+0.15);
        set(sub_h,'YTickLabel',1:num_tracks);
    end

    function label_comparison()
        % find the positions of the bouts of interest 
        % put some marker on those bouts
    end
    
    function update_bouts()     
        bouts = [];
        inds = [];
        for t=1:num_tracks
            if curr_type == 1       % gt bouts
                if curr_beh == 0
                    bouts = [bouts; gt_bouts{t}];
                else
                    inds = find(gt_bouts{t}(:,1) == curr_beh-1);
                    bouts = [bouts; gt_bouts{t}(inds,:)];
                end
            elseif curr_type == 2   % pred bouts
                if curr_beh == 0
                    bouts = [bouts; pred_bouts{t}];
                else
                    inds = find(pred_bouts{t}(:,1) == curr_beh-1);
                    bouts = [bouts; pred_bouts{t}(inds,:)];
                end
            elseif curr_type == 3   % fn bouts
                fn_inds = find(gt_bouts{t}(:,5) == 0);
                if curr_beh == 0
                    bouts = [bouts; gt_bouts{t}(fn_inds,:)];
                else
                    beh_inds = find(gt_bouts{t}(:,1) == curr_beh-1);
                    inds = intersect(fn_inds,beh_inds);
                    bouts = [bouts; gt_bouts{t}(inds,:)];
                end
            else                    % fp bouts
                fp_inds = find(pred_bouts{t}(:,5) == 0);
                if curr_beh == 0
                    bouts = [bouts; pred_bouts{t}(fp_inds,:)];
                else
                    beh_inds = find(pred_bouts{t}(:,1) == curr_beh-1);
                    inds = intersect(fp_inds,beh_inds);
                    bouts = [bouts; pred_bouts{t}(inds,:)];
                end
            end
        end
     
        update_curr_bout(1);
        if size(bouts,1) == 0
            set(butt_h,'Enable','off')
        else
            set(butt_h,'Enable','on')
            do_play = 0;
            set(butt_h,'String','Play')            
        end   
        set(num_bouts_h,'String',['Total number of bouts = ' num2str(size(bouts,1))]);
    end

    function update_curr_bout(new_curr_bout)
        curr_bout = min(max(new_curr_bout,1),size(bouts,1));
        set(curr_bout_h,'String',num2str(curr_bout));
        update_player();
        %plot_comparison();
    end

    function update_player()
        if curr_bout == 0
            return
        end
        
        trk_id = bouts(curr_bout,2);
        idx = find(strcmp(movies, gt_movie_names{trk_id}));
        fly_id = bouts(curr_bout,6);
        num_this = size(features{idx}.features{mod(fly_id-1,2)+1},1);
        num_other = size(features{idx}.features{mod(fly_id,2)+1},1);        
        frame_offset = 0;
        if num_this <= num_other
            frame_offset = num_other - num_this;
        end        
        start_frame = gt_start_frames(trk_id);
            
        f = bouts(curr_bout,3);
        f_limit = bouts(curr_bout,4);
        
        update_display();
    end

    function update_display()
        % plot the video frame
        sub_h = subplot('Position',[0.1 0.35 0.45 0.58]);
        if start_frame+f == vinfos{idx}.n_frames
            f = f-1;
        end
        img = video_read_frame(vinfos{idx},start_frame+f);
        imagesc(img)
        colormap(gray)
        hold on
        if use_tracking
            % plot the fly itself
            x = trks{idx}.trk.sequences{fly_id}.pos.x(f + start_frame+1-frame_offset);
            y = trks{idx}.trk.sequences{fly_id}.pos.y(f + start_frame+1-frame_offset);
            ori = trks{idx}.trk.sequences{fly_id}.pos.ori(f + start_frame+1-frame_offset);
            wings = trks{idx}.trk.sequences{fly_id}.pos.wings{f + start_frame+1-frame_offset};
            tip = [x+cos(ori)*20 y-sin(ori)*20];
            plot([x tip(1)],[y tip(2)],'-b','LineWidth',2);
            for w=1:numel(wings)
                plot([x wings{w}(1)],[y wings{w}(2)],'-r');
                plot(wings{w}(1),wings{w}(2),'.r', 'MarkerSize', 15)
            end
            plot(x,y,'.','Color',colors(gt_frames{trk_id}(f)+1,:),'MarkerSize',20);
        end
        
        set(sub_h,'XTick',[]);
        set(sub_h,'YTick',[]);   
        title(vinfos{idx}.filename)
        hold off

        % display segmentations
        f_start = max(1,f-f_radius);
        f_end = min(numel(pred_frames{trk_id}),f+f_radius);

        sub_h = subplot('Position',[0.1 0.24 0.45 0.06]);
        b_incl = zeros(1,num_behs);
        handles = [];
        seg_frames = f_start:f_end;
        plot([f f], [-5 5], '-k');
        hold on
        for jj=1:num_behs
            inds = find(gt_frames{trk_id}(seg_frames) == jj-1);
            if numel(inds)
                b_incl(jj) = 1;
                h = plot(f_start+inds-1, -ones(numel(inds),1), '.', 'Color', colors(jj,:),'MarkerSize',20);
                handles = [handles h];
            end
            inds = find(pred_frames{trk_id}(f_start:f_end) == jj-1);
            if numel(inds)
                h = plot(f_start+inds-1, ones(numel(inds),1), '.', 'Color', colors(jj,:),'MarkerSize',20);              
                if b_incl(jj) == 0
                    b_incl(jj) = 1;
                    handles = [handles h];
                end
            end
        end
        hold off
        xlim([f-f_radius f+f_radius]);
        ylim([-5 5]);
        ylabel('gt  |  prediction')
        set(sub_h,'XTick',[]);
        set(sub_h,'YTick',0);

        % display feature         
        sub_h = subplot('Position',[0.1 0.1 0.45 0.1]);
        feat_frames = seg_frames + start_frame-frame_offset; 
        feat_vec = features{idx}.features{fly_id}(feat_frames,feature_id);
        plot(f_start:f_end,feat_vec);
        hold on
        plot([f f], [min(feat_vec)-1 max(feat_vec)+1], '-k')
        hold off
        xlim([f-f_radius f+f_radius]); 
        range = max(feat_vec)-min(feat_vec);
        buff = range*.1;
        ylim([min(feat_vec)-buff max(feat_vec)+buff])
        xlabel('frame id')
        ylabel(feat_name)
        set(sub_h,'XTick',f);
        set(sub_h,'XTickLabel',num2str(start_frame+f));

        if ~ishandle(fig_h)
            return
        end 
        drawnow
        
        if f == f_limit
            do_play = 0;
        end
    end

    %% CALLBACK FUNCTIONS
    function set_feature(hObj,event) %#ok<INUSD>
        val = get(hObj,'Value');
        feature_id = val;        
        % display feature           
        subplot('Position',[0.1 0.1 0.45 0.1])
        feat_vec_ = features{idx}.features{fly_id}(feat_frames,feature_id);
        plot(seg_frames,feat_vec_);
        hold on
        plot([f f], [min(feat_vec_)-1 max(feat_vec_)+1], '-k')
        hold off
        xlim([f-f_radius f+f_radius]);
        range = max(feat_vec_)-min(feat_vec_);
        buff = range*.1;
        ylim([min(feat_vec_)-buff max(feat_vec_)+buff])
        
        rest = feat_names{feature_id};
        [tok,rest] = strtok(rest,'_');
        feat_name = '';
        while (numel(tok))
            feat_name = [feat_name ' ' tok];
            [tok,rest] = strtok(rest,'_');
        end
        ylabel(feat_name)    
    end
    
    function play(hObj,event) %#ok<INUSD>
        if do_play == 1
            do_play = 0;
            set(butt_h,'String','Play')
        elseif do_play == 2
            do_play = 1;
            set(butt_h,'String','Pause')
            update_curr_bout(curr_bout+1);
        else
            do_play = 1;
            set(butt_h,'String','Pause')
        end
        
        while do_play == 1 
           f = f+1;
           update_display();
        end

        if (f == f_limit)
            do_play = 2;
            set(butt_h,'String','Play next')
            if curr_bout == size(bouts,1)
                set(butt_h,'Enable','off')
            end
        end
    end

    function select_behavior(hObj,event) %#ok<INUSD>
        new_curr_beh = get(hObj,'Value')-1;
        if new_curr_beh == curr_beh
            return
        end
        curr_beh = new_curr_beh;
        plot_comparison()
        update_bouts();
    end

    function select_interval_type(hObj,event) %#ok<INUSD>
        new_curr_type = get(hObj,'Value');
        if new_curr_type == curr_type
            return
        end
        curr_type = new_curr_type;
        label_comparison()
        update_bouts();
    end

    function ui_edit_text_box(hObj,event) %#ok<INUSD>
        value = get(hObj,'String');
        update_curr_bout(str2num(value))
        do_play = 0;
        set(butt_h,'String','Play')
        set(butt_h,'Enable','on');
    end

    function ui_down(hObj,event) %#ok<INUSD>
        update_curr_bout(curr_bout-1);
        do_play = 0;
        set(butt_h,'String','Play')
        set(butt_h,'Enable','on');
    end

    function ui_up(hObj,event) %#ok<INUSD>
        update_curr_bout(curr_bout+1);
        do_play = 0;
        set(butt_h,'String','Play')
        set(butt_h,'Enable','on');
    end
end
