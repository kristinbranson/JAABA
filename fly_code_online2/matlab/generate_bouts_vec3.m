
function [trainvals, boutfeat] = ...
                 generate_bouts_vec(Xs,Bouts,Behs,first_valids,savenames,FPS,p,normalize,trainvals)
    % Usage:
    % boutfeat = generate_bouts(X,bouts,[p],[normalize],[feat_mu],[feat_sigma])
    %
    % X:         n_frames x n_features matrix
    % bouts:     n_bouts x 2 (start,end) matrix
    % FPS:       framerate (used to determine size of haar features)
    % p:         parameters for determining which bout features are computed
    % normalize: normalize bout featues (default 0)
    % trainvals: bout parameters gathered from training data set
    %   trainvals.histogram_thresholds
    %   trainvals.ave_feature_responses
    %   trainvals.min_feature_responses
    %   trainvals.max_feature_responses
    %   trainvals.
    %   trainvals.feat_mu
    %   trainvals.feat_gamma
    
    if nargin < 7 || isempty(p)
        p.add_wavelets                                = false;
        p.num_histogram_bins                          = 8;
        p.num_temporal_levels                         = 3;
        p.use_histogram_ave_features                  = true;
        p.use_histogram_sum_features                  = false;
        p.use_bout_ave_absolute_features              = true;
        p.use_bout_ave_features                       = true;
        p.use_bout_max_feature                        = true;
        p.use_bout_min_feature                        = true;
        p.use_bout_sum_absolute_features              = false;
        p.use_bout_sum_features                       = false;
        p.use_bout_standard_deviation                 = true;
        p.use_bout_sum_variance                       = false;
        p.use_ave_absolute_harmonic_features          = true;
        p.use_ave_harmonic_features                   = true;
        p.use_sum_absolute_harmonic_features          = false;
        p.use_sum_harmonic_features                   = false;
        p.use_bout_change                             = true;        
        p.use_bout_absolute_change                    = true;
        p.use_start_ave_absolute_diff_haar_features   = true;
        p.use_start_ave_diff_haar_features            = true;
        p.use_start_sum_absolute_diff_haar_features   = false;
        p.use_start_sum_diff_haar_features            = false;        
        p.use_end_ave_absolute_diff_haar_features     = true;
        p.use_end_ave_diff_haar_features              = true;
        p.use_end_sum_absolute_diff_haar_features     = false;
        p.use_end_sum_diff_haar_features              = false;
        p.use_global_difference_ave_ave_features      = false;
        p.use_global_difference_ave_sum_features      = false;
        p.use_global_difference_max_ave_features      = false;
        p.use_global_difference_max_sum_features      = false;
        p.use_global_difference_min_ave_features      = false;
        p.use_global_difference_min_sum_features      = false; 
        p.use_duration_feature                        = true;        
        
        
%         p.num_histogram_bins                          = 8;
%         p.num_temporal_levels                         = 3;
%         p.use_histogram_ave_features                  = true;
%         p.use_histogram_sum_features                  = false;
%         p.use_bout_ave_absolute_features              = false;
%         p.use_bout_ave_features                       = true;
%         p.use_bout_max_feature                        = true;
%         p.use_bout_min_feature                        = true;
%         p.use_bout_sum_absolute_features              = false;
%         p.use_bout_sum_features                       = false;
%         p.use_bout_standard_deviation                 = true;
%         p.use_bout_sum_variance                       = false;
%         p.use_ave_absolute_harmonic_features          = false;
%         p.use_ave_harmonic_features                   = false;
%         p.use_sum_absolute_harmonic_features          = false;
%         p.use_sum_harmonic_features                   = false;
%         p.use_bout_change                             = false;        
%         p.use_bout_absolute_change                    = false;
%         p.use_start_ave_absolute_diff_haar_features   = false;
%         p.use_start_ave_diff_haar_features            = false;
%         p.use_start_sum_absolute_diff_haar_features   = false;
%         p.use_start_sum_diff_haar_features            = false;        
%         p.use_end_ave_absolute_diff_haar_features     = false;
%         p.use_end_ave_diff_haar_features              = false;
%         p.use_end_sum_absolute_diff_haar_features     = false;
%         p.use_end_sum_diff_haar_features              = false;
%         p.use_global_difference_ave_ave_features      = false;
%         p.use_global_difference_ave_sum_features      = false;
%         p.use_global_difference_max_ave_features      = false;
%         p.use_global_difference_max_sum_features      = false;
%         p.use_global_difference_min_ave_features      = false;
%         p.use_global_difference_min_sum_features      = false; 
%         p.use_duration_feature                        = true;
    end

    % compute n_bout_feat
    n_bout_feat = get_num_feat(p);
    n_base_feat = size(Xs{1},2);
    n_total_feat = n_bout_feat*n_base_feat+p.use_duration_feature;
    td = max(round(1/30*FPS),1); % time delta used for haar like features

    % initialize temporal levels
    temporal_levels = 1:p.num_temporal_levels;
    temporal_grid_size = zeros(size(temporal_levels));    

    % apply wavelets to each feature vector
    if p.add_wavelets
        if FPS == 30
            wavelets(1,:) = [0 0.5 0 -0.5 0];
            wavelets(2,:) = [0 0.25 -0.5 0.25 0];
        elseif FPS == 200
            h = fspecial('gaussian', [1 34],4.3);
            wavelets(1,:) = diff(h);
            h = fspecial('gaussian', [1 35],4.3);
            wavelets(2,:) = diff(diff(h));
        else
            disp(['No wavelet option implemented for FPS = ' num2str(FPS)])
            return
        end
        for jj=1:numel(Xs)            
            X = Xs{jj};
            X_diff  = zeros(size(X));
            X_ddiff = zeros(size(X));
            first_valid_fr = first_valids(jj);
            n_frames = size(X,1);
            for ii=1:n_base_feat
                responses = conv(X(:,ii),wavelets(1,:),'valid');
                buff_left = floor((n_frames-numel(responses))/2);
                buff_right = ceil((n_frames-numel(responses))/2);
                X_diff(first_valid_fr:first_valid_fr+buff_left-1,ii) = responses(1);
                X_diff(first_valid_fr+buff_left:end-buff_right,ii) = responses;
                X_diff(end-buff_right+1:end,ii) = responses(end);
                responses = conv(X(:,ii),wavelets(2,:),'valid');
                buff_left = floor((n_frames-numel(responses))/2);
                buff_right = ceil((n_frames-numel(responses))/2);
                X_ddiff(first_valid_fr:first_valid_fr+buff_left-1,ii) = responses(1);
                X_ddiff(first_valid_fr+buff_left:end-buff_right,ii) = responses;
                X_ddiff(end-buff_right+1:end,ii) = responses(end);
            end
            Xs{jj} = [X X_diff X_ddiff];
        end
    end
    
    if nargin < 9 || isempty(trainvals)
        disp('computing histogram thresholds') 
        X = [];
        for jj=1:numel(Xs);
            X = [X; Xs{jj}];
        end
        tic
        % compute histogram thresholds
        histogram_thresholds = zeros(n_base_feat,p.num_histogram_bins-1);
        for ii=1:n_base_feat
            vals = X(:,ii);
            vals = sort(vals(~isnan(vals)));
            chunk_hist = numel(vals)/(p.num_histogram_bins);
            % - set histogram thresholds to be such that there is an equal number over
            %     values in between all thresholds
            for jj=1:p.num_histogram_bins-1
                target = jj*chunk_hist;
                ind = floor(target);
                w = 1-abs(target-ind);
                % if histogram thresholds same as previous one, use the next
                % value up
                while jj>1 && vals(ind) == histogram_thresholds(ii,jj-1)
                    ind = ind+1;
                    w = 1;
                end
                histogram_thresholds(ii,jj) = vals(min(ind,numel(vals)))*w + ...
                                              vals(min(ind+1,numel(vals)))*(1-w);
            end        
        end
        toc
        % compute ave_feature_responses, min_vals, and max_vals
        ave_feature_responses = nanmean(X);
        min_feature_responses = nanmin(X);
        max_feature_responses = nanmax(X);
    else
        histogram_thresholds = trainvals.histogram_thresholds;
        ave_feature_responses = trainvals.ave_feature_responses;
        min_feature_responses = trainvals.min_feature_responses;
        max_feature_responses = trainvals.max_feature_responses;
    end    
    
    feat_mus = cell(1,numel(Xs));
    for jj=1:numel(Xs)
        X = Xs{jj};
        bouts = Bouts{jj};
        behs = Behs{jj};
        T = size(X,1);        
        first_valid_fr = first_valids(jj);

        disp('computing integral features')
        tic
        % compute integral features
        integral_features = zeros(size(X,1)+1,n_base_feat);
        integral_features(first_valid_fr+1:end,:) = cumsum(X(first_valid_fr:end,:));

        % compute integral squared features
        %integral_sqr_features = zeros(size(X,1)+1,n_base_feat);
        %integral_sqr_features(first_valid_fr+1:end,:) = cumsum(X(first_valid_fr:end,:).^2);
        toc


        disp('computing histogram bins')
        %compute histogram bins
    % 	// Compute histogram features.  We discretize the space of values for the i_th feature into num_histogram bins
    % 	// of equal size between min_feature_responses[i] and max_feature_responses[i].  histogram_bins[i][j] stores the
    % 	// index of the assigned bin for the i_th feature at the j_th time step
        histogram_bins = zeros(T,n_base_feat);
        tic
        for ii=1:n_base_feat
            f_vals = X(:,ii);
            histogram_bins(:,ii) = p.num_histogram_bins;
            for kk=p.num_histogram_bins-1:-1:1;
                histogram_bins(f_vals<=histogram_thresholds(ii,kk),ii) = kk;
            end
        end    
        toc

        disp('computing integral histogram features')
        tic
        % compute integral histogram features
    % 	// Compute integral features on the histogram features (e.g. using the same method as for integral_features),
    % 	// such that the total count of the j_th histogram bin for the i_th feature can be computed using
    % 	//   integral_histogram_features[i][j][e]-integral_histogram_features[i][j][s]
        integral_histogram_features = zeros(T+1,p.num_histogram_bins,n_base_feat);
        tic
        for kk=1:p.num_histogram_bins
            integral_histogram_features(:,kk,:) = 0;
            tmp = histogram_bins == kk;
            integral_histogram_features(first_valid_fr+1:end,kk,:) = cumsum(tmp(first_valid_fr:end,:));
        end
        toc

        disp(['computing features for ' num2str(size(bouts,1)) 'bouts'])
        tic
        % compute features for all bouts
        boutfeat = zeros(size(bouts,1),n_total_feat);
        for b=1:size(bouts,1)
            if mod(b,10000)==0
                disp([num2str(b/size(bouts,1)*100) ' %'])
            end
            %boutfeat(b,:) = compute_bout_feat(bouts(b,1),bouts(b,2));
            t_start = bouts(b,1);
            t_end = bouts(b,2);

            feat = zeros(1,n_total_feat);
            dur = t_end-t_start+1;
            if dur < 1
                disp('Error: t_start must be smaller than t_end')
                return
            end
            inv = 1/dur;

            %% compute temporal grid size and harmonic grid size
            for i=1:p.num_temporal_levels
                temporal_grid_size(i) = dur/temporal_levels(i);
            end

            %% loop through base features
            ind = 0;
            %% compute sum, ave, std, min and max over temporal regions    
            %     // Sum feature response in the interval (t_start,t_end).  Can also divide into multiple temporal
            %     // levels, and compute the average response in grids of size (t_end-t_start)/2, (t_end-t_start)/4, 
            %     // (t_end-t_start)/8, ... This is an average feature response in temporal intervals of different size,
            %     // where the outputted feature vector is the concatenation of features over each temporal interval
            %     // The output feature space has size 
            %     //   F_NUM_GOOD_GLOBAL_FEATURES*(sum_{l=0}^num_temporal_levels(2^i))            
            for j=p.num_temporal_levels:-1:1
                bout = X(t_start:t_end,:);

                f_harm = 0;
                start = 1;
                for k=1:temporal_levels(j)
                    temp_dur = temporal_grid_size(j);
                    temp_inv = 1/temp_dur;
                    finish = k*temporal_grid_size(j);
                    
                    % compute sum, mean, min, max, variance
                    
                    % first restrict consideration to frames fully within
                    % the applicable temporal region
                    f_sum = sum(bout(ceil(start):floor(finish),:));
                    sum_sqr = sum(bout(ceil(start):floor(finish),:).^2);
                    f_min = min(bout(ceil(start):floor(finish),:));
                    f_max = max(bout(ceil(start):floor(finish),:));
                    
                    w = start-floor(start);
                    if w > 0,
                        % Add the applicable fraction of a frame at the
                        % start of the region into sums, while for min/max, 
                        % interpolate to find feature values at the bout
                        % start
                        frac = (1-w)*(bout(floor(start),:));
                        frac_sqr = (1-w)*(bout(floor(start),:).^2);
                        interp = frac+w*bout(ceil(start),:);
                        f_sum = f_sum + frac;
                        sum_sqr = sum_sqr + frac_sqr;
                        f_min = min(f_min,interp);
                        f_max = max(f_max,interp);
                    end
                    w = finish-floor(finish);
                    if w > 0,
                        % Add the applicable fraction of a frame at the
                        % end of the region into sums, while for min/max, 
                        % interpolate to find feature values at the bout
                        % end
                        frac = w*(bout(ceil(finish),:));
                        frac_sqr = w*(bout(ceil(finish),:).^2);
                        interp = frac+(1-w)*bout(floor(finish),:);
                        f_sum = f_sum + frac;
                        sum_sqr = sum_sqr + frac_sqr;
                        f_min = min(f_min,interp);
                        f_max = max(f_max,interp);
                    end
                    f_ave = f_sum*temp_inv;
                    f_harm = f_harm + (mod(k,2)*2-1) * f_sum;
                    f_var = sum_sqr - 2*f_ave.*f_sum + temp_dur*f_ave.*f_ave;
                    f_var(f_var<0) = 0; % in case value is slightly below 0
                    f_std = sqrt(f_var*temp_inv);
                    
                    %f_sum = sum(bout(start:finish,:))*temp_dur/ceil(temporal_grid_size(j));
                    %sum_sqr = sum(bout(start:finish,:).^2)*temp_dur/ceil(temporal_grid_size(j));
                    if p.use_bout_sum_features
                        feat(ind+1:ind+n_base_feat) = f_sum;
                        ind = ind + n_base_feat;
                    end
                    if p.use_bout_ave_features
                        feat(ind+1:ind+n_base_feat) = f_ave;
                        ind = ind + n_base_feat;
                    end
                    if p.use_bout_sum_absolute_features
                        feat(ind+1:ind+n_base_feat) = abs(f_sum);
                        ind = ind + n_base_feat;
                    end
                    if p.use_bout_ave_absolute_features
                        feat(ind+1:ind+n_base_feat) = abs(f_ave);
                        ind = ind + n_base_feat;
                    end
                    % compute std
                    if p.use_bout_sum_variance
                        feat(ind+1:ind+n_base_feat) = f_var;
                        ind = ind + n_base_feat;
                    end
                    if p.use_bout_standard_deviation
                        feat(ind+1:ind+n_base_feat) = f_std;
                        ind = ind + n_base_feat;
                    end
                    % compute min and max
                    if p.use_bout_max_feature
                        feat(ind+1:ind+n_base_feat) = f_max;
                        ind = ind + n_base_feat;
                    end
                    if p.use_bout_min_feature
                        feat(ind+1:ind+n_base_feat) = f_min;
                        ind = ind + n_base_feat;
                    end
                    start = finish+1;
                end

                % harmonic features    
                %     // An extension of the 1D Haar-like features, these capture changes or harmonic motion within the
                %     // bout.  When j=0, the response is the 1st half of the bout minus the end half of the bout.  When
                %     // j=1, the response is the 1st 3rd minus the 2nd 3rd plus the 3rd 3rd. And so on.
                %     // The total num of features (if absolute features are also used) is
                %     //   F_NUM_GOOD_GLOBAL_FEATURES*num_harmonic_features*2
                if p.use_sum_harmonic_features && k > 1
                    feat(ind+1:ind+n_base_feat) = f_harm;              
                    ind = ind + n_base_feat;
                end
                if p.use_ave_harmonic_features && k > 1
                    feat(ind+1:ind+n_base_feat) = f_harm*inv;              
                    ind = ind + n_base_feat;
                end
                if p.use_sum_absolute_harmonic_features && k > 1
                    feat(ind+1:ind+n_base_feat) = abs(f_harm);              
                    ind = ind + n_base_feat;
                end
                if p.use_ave_absolute_harmonic_features && k > 1
                    feat(ind+1:ind+n_base_feat) = abs(f_harm*inv);              
                    ind = ind + n_base_feat;
                end            
            end

            %% compute histogram
            %     // Discretize the feature space into different bins, and then compute the count of the number
            %     // of frames in each bin.  This is a histogram of feature values.  The histograms can also be
            %     // computed over different temporal intervals, yielding an output feature space of size
            %     //    F_NUM_GOOD_GLOBAL_FEATURES*num_histogram_bins*(sum_{l=0}^num_histogram_bins_temporal_levels(2^i))        
            for l=1:p.num_histogram_bins
                f_hist = integral_histogram_features(t_end+1,l,:) - ...
                         integral_histogram_features(t_start,l,:);
                if p.use_histogram_sum_features
                    feat(ind+1:ind+n_base_feat) = f_hist;
                    ind = ind + n_base_feat;
                end
                if p.use_histogram_ave_features
                    feat(ind+1:ind+n_base_feat) = f_hist*inv;                
                    ind = ind + n_base_feat;
                end
            end        

            %% compute global difference features
            %     // These are differences in the bout average or sum response for the bout as compared to the global max, min,
            %     // and average over the entire behavioral sequence
            %     // Should be compared to the min/max/ave of this behavior sequence, not the training one.
            %     // Mostly useful if features are not normalized beforehand
            if p.use_global_difference_max_sum_features
                feat(ind+1:ind+n_base_feat) = f_sum - max_feature_responses*dur;
                ind = ind + n_base_feat;
            end
            if p.use_global_difference_max_ave_features
                feat(ind+1:ind+n_base_feat) = f_ave - max_feature_responses;
                ind = ind + n_base_feat;
            end
            if p.use_global_difference_min_sum_features 
                feat(ind+1:ind+n_base_feat) = f_sum - min_feature_responses*dur;
                ind = ind + n_base_feat;
            end
            if p.use_global_difference_min_ave_features
                feat(ind+1:ind+n_base_feat) = f_ave - min_feature_responses;
                ind = ind + n_base_feat;
            end
            if p.use_global_difference_ave_sum_features
                feat(ind+1:ind+n_base_feat) = f_sum - ave_feature_responses*dur;
                ind = ind + n_base_feat;
            end
            if p.use_global_difference_ave_ave_features
                feat(ind+1:ind+n_base_feat) = f_ave - ave_feature_responses;
                ind = ind + n_base_feat;
            end

            %% compute bout change
            %     // The total change from the beginning of the bout to the end of the bout
            if p.use_bout_change || p.use_bout_absolute_change
                f_change = X(t_end,:)-X(t_start,:);
                if p.use_bout_change
                    feat(ind+1:ind+n_base_feat) = f_change;
                    ind = ind + n_base_feat;
                end
                if p.use_bout_absolute_change
                    feat(ind+1:ind+n_base_feat) = abs(f_change);
                    ind = ind + n_base_feat;
                end
            end

            %% compute haar features    
            %     // Kind of like 1D Haar-like features, these are the difference in sum feature response
            %     // in some temporal interval minus the difference in sum feature response of some other interval.
            %     // Can also use features of the absolute value of these feature responses
            %     // The total num of features (if absolute features are also used) is
            %     //   F_NUM_GOOD_GLOBAL_FEATURES*num_difference_temporal_levels*2*2

            %       // Difference in sum feature response in the region inside the bout (t_start,t_end) and 
            %       // the region of duration (t_end-t_start)/(2^j) immediately before t_start
            inv = 1/(td+1);
            if p.use_start_sum_diff_haar_features || p.use_start_ave_diff_haar_features || ...
                    p.use_start_sum_absolute_diff_haar_features || p.use_start_ave_absolute_diff_haar_features

                f_diff = integral_features(min(t_start+td+1,t_end),:) - integral_features(t_start,:) - ...
                         integral_features(t_start+1,:) + integral_features(max(t_start-td,1),:);
                if p.use_start_sum_diff_haar_features
                    feat(ind+1:ind+n_base_feat) = f_diff;               
                    ind = ind + n_base_feat;
                end
                if p.use_start_ave_diff_haar_features
                    feat(ind+1:ind+n_base_feat) = f_diff*inv;               
                    ind = ind + n_base_feat;
                end
                if p.use_start_sum_absolute_diff_haar_features
                    feat(ind+1:ind+n_base_feat) = abs(f_diff);               
                    ind = ind + n_base_feat;
                end
                if p.use_start_ave_absolute_diff_haar_features
                    feat(ind+1:ind+n_base_feat) = abs(f_diff*inv);               
                    ind = ind + n_base_feat;
                end
            end
            if p.use_end_sum_diff_haar_features || p.use_end_ave_diff_haar_features || ...
                    p.use_end_sum_absolute_diff_haar_features || p.use_end_ave_absolute_diff_haar_features
                %       // Difference in average feature response in the region inside the bout (t_start,t_end) and 
                %       // the region of duration (t_end-t_start)/(2^j) immediately after t_end
                f_diff = -integral_features(min(t_end+td+1,T+1),:) + integral_features(t_end,:) + ...
                          integral_features(t_end+1,:) - integral_features(max(t_end-td,t_start),:);
                if p.use_end_sum_diff_haar_features
                    feat(ind+1:ind+n_base_feat) = f_diff;               
                    ind = ind + n_base_feat;
                end
                if p.use_end_ave_diff_haar_features
                    feat(ind+1:ind+n_base_feat) = f_diff*inv;               
                    ind = ind + n_base_feat;
                end
                if p.use_end_sum_absolute_diff_haar_features
                    feat(ind+1:ind+n_base_feat) = abs(f_diff);               
                    ind = ind + n_base_feat;
                end
                if p.use_end_ave_absolute_diff_haar_features
                    feat(ind+1:ind+n_base_feat) = abs(f_diff*inv);               
                    ind = ind + n_base_feat;
                end
            end

            if p.use_duration_feature
                ind = ind+1;
                feat(ind) = dur/FPS;
            end

            if ind ~= n_total_feat
                disp('Warning: n_bout_feat ~= ind')
            end          
            boutfeat(b,:) = feat;
        end
        feat_mus{jj} = mean(boutfeat);
        toc
        % save results for current sequences
        save(savenames{jj},'X','bouts','behs','boutfeat','first_valid_fr','-v7.3');
        clear boutfeat
    end
    
    % normalize features if requested
    disp('computing normalizing factors')
    tic
    if nargin < 9 || isempty(trainvals)
        % compute normalizing parameters
        feat_mu = zeros(1,n_total_feat);
        total_bouts = 0;
        for jj=1:numel(Xs)
            feat_mu = feat_mu + feat_mus{jj}*size(Bouts{jj},1);
            total_bouts = total_bouts + size(Bouts{jj},1);
        end
        feat_mu = feat_mu/total_bouts;
        vars = zeros(numel(Xs),n_total_feat);
        for jj=1:numel(Xs)
            % load features
            D = load(savenames{jj});
            for b=1:size(Bouts{jj},1)
                vars(jj,:) = vars(jj,:) + (D.boutfeat(b,:) - feat_mu).^2;
            end
            clear D
        end
        feat_sigma = sqrt(sum(vars)/total_bouts);
        feat_gamma = 1./feat_sigma;
    else
        feat_mu = trainvals.feat_mu;
        feat_gamma = trainvals.feat_gamma;
    end
    toc
    if nargin >= 7 && normalize
        disp('normalizing bouts')
        tic
        for jj=1:numel(Xs)
            D = load(savenames{jj});
            X = D.X; bouts = D.bouts; behs = D.behs; boutfeat = D.boutfeat; first_valid_fr = D.first_valid_fr;
            clear D
            for b=1:size(Bouts{jj},1)
                boutfeat(b,:) = (boutfeat(b,:)-feat_mu) .* feat_gamma;  
            end
            save(savenames{jj},'X','bouts','behs','boutfeat','first_valid_fr','-v7.3');                             
            clear boutfeat
        end
        toc
    end
    
    if nargin < 9 || isempty(trainvals) && nargout > 1
        trainvals.histogram_thresholds  = histogram_thresholds;
        trainvals.ave_feature_responses = ave_feature_responses;
        trainvals.min_feature_responses = min_feature_responses;
        trainvals.max_feature_responses = max_feature_responses;    
        trainvals.feat_mu               = feat_mu;
        trainvals.feat_gamma            = feat_gamma;  
        trainvals.p = p;
    end
    disp('finished!')
end