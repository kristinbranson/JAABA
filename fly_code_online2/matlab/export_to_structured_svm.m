function params = export_to_structured_svm(p, trainvals, FPS, datadir, export_dir, do_feature_selection, train, test, detect_percent_overlap_train, max_frames_train, max_frames_test)

% p defines the set of bout features to use
if nargin < 1 || isempty(p), p = default_params(); end

% If non-empty, trainval encodes pre-learned histogram thresholds and bout
% feature normalization factors (output of generate_bouts_vec.m)
if nargin < 2, trainvals = {}; end

% The framerate of the video
if nargin < 3 || isempty(FPS), FPS = 200; end

% Folder with the movies to import (Eyrun's format)
if nargin < 4 || isempty(datadir), datadir = '/scratch/Datasets/FlyPairs/Jon_scored_movies/'; end
%if nargin < 4 || isempty(datadir), datadir = '/Users/eyrun/Caltech/Data/Jon_scored_movies/'; end

% Output folder (structured svm format)
if nargin < 5 || isempty(export_dir), export_dir = 'data_ssvm'; end

if nargin < 6 || isempty(do_feature_selection), do_feature_selection = 'false'; end

% Training movies
if nargin < 7 || isempty(train), train = { 'movie1', 'movie2', 'movie3' }; end

% Test movies
if nargin < 8 || isempty(test), test = { 'movie4', 'movie5', 'movie6' }; end

% If non-zero, purposefully leave unlabeled gaps.  The structured svm code
% will ignore these gaps, effectively not penalizing predictions that move
% the start/end of the bout by no more than detect_percent_overlap percent.
% This might give better results when the start/end of bouts aren't precisely 
% labeled or when
if nargin < 9 || isempty(detect_percent_overlap_train), detect_percent_overlap_train = 0; end
detect_percent_overlap_test = 0;

% Optionally break apart training movies into multiple segments for better
% computational tractability
if nargin < 10 || isempty(max_frames_train), max_frames_train = 100000; end

if nargin < 11 || isempty(max_frames_test), max_frames_test = 100000; end
    

[~,behs,feat]=loadmovie(fullfile(datadir, train{1}));

% Save behavior and feature parameters
disp('computing feature definitions...')
params = {};
params.behaviors = behavior_parameters(behs);
params.frame_feature_params = frame_feature_parameters(feat, p, trainvals);
%if do_feature_selection,
%    params.maxDual = .03;
%    params.bout_expansion_feature_params = bout_feature_parameters(params.frame_feature_params, p, trainvals, FPS);
%else
    params.bout_feature_params = bout_feature_parameters(params.frame_feature_params, p, trainvals, FPS);
%end
if ~isempty(trainvals) && ~isempty(trainvals.feat_mu) && ~isempty(trainvals.feat_gamma) && ~isempty(trainvals.histogram_thresholds) 
    params.dontComputeFeaturMeanVarianceMedianStatistics = true;
end

mkdir(export_dir);


% export train/test data
disp('exporting frame-level features...')
movies = [train(:);test(:)]';

    
% Save features and bout definitions
basenames = cell(numel(movies),1);
for ii=1:(numel(behs)+1),
    trainsets{ii} = {};
    testsets{ii} = {};
end
for ii=1:numel(movies),
    [bouts,behs,feat,basenames{ii},moviename] = loadmovie(fullfile(datadir, movies{ii}));
    for fly_id=1:size(feat.data,1),
        if ii <= numel(train),
            [valid_beh, trainsets]=save_example(export_dir, basenames{ii}, moviename, bouts, behs, feat, fly_id, FPS, params.frame_feature_params, trainsets, detect_percent_overlap_train, max_frames_train);
        else
            try
            [valid_beh, testsets]=save_example(export_dir, basenames{ii}, moviename, bouts, behs, feat, fly_id, FPS, params.frame_feature_params, testsets, detect_percent_overlap_test, max_frames_test);
            catch 
                a=1;
            end
        end
    end
end


% save different train/test sets for each binary behavior as well as a
% multiclass version
disp('saving behavior and feature definitions...')
for ii=1:(numel(behs)+1),
    params2 = params;
    if ii==(numel(behs)+1), 
        beh = 'multiclass'; 
        params2.behaviors = [params.behaviors(1),params.behaviors(valid_beh+1)'];
    else 
        beh = behs{ii}; 
        params2.behaviors = [params.behaviors(1),params.behaviors(ii+1)];
    end
    
    % Save behavior and feature definitions
    fout = fopen(sprintf('%s/params_%s.txt', export_dir, beh), 'w');
    fprintf(fout, '%s', savejson(params2));
    fclose(fout);
    
    fout = fopen(sprintf('%s/params_%s_fs.txt', export_dir, beh), 'w');
    params3 = rmfield(params2,'bout_feature_params');
    params3.maxDual = .03;
    params3.bout_expansion_feature_params = params2.bout_feature_params;
    fprintf(fout, '%s', savejson(params3));
    fclose(fout);
    
    % Save train/test file lists
    fout=fopen(fullfile(export_dir, sprintf('train_%s.txt',beh)),'w');
    for jj=1:length(trainsets{ii}),
        fprintf(fout, '%s\n', strrep(strrep(savejson(trainsets{ii}(jj)),sprintf('\n'),''),sprintf('\t'),''));
    end
    fclose(fout);
    
    fout=fopen(fullfile(export_dir, sprintf('test_%s.txt',beh)),'w');
    for jj=1:length(testsets{ii}),
        fprintf(fout, '%s\n', strrep(strrep(savejson(testsets{ii}(jj)),sprintf('\n'),''),sprintf('\t'),''));
    end
    fclose(fout);
end



function [bouts,behs,feat,basename,moviename] = loadmovie(moviedir)
    f = dir(fullfile(moviedir, '*updated_actions.mat'));
    b = load(fullfile(moviedir,f.name));
    bouts = b.bouts;
    behs = b.behs;
    
    
    if isempty(strfind(f.name,'-updated_actions.mat')),
        basename = f.name(1:(length(f.name)-length('_actions.mat')));
    else
        basename = f.name(1:(length(f.name)-length('-updated_actions.mat')));
    end
    b = load(fullfile(moviedir, basename, strcat(basename,'-feat.mat')));
    feat = b.feat;

    moviename = fullfile(moviedir,strcat(basename,'.seq'));
end


% Initialize behavior definitions
function behaviors = behavior_parameters(behs)
    behaviors = cell(numel(behs)+1,1);
    behaviors{1}.name = 'none';
    for ii=1:numel(behs),
        behaviors{ii+1}.name = behs{ii};
    end
end

% Initialize frame feature definitions
function base_feat = frame_feature_parameters(feat, p, trainvals)
    base_feat = cell(numel(feat.names),1); 
    for ii=1:numel(feat.names)
        base_feat{ii}.name = feat.names{ii}; 
        base_feat{ii}.num_histogram_bins = p.num_histogram_bins;
        if ~isempty(trainvals) && ~isempty(trainvals.histogram_thresholds)
            base_feat{ii}.histogram_thresholds = trainvals.histogram_thresholds(ii,:);
        end
    end
end


% Initialize bout-level features
function bout_feat = bout_feature_parameters(base_feat, p, trainvals, FPS)
    n_base_feat = numel(base_feat);
    bout_feat = {}; 
    ind = 1;
    td = max(round(1/30*FPS),1); % time delta used for haar like features
    
    temporal_levels = 1:p.num_temporal_levels;
    temporal_grid_size = zeros(size(temporal_levels));  
    for i=1:p.num_temporal_levels
        temporal_grid_size(i) = 1/temporal_levels(i);
    end          
    for j=p.num_temporal_levels:-1:1
        for k=1:temporal_levels(j)
            start = (k-1)*temporal_grid_size(j);
            finish = k*temporal_grid_size(j);
            
            if p.use_bout_sum_features
                for ii=1:n_base_feat,
                    r = struct('op', 'sum', 'frame_feature', ii-1, 't_start', start, 't_end', finish);
                    f = struct('regions', {r}, 'weights', {1}, 'name', ...
                            sprintf('%s-sum(%.3f,%.3f)', base_feat{ii}.name, start, finish));
                    bout_feat{ind} = f;
                    ind = ind+1;
                end
            end
            if p.use_bout_ave_features
                for ii=1:n_base_feat,
                    r = struct('op', 'ave', 'frame_feature', ii-1, 't_start', start, 't_end', finish);
                    f = struct('regions', {r}, 'weights', {1}, 'name', ...
                            sprintf('%s-ave(%.3f,%.3f)', base_feat{ii}.name, start, finish));
                    bout_feat{ind} = f;
                    ind = ind+1;
                end
            end
            if p.use_bout_sum_absolute_features
                for ii=1:n_base_feat,
                    r = struct('op', 'sum', 'frame_feature', ii-1, 't_start', start, 't_end', finish);
                    f = struct('regions', {r}, 'weights', {1}, 'absolute_value', true, 'name', ...
                            sprintf('%s-abs_sum(%.3f,%.3f)', base_feat{ii}.name, start, finish));
                    bout_feat{ind} = f;
                    ind = ind+1;
                end
            end
            if p.use_bout_ave_absolute_features
                for ii=1:n_base_feat,
                    r = struct('op', 'ave', 'frame_feature', ii-1, 't_start', start, 't_end', finish);
                    f = struct('regions', {r}, 'weights', {1}, 'absolute_value', true, 'name', ...
                            sprintf('%s-abs_ave(%.3f,%.3f)', base_feat{ii}.name, start, finish));
                    bout_feat{ind} = f;
                    ind = ind+1;
                end
            end
            if p.use_bout_sum_variance
                for ii=1:n_base_feat,
                    r = struct('op', 'var', 'frame_feature', ii-1, 't_start', start, 't_end', finish);
                    f = struct('regions', {r}, 'weights', {1}, 'name', ...
                            sprintf('%s-sum_var(%.3f,%.3f)', base_feat{ii}.name, start, finish));
                    bout_feat{ind} = f;
                    ind = ind+1;
                end
            end
            if p.use_bout_standard_deviation
                for ii=1:n_base_feat,
                    r = struct('op', 'dev', 'frame_feature', ii-1, 't_start', start, 't_end', finish);
                    f = struct('regions', {r}, 'weights', {1}, 'name', ...
                            sprintf('%s-dev(%.3f,%.3f)', base_feat{ii}.name, start, finish));
                    bout_feat{ind} = f;
                    ind = ind+1;
                end
            end
            if p.use_bout_max_feature
                for ii=1:n_base_feat,
                    r = struct('op', 'max', 'frame_feature', ii-1, 't_start', start, 't_end', finish);
                    f = struct('regions', {r}, 'weights', {1}, 'name', ...
                            sprintf('%s-max(%.3f,%.3f)', base_feat{ii}.name, start, finish));
                    bout_feat{ind} = f;
                    ind = ind+1;
                end
            end
            if p.use_bout_min_feature
                for ii=1:n_base_feat,
                    r = struct('op', 'min', 'frame_feature', ii-1, 't_start', start, 't_end', finish);
                    f = struct('regions', {r}, 'weights', {1}, 'name', ...
                            sprintf('%s-min(%.3f,%.3f)', base_feat{ii}.name, start, finish));
                    bout_feat{ind} = f;
                    ind = ind+1;
                end
            end
        end
        if p.use_sum_harmonic_features && temporal_levels(j) > 1,
            for ii=1:n_base_feat,
                f = {};
                f.weights = zeros(1,temporal_levels(j));
                f.regions = cell(1,temporal_levels(j));
                f.name = sprintf('%s-sum_harmonic_%d', base_feat{ii}.name, temporal_levels(j));
                for k=1:temporal_levels(j),
                    start = (k-1)*temporal_grid_size(j);
                    finish = (k)*temporal_grid_size(j);
                    f.regions{k} = struct('op', 'sum', 'frame_feature', ii-1, 't_start', start, 't_end', finish);
                    f.weights(k) = mod(k,2)*2-1;
                end
                bout_feat{ind} = f;
                ind = ind+1;
            end
        end
        if p.use_ave_harmonic_features && temporal_levels(j) > 1,
            for ii=1:n_base_feat,
                f = {};
                f.weights = zeros(1,temporal_levels(j));
                f.regions = cell(1,temporal_levels(j));
                f.name = sprintf('%s-ave_harmonic_%d', base_feat{ii}.name, temporal_levels(j));
                f.time_normalize = true;  % normalization occurs on the total response
                for k=1:temporal_levels(j),
                    start = (k-1)*temporal_grid_size(j);
                    finish = (k)*temporal_grid_size(j);
                    f.regions{k} = struct('op', 'sum', 'frame_feature', ii-1, 't_start', start, 't_end', finish);
                    f.weights(k) = mod(k,2)*2-1;
                end
                bout_feat{ind} = f;
                ind = ind+1;
            end
        end
        if p.use_sum_absolute_harmonic_features && temporal_levels(j) > 1,
            for ii=1:n_base_feat,
                f = {};
                f.weights = zeros(1,temporal_levels(j));
                f.regions = cell(1,temporal_levels(j));
                f.name = sprintf('%s-abs_sum_harmonic_%d', base_feat{ii}.name, temporal_levels(j));
                f.absolute_value = true;
                for k=1:temporal_levels(j),
                    start = (k-1)*temporal_grid_size(j);
                    finish = (k)*temporal_grid_size(j);
                    f.regions{k} = struct('op', 'sum', 'frame_feature', ii-1, 't_start', start, 't_end', finish);
                    f.weights(k) = mod(k,2)*2-1;
                end
                bout_feat{ind} = f;
                ind = ind+1;
            end
        end
        if p.use_ave_absolute_harmonic_features && temporal_levels(j) > 1,
            for ii=1:n_base_feat,
                f = {};
                f.weights = zeros(1,temporal_levels(j));
                f.regions = cell(1,temporal_levels(j));
                f.name = sprintf('%s-abs_ave_harmonic_%d', base_feat{ii}.name, temporal_levels(j));
                f.time_normalize = true;  % normalization occurs on the total response
                f.absolute_value = true;
                for k=1:temporal_levels(j),
                    start = (k-1)*temporal_grid_size(j);
                    finish = (k)*temporal_grid_size(j);
                    f.regions{k} = struct('op', 'sum', 'frame_feature', ii-1, 't_start', start, 't_end', finish);
                    f.weights(k) = mod(k,2)*2-1;
                end
                bout_feat{ind} = f;
                ind = ind+1;
            end
        end
    end
    
    for l=1:p.num_histogram_bins
        if p.use_histogram_sum_features
            for ii=1:n_base_feat,
                r = struct('op', 'sum', 'hist_bin', l-1, 'frame_feature', ii-1, 't_start', 0, 't_end', 1);
                f = struct('regions', {r}, 'weights', {1}, 'name', ...
                            sprintf('%s-hist%d_sum(%.3f,%.3f)', base_feat{ii}.name, l, start, finish));
                bout_feat{ind} = f;
                ind = ind+1;
            end
        end
        if p.use_histogram_ave_features
            for ii=1:n_base_feat,
                r = struct('op', 'ave', 'hist_bin', l-1, 'frame_feature', ii-1, 't_start', 0, 't_end', 1);
                f = struct('regions', {r}, 'weights', {1}, 'name', ...
                            sprintf('%s-hist%d_ave(%.3f,%.3f)', base_feat{ii}.name, l, start, finish));
                bout_feat{ind} = f;
                ind = ind+1;
            end
        end
    end
    if p.use_global_difference_max_sum_features
        for ii=1:n_base_feat,
            f = {};
            f.weights = [1,-1];
            f.regions = cell(1,2);
            f.name = sprintf('%s-global_max_sum_diff', base_feat{ii}.name);
            f.regions{1} = struct('op', 'sum', 'frame_feature', ii-1, 't_start', 0, 't_end', 1);
            f.regions{2} = struct('op', 'max', 'frame_feature', ii-1,  'is_global', true, 'mult_dur', true);
            bout_feat{ind} = f;
            ind = ind+1;
        end
    end
    if p.use_global_difference_max_ave_features
        for ii=1:n_base_feat,
            f = {};
            f.weights = [1,-1];
            f.regions = cell(1,2);
            f.name = sprintf('%s-global_max_ave_diff', base_feat{ii}.name);
            f.regions{1} = struct('op', 'ave', 'frame_feature', ii-1, 't_start', 0, 't_end', 1);
            f.regions{2} = struct('op', 'max', 'frame_feature', ii-1,  'is_global', true);
            bout_feat{ind} = f;
            ind = ind+1;
        end
    end
    if p.use_global_difference_min_sum_features
        for ii=1:n_base_feat,
            f = {};
            f.weights = [1,-1];
            f.regions = cell(1,2);
            f.name = sprintf('%s-global_min_sum_diff', base_feat{ii}.name);
            f.regions{1} = struct('op', 'sum', 'frame_feature', ii-1, 't_start', 0, 't_end', 1);
            f.regions{2} = struct('op', 'min', 'frame_feature', ii-1,  'is_global', true, 'mult_dur', true);
            bout_feat{ind} = f;
            ind = ind+1;
        end
    end
    if p.use_global_difference_min_ave_features
        for ii=1:n_base_feat,
            f = {};
            f.weights = [1,-1];
            f.regions = cell(1,2);
            f.name = sprintf('%s-global_min_ave_diff', base_feat{ii}.name);
            f.regions{1} = struct('op', 'ave', 'frame_feature', ii-1, 't_start', 0, 't_end', 1);
            f.regions{2} = struct('op', 'min', 'frame_feature', ii-1,  'is_global', true);
            bout_feat{ind} = f;
            ind = ind+1;
        end
    end
    if p.use_global_difference_ave_sum_features
        for ii=1:n_base_feat,
            f = {};
            f.weights = [1,-1];
            f.regions = cell(1,2);
            f.name = sprintf('%s-global_ave_sum_diff', base_feat{ii}.name);
            f.regions{1} = struct('op', 'sum', 'frame_feature', ii-1, 't_start', 0, 't_end', 1);
            f.regions{2} = struct('op', 'ave', 'frame_feature', ii-1,  'is_global', true, 'mult_dur', true);
            bout_feat{ind} = f;
            ind = ind+1;
        end
    end
    if p.use_global_difference_ave_ave_features
        for ii=1:n_base_feat,
            f = {};
            f.weights = [1,-1];
            f.regions = cell(1,2);
            f.name = sprintf('%s-global_ave_ave_diff', base_feat{ii}.name);
            f.regions{1} = struct('op', 'ave', 'frame_feature', ii-1, 't_start', 0, 't_end', 1);
            f.regions{2} = struct('op', 'ave', 'frame_feature', ii-1,  'is_global', true);
            bout_feat{ind} = f;
            ind = ind+1;
        end
    end
    if p.use_bout_change
        for ii=1:n_base_feat,
            r = struct('op', 'raw', 'frame_feature', ii-1, 't_start', 0, 't_end', 1);
            f = struct('regions', {r}, 'weights', {1}, 'name', ...
                    sprintf('%s-bout_change(%.3f,%.3f)', base_feat{ii}.name, start, finish));
            bout_feat{ind} = f;
            ind = ind+1;
        end
    end
    if p.use_bout_absolute_change
        for ii=1:n_base_feat,
            r = struct('op', 'raw', 'frame_feature', ii-1, 't_start', 0, 't_end', 1);
            f = struct('regions', {r}, 'weights', {1}, 'absolute_value', true, 'name', ...
                    sprintf('%s-bout_change(%.3f,%.3f)', base_feat{ii}.name, start, finish));
            bout_feat{ind} = f;
            ind = ind+1;
        end
    end
    
    if p.use_start_sum_diff_haar_features
        for ii=1:n_base_feat,
            f = {};
            f.weights = [1,-1];
            f.regions = cell(1,2);
            f.name = sprintf('%s-start_sum_diff_haar', base_feat{ii}.name);
            %f.regions{1} = struct('op', 'sum', 'frame_feature', ii-1, 'frame_coords', 1, 't_start', 0, 't_end', td);
            %f.regions{2} = struct('op', 'sum', 'frame_feature', ii-1, 'frame_coords', 1, 't_start', -td, 't_end', 0);
            f.regions{1} = struct('op', 'sum', 'frame_feature', ii-1, 'frame_coords', 1, 't_start', 0, 't_end', td+1);
            f.regions{2} = struct('op', 'sum', 'frame_feature', ii-1, 'frame_coords', 1, 't_start', -td, 't_end', 1);
            bout_feat{ind} = f;
            ind = ind+1;
        end
    end
    if p.use_start_ave_diff_haar_features
        for ii=1:n_base_feat,
            f = {};
            f.weights = [1,-1];
            f.regions = cell(1,2);
            f.name = sprintf('%s-start_ave_diff_haar', base_feat{ii}.name);
            %f.regions{1} = struct('op', 'ave', 'frame_feature', ii-1, 'frame_coords', 1, 't_start', 0, 't_end', td);
            %f.regions{2} = struct('op', 'ave', 'frame_feature', ii-1, 'frame_coords', 1, 't_start', -td, 't_end', 0);
            f.regions{1} = struct('op', 'ave', 'frame_feature', ii-1, 'frame_coords', 1, 't_start', 0, 't_end', td+1);
            f.regions{2} = struct('op', 'ave', 'frame_feature', ii-1, 'frame_coords', 1, 't_start', -td, 't_end', 1);
            bout_feat{ind} = f;
            ind = ind+1;
        end
    end
    if p.use_start_sum_absolute_diff_haar_features
        for ii=1:n_base_feat,
            f = {};
            f.weights = [1,-1];
            f.regions = cell(1,2);
            f.name = sprintf('%s-start_abs_sum_diff_haar', base_feat{ii}.name);
            f.absolute_value = true;
            %f.regions{1} = struct('op', 'sum', 'frame_feature', ii-1, 'frame_coords', 1, 't_start', 0, 't_end', td);
            %f.regions{2} = struct('op', 'sum', 'frame_feature', ii-1, 'frame_coords', 1, 't_start', -td, 't_end', 0);
            f.regions{1} = struct('op', 'sum', 'frame_feature', ii-1, 'frame_coords', 1, 't_start', 0, 't_end', td+1);
            f.regions{2} = struct('op', 'sum', 'frame_feature', ii-1, 'frame_coords', 1, 't_start', -td, 't_end', 1);
            bout_feat{ind} = f;
            ind = ind+1;
        end
    end
    if p.use_start_ave_absolute_diff_haar_features
        for ii=1:n_base_feat,
            f = {};
            f.weights = [1,-1];
            f.regions = cell(1,2);
            f.name = sprintf('%s-start_abs_ave_diff_haar', base_feat{ii}.name);
            f.absolute_value = true;
            %f.regions{1} = struct('op', 'ave', 'frame_feature', ii-1, 'frame_coords', 1, 't_start', 0, 't_end', td);
            %f.regions{2} = struct('op', 'ave', 'frame_feature', ii-1, 'frame_coords', 1, 't_start', -td, 't_end', 0);
            f.regions{1} = struct('op', 'ave', 'frame_feature', ii-1, 'frame_coords', 1, 't_start', 0, 't_end', td+1);
            f.regions{2} = struct('op', 'ave', 'frame_feature', ii-1, 'frame_coords', 1, 't_start', -td, 't_end', 1);
            bout_feat{ind} = f;
            ind = ind+1;
        end
    end
    if p.use_end_sum_diff_haar_features
        for ii=1:n_base_feat,
            f = {};
            f.weights = [-1,1];
            f.regions = cell(1,2);
            f.name = sprintf('%s-end_sum_diff_haar', base_feat{ii}.name);
            %f.regions{1} = struct('op', 'sum', 'frame_feature', ii-1, 'frame_coords', 2, 't_start', 0, 't_end', td);
            %f.regions{2} = struct('op', 'sum', 'frame_feature', ii-1, 'frame_coords', 2, 't_start', -td, 't_end', 0);
            f.regions{1} = struct('op', 'sum', 'frame_feature', ii-1, 'frame_coords', 2, 't_start', 0, 't_end', td+1);
            f.regions{2} = struct('op', 'sum', 'frame_feature', ii-1, 'frame_coords', 2, 't_start', -td, 't_end', 1);
            bout_feat{ind} = f;
            ind = ind+1;
        end
    end
    if p.use_end_ave_diff_haar_features
        for ii=1:n_base_feat,
            f = {};
            f.weights = [-1,1];
            f.regions = cell(1,2);
            f.name = sprintf('%s-end_ave_diff_haar', base_feat{ii}.name);
            %f.regions{1} = struct('op', 'ave', 'frame_feature', ii-1, 'frame_coords', 2, 't_start', 0, 't_end', td);
            %f.regions{2} = struct('op', 'ave', 'frame_feature', ii-1, 'frame_coords', 2, 't_start', -td, 't_end', 0);
            f.regions{1} = struct('op', 'ave', 'frame_feature', ii-1, 'frame_coords', 2, 't_start', 0, 't_end', td+1);
            f.regions{2} = struct('op', 'ave', 'frame_feature', ii-1, 'frame_coords', 2, 't_start', -td, 't_end', 1);
            bout_feat{ind} = f;
            ind = ind+1;
        end
    end
    if p.use_end_sum_absolute_diff_haar_features
        for ii=1:n_base_feat,
            f = {};
            f.weights = [-1,1];
            f.regions = cell(1,2);
            f.name = sprintf('%s-end_abs_sum_diff_haar', base_feat{ii}.name);
            f.absolute_value = true;
            %f.regions{1} = struct('op', 'sum', 'frame_feature', ii-1, 'frame_coords', 2, 't_start', 0, 't_end', td);
            %f.regions{2} = struct('op', 'sum', 'frame_feature', ii-1, 'frame_coords', 2, 't_start', -td, 't_end', 0);
            f.regions{1} = struct('op', 'sum', 'frame_feature', ii-1, 'frame_coords', 2, 't_start', 0, 't_end', td+1);
            f.regions{2} = struct('op', 'sum', 'frame_feature', ii-1, 'frame_coords', 2, 't_start', -td, 't_end', 1);
            bout_feat{ind} = f;
            ind = ind+1;
        end
    end
    if p.use_end_ave_absolute_diff_haar_features
        for ii=1:n_base_feat,
            f = {};
            f.weights = [-1,1];
            f.regions = cell(1,2);
            f.name = sprintf('%s-end_abs_ave_diff_haar', base_feat{ii}.name);
            f.absolute_value = true;
            %f.regions{1} = struct('op', 'ave', 'frame_feature', ii-1, 'frame_coords', 2, 't_start', 0, 't_end', td);
            %f.regions{2} = struct('op', 'ave', 'frame_feature', ii-1, 'frame_coords', 2, 't_start', -td, 't_end', 0);
            f.regions{1} = struct('op', 'ave', 'frame_feature', ii-1, 'frame_coords', 2, 't_start', 0, 't_end', td+1);
            f.regions{2} = struct('op', 'ave', 'frame_feature', ii-1, 'frame_coords', 2, 't_start', -td, 't_end', 1);
            bout_feat{ind} = f;
            ind = ind+1;
        end
    end 
    
    if p.use_duration_feature
        % note: B_DUR feature divides by the fps parameter in trx file
        r = struct('op', 'dur', 'frame_feature', -1, 't_start', 0, 't_end', 1);
        f = struct('regions', {r}, 'weights', {1}, 'name', ...
                sprintf('%s-duration', base_feat{ii}.name));
        bout_feat{ind} = f;
        ind = ind+1;
    end
    if p.num_duration_hist_bins
        % note: B_DUR feature divides by the fps parameter in trx file
        r = struct('op', 'dur', 'frame_feature', -1, 't_start', 0, 't_end', 1);
        f = struct('regions', {r}, 'weights', {1}, 'name', ...
                sprintf('%s-duration', base_feat{ii}.name), 'num_thresholds', p.num_duration_hist_bins);
        for jj=1:p.num_duration_hist_bins,
            bout_feat{ind} = f;
            ind = ind+1;
        end
    end  
    
    % set feature normalization factors
    if ~isempty(trainvals) && ~isempty(trainvals.feat_mu) && ~isempty(trainvals.feat_gamma)
        for jj=1:length(trainvals.feat_mu)
            bout_feat{jj}.mu = trainvals.feat_mu(jj);
            bout_feat{jj}.gamma = trainvals.feat_gamma(jj);
        end
    end
end



% Export frame-level features into Kristin's format, which is readable 
% by the SVM struct code
function [valid_beh, datasets] = save_example(export_dir, basename, moviename, bouts, behs, feat, fly_id, FPS, base_feat, datasets, detect_percent_overlap, max_frames)
    name = sprintf('%s_%d', basename, fly_id);
    dirname = fullfile(export_dir, name);
    mkdir(dirname);
        
    T = size(feat.data,2);
    f=reshape(feat.data(fly_id,:,:),size(feat.data,2),size(feat.data,3));
    trx = {};
    
    % set first frame to be the first frame where both flies present
    tmp = find(~isnan(feat.data(1,:,1)));
    first1 = tmp(1);
    tmp = find(~isnan(feat.data(2,:,1)));
    first2 = tmp(1);
    trx.firstframe = max(first1,first2);
    trx.endframe = T;
    %trx.firstframe = 1;
    %trx.endframe = feat.end_frame;
    
    % handle in-between and edge nan values
    for j=3:size(f,2)
        tmp = find(~isnan(f(:,j)));
        % interpolate in-between nan values
        cc = bwconncomp(~isnan(f(:,j)));
        for c=1:cc.NumObjects
            prev_fr = cc.PixelIdxList{c}(1);
            next_fr = cc.PixelIdxList{c}(end);
            if prev_fr < 1 || next_fr > trx.endframe
                continue
            end
            prev_val = f(prev_fr,j);
            next_val = f(next_fr,j);
            n_fr = next_fr-prev_fr+1;
            d_val = (next_val-prev_val)/(n_fr-1);
            f(prev_fr:next_fr,j) = prev_val + (0:n_fr-1)*d_val;
        end
        % extrapolate for edge nan values
        last_invalid_fr = tmp(1)-1;
        if last_invalid_fr+1 > trx.firstframe
            f(trx.firstframe:last_invalid_fr,j) = f(last_invalid_fr+1,j);
        end
        last_valid_fr = tmp(end);
        if last_valid_fr < size(f,1)
            f(last_valid_fr+1:end,j) = f(last_valid_fr,j);
        end
    end
    
    [m, ~]=ind2sub(size(f),find(isnan(f)));
    if ~isempty(m), 
        if feat.end_frame-min(m) < max(m)-1,
            trx.endframe = min(m)-1;
        else
            trx.firstframe = max(m)+1;
        end
        discarded = feat.end_frame-1-(trx.endframe-trx.firstframe);
        disp(sprintf('Discarding %d frames in %s due to nan valued features', discarded, dirname));
        %assert(discarded < 100);
    end
    trx.id = fly_id;
    trx.fps = FPS;
    trx.moviefile = '';
    trx.moviename = moviename;
    trx.sex = 'm';
    %trx.timestamps = (0:(T-1))/FPS;
            
    % Save features
    mkdir(fullfile(dirname, 'perframe'));
    for kk=1:size(f,2)
        data = {f(:,kk)};
        save(sprintf('%s/perframe/%s.mat',dirname,base_feat{kk}.name), 'data');
    end
    
    % Save bouts
    all_t0s=[];  all_t1s=[];  all_beh=[];   valid_beh = [];
    for kk=1:numel(behs)
        % Save for binary training
        if isempty(bouts{fly_id,kk}),
            t0s = [];   t1s = [];  be = [];
        else
            t0s = max(bouts{fly_id,kk}(:,1), trx.firstframe);
            t1s = min(bouts{fly_id,kk}(:,2), trx.endframe);
            be = ones(length(t0s),1)*(kk+1);
        end
        if ~strcmp(behs{kk}, 'unknown')
            all_t0s = [all_t0s; t0s];
            all_t1s = [all_t1s; t1s];
            all_beh = [all_beh; be];
            valid_beh = [valid_beh; kk];
        end
        d = save_labeled_bouts(t0s, t1s, be, export_dir, name, behs{kk}, behs, trx, detect_percent_overlap, max_frames);
        if isempty(datasets{kk}), 
            datasets{kk} = d; 
        else
            datasets{kk}((end+1):(end+length(d))) = d;
        end
    end
    
    % Save bouts for multiclass training
    d=save_labeled_bouts(all_t0s, all_t1s, all_beh, export_dir, name, 'multiclass', behs, trx, detect_percent_overlap, max_frames);
    if isempty(datasets{end}), 
            datasets{end} = d; 
        else
            datasets{end}((end+1):(end+length(d))) = d;
    end
end

% The structured SVM code assumes the set of labeled bouts spans all time
% frames and is sorted by time.  Any portion of the video that is unlabeled
% is ignored.  This function currently fills all unlabeled regions with the
% none behavior.
function dataset=save_labeled_bouts(t0s, t1s, beh, export_dir, name, beh_name, behs, trx, detect_percent_overlap, max_frames)
    behs2 = [ 'none', behs];
    [~, inds] = sort(t0s);
    t0s = t0s(inds);
    t1s = t1s(inds);
    beh = beh(inds);
    prev = trx.firstframe;
    dataset = {};
    n = 1;
    b = {};
    b(n).r_t0s = [];  b(n).r_t1s = [];  b(n).r_beh = [];   b(n).trx = trx;
    for i=1:(length(t0s)+1)
        if i <= length(t0s), 
            o = detect_percent_overlap*(t1s(i)-t0s(i)+1)/2;
            s = ceil(t0s(i)-o); 
        else
            o = 0;
            s = trx.endframe; 
        end
        if s > prev
            % If total movie length exceeds some maximum number of frames, break
            % the movie apart in the middle of frames with 'none' behavior
            while s-trx.firstframe > max_frames,
                if trx.firstframe+max_frames >= prev, 
                    b(n).r_t0s(end+1) = prev;
                    b(n).r_t1s(end+1) = trx.firstframe+max_frames;
                    b(n).r_beh(end+1) = 1;
                    b(n).trx.endframe = trx.firstframe+max_frames;
                    prev = trx.firstframe+max_frames;
                    trx.firstframe = trx.firstframe+max_frames;
                else
                    b(n).trx.endframe =  b(n).r_t1s(end);
                    trx.firstframe = prev;
                end
                
                n = n + 1;
                b(n).r_t0s = [];  b(n).r_t1s = [];  b(n).r_beh = [];   b(n).trx = trx;
            end
            
            % append 'none' behaviors inside unlabeled time intervals
            b(n).r_t0s(end+1) = prev;
            b(n).r_t1s(end+1) = s;
            b(n).r_beh(end+1) = 1;
        end
        if i <= length(t0s),
            b(n).r_t0s(end+1) = t0s(i);
            b(n).r_t1s(end+1) = t1s(i)+1;
            b(n).r_beh(end+1) = beh(i);
            prev = floor(t1s(i)+1+o);
        end
    end
        
    for i=1:numel(b),
        for tt=2:(length(b(i).r_t0s))
            if b(i).r_t0s(tt) < b(i).r_t1s(tt-1),
                display(sprintf('Bouts (%s,%d,%d) and (%s,%d,%d) overlap', behs2{b(i).r_beh(tt-1)}, b(i).r_t0s(tt-1), b(i).r_t1s(tt-1), behs2{b(i).r_beh(tt)}, b(i).r_t0s(tt), b(i).r_t1s(tt)));
                b(i).r_t0s(tt) = round((b(i).r_t0s(tt)+b(i).r_t1s(tt-1))/2);
                b(i).r_t1s(tt-1) = b(i).r_t0s(tt);
            end
        end
        t0s = { b(i).r_t0s };
        t1s = { b(i).r_t1s };
        names = behs2(b(i).r_beh);
        trx = b(i).trx;
        e = {};  e.x = {};  e.y = {};
        e.y.fname = fullfile(export_dir,name,sprintf('labeled_gt_%s_%d.mat',beh_name,i));
        e.x.fname = fullfile(export_dir,name,sprintf('trx_%s_%d.mat',beh_name,i));
        display(sprintf('%s_%d_%d: %s %s', beh_name, fly_id, i, e.x.fname, e.y.fname));
        save(e.x.fname, 'trx');
        save(e.y.fname, 't0s', 't1s', 'names');
        if i==1, dataset=e; else, dataset(i) = e; end
    end
end

function p=default_params()
        p = {};
        p.add_wavelets                                = false;
        p.num_histogram_bins                          = 8;
        p.num_temporal_levels                         = 3;
        p.use_histogram_ave_features                  = true;
        p.use_histogram_sum_features                  = false; %true; %
        p.use_bout_ave_absolute_features              = true;
        p.use_bout_ave_features                       = true;
        p.use_bout_max_feature                        = true;
        p.use_bout_min_feature                        = true;
        p.use_bout_sum_absolute_features              = false; %true; %
        p.use_bout_sum_features                       = false; %true; %
        p.use_bout_standard_deviation                 = true;
        p.use_bout_sum_variance                       = false; %true; %
        p.use_ave_absolute_harmonic_features          = true;
        p.use_ave_harmonic_features                   = true;
        p.use_sum_absolute_harmonic_features          = false; %true; %
        p.use_sum_harmonic_features                   = false; %true; %
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
        p.num_duration_hist_bins                     = 0; %20; %
end

end
