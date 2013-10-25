
function pred_to_actions(pred_path, traintest, behs, extn)
    % set default variables
    if nargin < 1 || isempty(pred_path)
        pred_path = '/Users/eyrun/caltech/Code/Jdetect/fly_code_online2/data_ssvm';
    end
    if nargin < 2 || isempty(traintest)
        traintest = 'train';
    end
    if nargin < 3 || isempty(behs)
        behs = {'touch','lunge','wing_threat','charge','wing_extension'};
    end
    if nargin < 4 || isempty(extn)
        extn = '.16';
    elseif extn(1)~='.'
        extn = ['.' extn];
    end
    movies = {};
    all_bouts = {};
    
    % loop through all behaviors
    for b=1:numel(behs)
        % load prediction file
        pred_file = fullfile(pred_path,[traintest '_' behs{b} '.txt.pred' extn]);        
        if ~exist(pred_file,'file')
            disp(['File does not exist: ' pred_file])
            continue
        end
        data = loadjson(pred_file);
        % assign predictions to movie and behavior
        for d=1:numel(data)
            pred = data{d}.predicted.y;
            fly_id = data{d}.ground_truth.y.fly_id;
            C = strsplit(pred.moviename,'/');
            moviename = fullfile(C{end-2},C{end-1});
            if ismember(moviename,movies)
                movie_id = find(strcmp(moviename,movies));
            else
                movies{end+1} = moviename;
                all_bouts{end+1} = cell(2,numel(behs));
                movie_id = numel(movies);
            end
            tmp_bouts = zeros(floor(numel(pred.bouts)/2),4);
            count = 0;
            for bb=1:numel(pred.bouts)
                bout = pred.bouts(bb);
                if bout.behavior == 0
                    continue
                end
                count = count + 1;
                tmp_bouts(count,1) = bout.start_frame;
                tmp_bouts(count,2) = bout.end_frame;
                tmp_bouts(count,4) = bout.bout_score;
            end
            all_bouts{movie_id}{fly_id,b} = [all_bouts{movie_id}{fly_id,b}; tmp_bouts];
        end
    end
    
    parent_folders = cell(1,numel(movies));
    for m=1:numel(movies)
        parent_folders{m} = fileparts(movies{m});
    end
    parent_folders = unique(parent_folders);
    if numel(parent_folders) > 1
        parent_folder = C{end-3};
    else
        parent_folder = '';
    end
    
    % save consolidated behavior bouts into action files 
    % (readable by visualizer)
    savedir = fullfile(pred_path,'actions');
    for i=1:numel(parent_folders)
        curr_savedir = fullfile(savedir,parent_folders{i});
        if ~exist(curr_savedir,'dir')
            mkdir(curr_savedir)
        end    
    end
    for m=1:numel(movies)
        savename = fullfile(savedir,parent_folder,[movies{m} '_ssvm_actions.mat']);
        bouts = all_bouts{m};
        save(savename,'behs','bouts')
    end
end