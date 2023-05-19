function copyhoghofParamstoJson(exp_file, params_file, json_dir_file,varargin)

    % file name defaults 
    [output_file_HOG, output_file_HOF, ...
     output_file_CropSide, output_file_CropFront] = myparse(varargin,...
    'output_file_HOG', 'HOGParam.json', ...
    'output_file_HOF', 'HOFParam.json', ...
    'output_file_CropSide', 'Cropsde_param.json', ...
    'output_file_CropFront', 'Cropfrt_param.json');

    % get params from parameter file
    params_list = importdata(params_file);
    params = params_list.data;
    params_names = params_list.textdata;
    nparams = size(params_names,1);
    params_mapping = struct(params_names{1}, params(1));
    for i=1:nparams
        name = params_names{i};
        params_mapping.(name) = params(i);
    end

    % get interest points from trx
    exp_list = importdata(exp_file);
    nexp = size(exp_list,1);
    disp(exp_list)
    %open trx file for first experiment
    trx = load(fullfile(exp_list(1), "trx.mat"));
    trx = trx.trx;
    % side ips
    trx1 = [trx(1).x(1), trx(1).y(1)];
    trx2 = [trx(2).x(1), trx(2).y(1)];
    trx3 = [trx(3).x(1), trx(3).y(1)];
    %front ips
    trx4 = [trx(4).x(1), trx(4).y(1)];
    trx5 = [trx(5).x(1), trx(5).y(1)];

    HOGParams = struct('nbins',string(params_mapping.nbins), ...
                       'cell', struct('w',string(params_mapping.cell_height),...
                                      'h', string(params_mapping.cell_height)));


    HOFParams = struct('nbins',string(params_mapping.nbins), ...
                       'cell', struct('w',string(params_mapping.cell_height),...
                                      'h', string(params_mapping.cell_height)),...
                        'lk', struct('threshold', string(params_mapping.threshold), ...
                                     'sigma',struct('derivative', string(double(params_mapping.sigma_derivative)), ...
                                                     'smoothing', string(double(params_mapping.sigma_smoothing)))));

    CropParams_side = struct('interest_pnts', struct('food', floor(trx1), ...
                              'mouth' ,floor(trx2), 'perch', floor(trx3)), ...
                              'npatches', string(params_mapping.npatches_side), ...
                              'crop_flag', string(params_mapping.crop_flag), ...
                              'ncells', string(params_mapping.ncells));
    CropParams_front = struct('interest_pnts', struct('food', floor(trx4), ...
                              'mouth' ,floor(trx5)), ...
                              'npatches', string(params_mapping.npatches_front), ...
                              'crop_flag', string(params_mapping.crop_flag), ...
                              'ncells', string(params_mapping.ncells));
    
    %%copy json params to exps
    dirname = fullfile(exp_list, json_dir_file);
    for i=1:nexp
        disp(dirname{i})    
        if ~exist(dirname{i}, 'dir')
            mkdir(dirname{i})
            saveJSONfile(HOGParams,fullfile(dirname{i},output_file_HOG));
            saveJSONfile(HOFParams,fullfile(dirname{i},output_file_HOF));
            saveJSONfile(CropParams_side,fullfile(dirname{i},output_file_CropSide));
            saveJSONfile(CropParams_front,fullfile(dirname{i},output_file_CropFront));
        else
           continue
        end
    end



end
    