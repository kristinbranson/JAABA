function convertFeaturestoMatlab(exp_file, varargin)

    exp_list = importdata(exp_file);
    nexp = size(exp_list,1);
    disp(exp_list)
    
    for i=1:nexp
        exp = exp_list{i};
        disp(fullfile(exp,'features.mat'))
        if exist(fullfile(exp,'features.mat'),'file')
            delete(fullfile(exp,'features.mat'))
        end
        cuda2matlab_features('indir',exp);
    end
end


