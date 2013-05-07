function [params,behs]=load_parameters(model_file)
    str = fileread(model_file);
    params = loadjson(str);
    behs = {};
    for ii=1:length(params.behaviors),
        behs{ii} = params.behaviors(ii).name;
    end
end

