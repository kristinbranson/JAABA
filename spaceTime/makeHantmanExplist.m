function makeHantmanExplist(rootdir , expdirs, varargin)

    % write explist to txt file with unix style line endings
    nexps = size(expdirs,1);
    new_exp_list = cell(nexps,1);
    fd = fopen(fullfile(rootdir,'hantman_list.txt'),'w');
    for i=1:nexps
        new_exp_list{i} = strrep(expdirs{i},'\','/');
        new_exp_list{i} = append(new_exp_list{i},'/');
        fprintf(fd, "%s\n", new_exp_list{i});
    end
    fclose(fd);
end
   