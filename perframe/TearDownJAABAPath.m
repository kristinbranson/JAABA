function TearDownJAABAPath()
    % Tear down the JAABA path
    % Currently (Feb 2021) only used to avoid name collisions with FlyBowlAnalysis compute_*
    % functions.

    original_warning_state = warning('query', 'MATLAB:rmpath:DirNotFound') ;
    warning('off', 'MATLAB:rmpath:DirNotFound') ;
    
    jlabelpath = fileparts(mfilename('fullpath'));
    % Initialize all the paths.
    %rmpath(jlabelpath);  % don't rm b/c need for FBA
    baseDir = fileparts(jlabelpath);
    %rmpath(fullfile(baseDir,'misc')); % don't rm b/c need for FBA
    %rmpath(fullfile(baseDir,'filehandling')); % don't rm b/c need for FBA
    rmpath(fullfile(jlabelpath,'larva_compute_perframe_features'));
    rmpath(fullfile(jlabelpath,'compute_perframe_features'));
    rmpath(fullfile(baseDir,'perframe','params'));
    rmpath(fullfile(baseDir,'tests'));
    st_dir = fullfile(baseDir,'spaceTime'); 
    rmpath(st_dir);
    rmpath(genpath(st_dir));
    %rmpath(genpath(st_dir));
    
    warning(original_warning_state) ;
end
