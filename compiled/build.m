function build()

% get the directory where this file lives
% figure out where the root of the Ohayon code is
thisScriptFileName=mfilename('fullpath');
thisScriptDirName=fileparts(thisScriptFileName);
thisScriptDirParts=split_on_filesep(thisScriptDirName);
  % a cell array with each dir an element
jaabaRootParts=thisScriptDirParts(1:end-1);
jaabaRootDirName=combine_with_filesep(jaabaRootParts);

% just put the executable in the did with the build script
exeDirName=thisScriptDirName;

% Invoke mcc to compile the MouseTrackProj executable.
% Doing this via mcc means we only use the Matlab Compiler license
% for the duration of the compile
% cmdLine= ...
%   sprintf(['mcc ' ...
%            '-o JAABA ' ...
%            '-m ' ...
%            '-d "%s" ' ....
%            '-I "%s/filehandling" ' ...
%            '-I "%s/misc" ' ...
%            '-I "%s/perframe" ' ...
%            '-I "%s/perframe/compute_perframe_features" ' ...
%            '-v ' ...
%            '"%s/perframe/StartJAABA_compiled.m" ' ...
%            '-a "%s/perframe/janelia_logo.png" ' ...
%            '-a "%s/perframe/params"'], ...
%           exeDirName, ...
%           jaabaRootDirName, ...
%           jaabaRootDirName, ...
%           jaabaRootDirName, ...
%           jaabaRootDirName, ...
%           jaabaRootDirName, ...
%           jaabaRootDirName, ...
%           jaabaRootDirName);
% %            '-w enable:specified_file_mismatch ' ...
% %            '-w enable:repeated_file ' ...
% %            '-w enable:switch_ignored ' ...
% %            '-w enable:missing_lib_sentinel ' ...
% %            '-w enable:demo_license ' ...        
% %           '-W main:JAABA -T link:exe ' ...
% % we have to -a perframe/params just to get
% % perframe/params/featureConfig.xml
% fprintf('%s\n',cmdLine);
% system(cmdLine);
fprintf('Invoking mcc...\n');
mcc('-o','JAABA', ...
    '-m', ...
    '-d',exeDirName, ...
    '-I',fullfile(jaabaRootDirName,'filehandling'), ...
    '-I',fullfile(jaabaRootDirName,'misc'), ...
    '-I',fullfile(jaabaRootDirName,'perframe'), ...
    '-I',fullfile(jaabaRootDirName,'perframe','compute_perframe_features'), ...
    '-v', ...
    fullfile(jaabaRootDirName,'perframe','StartJAABA_compiled.m'), ...
    '-a',fullfile(jaabaRootDirName,'perframe','janelia_logo.png'), ...
    '-a',fullfile(jaabaRootDirName,'perframe','version.txt'), ...
    '-a',fullfile(jaabaRootDirName,'perframe','params'));

fprintf('Clearing out intermediate/useless files...\n');

file_name=fullfile(exeDirName,'mccExcludedFiles.log');
if exist(file_name,'file')
  delete(file_name);
end

file_name=fullfile(exeDirName,'readme.txt');
if exist(file_name,'file')
  delete(file_name);
end

fprintf('Done building.\n');

end





function path=combine_with_filesep(path_as_array)
  % combine a cell array of dir names into a single path name
  n=length(path_as_array);
  if n>0
    path=path_as_array{1};
    for i=2:n
      path=[path filesep path_as_array{i}];  %#ok
    end
  end
end





function path_as_array=split_on_filesep(path)
  % split a path on filesep into a cell array of single dir names
  
  % check for empty input
  if isempty(path)
    path_as_array=cell(0,1);
    return;
  end
  
  % trim a trailing fileseparator
  if path(end)==filesep
    path=path(1:end-1);
  end
  
  i_pathsep=strfind(path,filesep);
  n=length(i_pathsep)+1;
  path_as_array=cell(n,1);
  if n>0
    if n==1
      path_as_array{1}=path;
    else
      % if here, n>=2
      path_as_array{1}=path(1:i_pathsep(1)-1);
      for i=2:(n-1)
        path_as_array{i}=path(i_pathsep(i-1)+1:i_pathsep(i)-1);
      end
      path_as_array{n}=path(i_pathsep(n-1)+1:end);
    end
  end
end

