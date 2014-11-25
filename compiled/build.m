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
    '-I',fullfile(jaabaRootDirName,'perframe','larva_compute_perframe_features'), ...
    '-v', ...
    fullfile(jaabaRootDirName,'perframe','StartJAABA_compiled.m'), ...
    '-a',fullfile(jaabaRootDirName,'perframe','janelia_logo.png'), ...
    '-a',fullfile(jaabaRootDirName,'perframe','version.txt'), ...
    '-a',fullfile(jaabaRootDirName,'perframe','params'),...
    '-a',fullfile(jaabaRootDirName,'perframe','compute_perframe_features'), ...
    '-a',fullfile(jaabaRootDirName,'perframe','larva_compute_perframe_features'), ...
    '-R','-logfile,JAABA.log',...
    '-R','-startmsg,"Starting JAABA. Please wait, this may take a while."'...
    );

fprintf('Clearing out intermediate/useless files...\n');

file_name=fullfile(exeDirName,'mccExcludedFiles.log');
if exist(file_name,'file')
  delete(file_name);
end

file_name=fullfile(exeDirName,'readme.txt');
if exist(file_name,'file')
  delete(file_name);
end

% make the JAABA app executable for mac
if ismac,
  ff = fopen('JAABA.app/Contents/MacOS/preprelaunch','w');
  fprintf(ff,'#!/bin/bash\n');
  fprintf(ff,'cd `dirname "$0"` \n');
  fprintf(ff,'cd ../../..\n');
  fprintf(ff,'exe_dir=`pwd -P`\n');
  fprintf(ff,'echo exe_dir is ${exe_dir};\n');
  fprintf(ff,'export JAABAROOT=$exe_dir\n');
  fprintf(ff,'exec "${exe_dir}"/JAABA.app/Contents/MacOS/prelaunch\n');
  fprintf(ff,'exit');
  cmd = 'sed -i "" ''s/prelaunch/preprelaunch/g'' JAABA.app/Contents/Info.plist';
  system(cmd);
  fclose(ff);
  fileattrib('JAABA.app/Contents/MacOS/preprelaunch','+x');

end



% JAABAPlot
mcc('-o','JAABAPlot', ...
    '-m', ...
    '-d',exeDirName, ...
    '-I',fullfile(jaabaRootDirName,'filehandling'), ...
    '-I',fullfile(jaabaRootDirName,'misc'), ...
    '-I',fullfile(jaabaRootDirName,'perframe'), ...
    '-I',fullfile(jaabaRootDirName,'perframe','compute_perframe_features'), ...
    '-I',fullfile(jaabaRootDirName,'perframe','larva_compute_perframe_features'), ...
    '-v', ...
    fullfile(jaabaRootDirName,'perframe','JAABAPlot_compiled.m'), ...
    '-a',fullfile(jaabaRootDirName,'perframe','janelia_logo.png'), ...
    '-a',fullfile(jaabaRootDirName,'perframe','version.txt'), ...
    '-a',fullfile(jaabaRootDirName,'misc','javasysmon-0.3.4.jar'),...
    '-R','-logfile,JAABAPlot.log',...
    '-R','-startmsg,"Starting JAABAPlot. Please wait, this may take a while."'...
    );

  
fprintf('Clearing out intermediate/useless files...\n');

file_name=fullfile(exeDirName,'mccExcludedFiles.log');
if exist(file_name,'file')
  delete(file_name);
end

file_name=fullfile(exeDirName,'readme.txt');
if exist(file_name,'file')
  delete(file_name);
end
if ismac,
  ff = fopen('JAABAPlot.app/Contents/MacOS/preprelaunch','w');
  fprintf(ff,'#!/bin/bash\n');
  fprintf(ff,'cd `dirname "$0"` \n');
  fprintf(ff,'cd ../../..\n');
  fprintf(ff,'exe_dir=`pwd -P`\n');
  fprintf(ff,'echo exe_dir is ${exe_dir};\n');
  fprintf(ff,'export JAABAROOT=$exe_dir\n');
  fprintf(ff,'exec "${exe_dir}"/JAABAPlot.app/Contents/MacOS/prelaunch\n');
  fprintf(ff,'exit');
  fclose(ff);
  cmd = 'sed -i "" ''s/prelaunch/preprelaunch/g'' JAABAPlot.app/Contents/Info.plist';
  system(cmd);
  fileattrib('JAABAPlot.app/Contents/MacOS/preprelaunch','+x');
end


% PrepareJAABAData
mcc('-o','PrepareJAABAData', ...
    '-m', ...
    '-d',exeDirName, ...
    '-I',fullfile(jaabaRootDirName,'filehandling'), ...
    '-I',fullfile(jaabaRootDirName,'misc'), ...
    '-I',fullfile(jaabaRootDirName,'perframe'), ...
    '-v', ...
    fullfile(jaabaRootDirName,'perframe','PrepareJAABAData_compiled.m'), ...
    '-a',fullfile(jaabaRootDirName,'perframe','janelia_logo.png'), ...
    '-a',fullfile(jaabaRootDirName,'perframe','version.txt'), ...
    '-R','-logfile,PrepareJAABAData.log',...
    '-R','-startmsg,"Starting PrepareJAABAData. Please wait, this may take a while."'...
    );

  
fprintf('Clearing out intermediate/useless files...\n');

file_name=fullfile(exeDirName,'mccExcludedFiles.log');
if exist(file_name,'file')
  delete(file_name);
end

file_name=fullfile(exeDirName,'readme.txt');
if exist(file_name,'file')
  delete(file_name);
end

if ismac,
  ff = fopen('PrepareJAABAData.app/Contents/MacOS/preprelaunch','w');
  fprintf(ff,'#!/bin/bash\n');
  fprintf(ff,'cd `dirname "$0"` \n');
  fprintf(ff,'cd ../../..\n');
  fprintf(ff,'exe_dir=`pwd -P`\n');
  fprintf(ff,'echo exe_dir is ${exe_dir};\n');
  fprintf(ff,'export JAABAROOT=$exe_dir\n');
  fprintf(ff,'exec "${exe_dir}"/PrepareJAABAData.app/Contents/MacOS/prelaunch\n');
  fprintf(ff,'exit');
  fclose(ff);
  cmd = 'sed -i "" ''s/prelaunch/preprelaunch/g'' PrepareJAABAData.app/Contents/Info.plist';
  system(cmd);
  fileattrib('PrepareJAABAData.app/Contents/MacOS/preprelaunch','+x');

end


if ismac,
  vid = fopen(fullfile(jaabaRootDirName,'perframe','version.txt'),'r');
  vv = textscan(vid,'%s');
  fclose(vid);
  vv = vv{1}{1};
  outdirname = sprintf('JAABA_MAC_%s',vv);
  if ~exist(outdirname,'dir')
    mkdir(outdirname);
  end
  movefile('JAABA.app',outdirname,'f');
  movefile('JAABAPlot.app',outdirname,'f');
  movefile('PrepareJAABAData.app',outdirname,'f');
  copyfile(fullfile(jaabaRootDirName,'LICENSE.txt'),outdirname);
  copyfile(fullfile(jaabaRootDirName,'README.txt'),outdirname);
  copyfile(fullfile(jaabaRootDirName,'misc','javasysmon-0.3.4.jar'),outdirname);
  
  
end

if isunix
  vid = fopen(fullfile(jaabaRootDirName,'perframe','version.txt'),'r');
  vv = textscan(vid,'%s');
  fclose(vid);
  vv = vv{1}{1};
  outdirname = sprintf('JAABA_LINUX_%s',vv);
  if ~exist(outdirname,'dir')
    mkdir(outdirname);
  end
  copyfile('JAABA',outdirname);
  copyfile('JAABAPlot',outdirname);
  copyfile('PrepareJAABAData',outdirname);
  copyfile('run_JAABA.sh',outdirname);
  copyfile('run_JAABAPlot.sh',outdirname);
  copyfile('run_PrepareJAABAData.sh',outdirname);
  copyfile(fullfile(jaabaRootDirName,'LICENSE.txt'),outdirname);
  copyfile(fullfile(jaabaRootDirName,'README.txt'),outdirname);
  copyfile(fullfile(jaabaRootDirName,'misc','javasysmon-0.3.4.jar'),outdirname);
  copyfile(fullfile(jaabaRootDirName,'perframe','JAABAParCompProfile.settings'),outdirname);
  
  
end

fprintf('Done Building\n');


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

