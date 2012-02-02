function  launchJobs( exePathInFwHome, scriptPathOnFw, nodeMatlabName,...
  batchId, nTask, varargin )
% Launch a job on the FwGrid cluster
% 
% Each node function must be defined as:
%   function fwNode( ind, batchId, savePath )
% and end with:
%   save( [ savePath batchId sprintf('-%05d',ind) ], 'var1' );
%   exit
%
% USAGE
%  launchJobs( exePathInFwHome, scriptPathOnFw, nodeMatlabName,...
%    batchId, nTask, varargin )
%
% INPUTS
%  exePathInFwHome - path on cluster where the processes will be executed
%                    It has to be such that [ '/home/' userName '/'
%                    exePathInFwHome ] is the execution path
%  scriptPathOnFw  - path where the script will be saved temporarily
%  nodeMatlabName  - matlab routine that will be executed on each node
%  batchId         - just a string indicating the id of the batch
%  nTask           - number of total tasks to execute
%  varargin        - any variables that need to be copied to cluster for
%                    processing. Their name will be preserved
%
% Modified version of fwGridLaunch.m from:
%
% Vincent's Structure From Motion Toolbox      Version NEW
% Copyright (C) 2008 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

% we don't need to do this because we are on the same file server
%system(['sshfs ' userName '@137.110.131.21:/home/' userName ...
%  ' /projects/fwHome']);
%system(['sshfs ' userName '@137.110.131.21:/scratch/' userName ...
%  ' /projects/fwScratch']);

% Save the variables in varargin so that each node on the cluster can use
% it
while 1
  try
    for i=1:length(varargin)
      varToCopy.(inputname(5+i))=varargin{i};
    end
    save([exePathInFwHome '/' batchId '.mat'],'-struct','varToCopy');
    
    [ exePathInFwHome '/' batchId '.mat']
  catch %#ok<CTCH>
    pause(1); continue;
  end
  break
end
savePath=[ exePathInFwHome ];

temp = clock; scriptName = [ batchId '-' int2str(temp(4)) '-' ...
  int2str(temp(5)) '-' int2str(temp(6)) ];

% Create the script file
fid=fopen('launchJob.sh','w');
fprintf(fid,'#!/bin/bash \n');
fprintf(fid,'nitr=%d\n', nTask );
fprintf(fid,'matlabCommand=''%s''\n', nodeMatlabName);
fprintf(fid,'executePath=''~/%s''\n',exePathInFwHome);
fprintf(fid,'matlabPath=''/opt/matlab/bin/matlab -nosplash -nodesktop -nojvm''\n');
fprintf(fid,'\n');
fprintf(fid,'find %s -name ''*'' | xargs rm\n', savePath );
fprintf(fid,'find ~/%s -name ''%s-*.o*'' | xargs rm\n', exePathInFwHome,batchId );
fprintf(fid,'find ~/%s -name ''%s-*.e*'' | xargs rm\n', exePathInFwHome,batchId );
fprintf(fid,'find %s -name ''%s-*'' | xargs rm\n', scriptPathOnFw,batchId );
%fprintf(fid,'find ~/temp -name ''temp*'' | xargs rm\n');
fprintf(fid,'\n');
fprintf(fid,'# For each script\n');
fprintf(fid,'i=1\n');
fprintf(fid,'while [ $i -le $nitr ]\n');
fprintf(fid,'do\n');
fprintf(fid,'  # Create the script\n');
fprintf(fid,'  name=%s/%s-$i\n', scriptPathOnFw, scriptName);
fprintf(fid,['  echo "cd ${executePath}; ${matlabPath} -r \\"addpath(genpath(''~/matlab'')), $' ...
  '{matlabCommand}(${i},''%s'',''%s'')\\" " >> $name\n'], scriptName,savePath);
fprintf(fid,'  chmod 755 $name\n');
fprintf(fid,'  \n');
fprintf(fid,'  # If no possible allocation, wait and resubmit\n');
fprintf(fid,'  while true\n');
fprintf(fid,'  do\n');
fprintf(fid,'    # Check the number of running MATLAB licenses\n');
fprintf(fid,'    nLicense=$(/opt/matlab/etc/lmstat | grep -o "Total of [0-9]* license[s]* in use" | grep -o "[0-9][0-9]*")\n');
fprintf(fid,'    \n');
fprintf(fid,'    if [ $nLicense -gt 45 ]\n');
fprintf(fid,'    then\n');
fprintf(fid,'      sleep 2\n');
fprintf(fid,'      continue\n');
fprintf(fid,'    fi\n');
fprintf(fid,'    \n');
fprintf(fid,'    # Check if the script can be submitted properly\n');
fprintf(fid,'    test=`qsub -l matlab=1 -l arch=lx26-x86 -cwd $name 2>&1 | grep -c submitted`\n');
fprintf(fid,'    if [ $test -gt 0 ]\n');
fprintf(fid,'    then\n');
fprintf(fid,'      sleep 5\n');
fprintf(fid,'      break\n');
fprintf(fid,'    fi\n');
fprintf(fid,'    \n');
fprintf(fid,'    sleep 5\n');
fprintf(fid,'  done\n');
fprintf(fid,'  ((i=i+1))\n');
fprintf(fid,'done\n');

fclose(fid);

% Make script executable
system('chmod 755 launchJob.sh');
% Put the script in the directory we will be executing in
system(['mv launchJob.sh ' exePathInFwHome '/.']);
% Touch it 
system(['touch ' exePathInFwHome '/launchJob.sh']);
% Launch it
%system(['cd ' exePathInFwHome '; ./launchJob.sh &']);

% Reconnect to the ssh shares
%system(['sshfs ' userName '@137.110.131.21:/home/' userName ...
%  ' /projects/fwHome']);
%system(['sshfs ' userName '@137.110.131.21:/scratch/' userName ...
%  ' /projects/fwScratch']);

