function [sshcmd,bsubcmd,info] = wrapjobcmd_bsub_ssh(jobcmd,varargin)

[info.jobname,info.resfile,info.reserrfile,info.ncores] = myparse(varargin,...
  'jobname',sprintf('job%s',datestr(now,'yyyymmddTHHMMSS')),...
  'resfile',[tempname,'.out'],...
  'reserrfile',[tempname,'.err'],...
  'ncores',1);

bsubcmd = sprintf('bsub -n %d -J %s -o %s -e %s "%s"',info.ncores,info.jobname,info.resfile,info.reserrfile,strrep(jobcmd,'"','\"'));
sshcmd = sprintf('ssh login1 "%s"',strrep(strrep(bsubcmd,'\','\\'),'"','\"'));
