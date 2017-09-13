function prepareDeepPerFrameFtrsCompiled(moviefilename,trackfilename,stationary,MCRloc,tempdir,maxcores)

% function prepareDeepPerFrameFtrsCompiled(moviefilename,trackfilename,stationary,maxcores)
setup;
usedeep = true;
method = 'deep-sup';
flowname = 'DS';

if ~exist(maxcores,'var')
  totcores = features('numcores');
  if totcores>10,
    maxcores = totcores-4;
  else
    maxcores = totcores-1;
  end
end
gparams = getParams;
[readfcn,nframes,fid,headerinfo] = get_readframe_fcn(moviename);

nblocks = ceil(nframes/gparams.blocksize);
expdir = fileparts(moviefilename);
[~,expname] = fileparts(expdir);
tempname = fullfile(tempdir,expname);
partcmd = sprintf('%s %s %s %s %d deep-sup ',...
  'computeFeaturesCompiled/for_testing/run_computeFeaturesCompiled.sh',...
  MCRloc,moviefilename,trackfilename,stationary);
