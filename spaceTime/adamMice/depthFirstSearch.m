function m = depthFirstSearch(rootdir,fcn)
% depth first search on rootdir, executing fcn at each dir encountered:
%
% out = fcn(dir)
%
% m is a containers.Map where the keys are the dirs encountered and the
% vals are the values of out arising from execution of fcn.

assert(exist(rootdir,'dir')==7);
assert(isa(fcn,'function_handle'));

m = containers.Map();
m = lclDFS(m,rootdir,fcn);

function m = lclDFS(m,rootdir,fcn)

% this dir
try
  val = fcn(rootdir);
catch ME
  fprintf(2,'Err executing fcn in dir %s:\n',rootdir);
  fprintf(2,'%s\n',ME.getReport());
  val = [];
end
m(rootdir) = val;

% recurse
dd = dir(rootdir);
for ndx = 1:numel(dd)
  if strcmp(dd(ndx).name(1),'.'), continue; end
  
  curd = fullfile(rootdir,dd(ndx).name);
  if isdir(curd)
    m = lclDFS(m,curd,fcn);
  end  
end