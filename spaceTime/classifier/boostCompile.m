% Compile boosting code.

%% OpenMp options, set as desired (compiler specific)
% num threads will not exceed number of processors or OMP_NUM_THREADS
setenv('OMP_NUM_THREADS','128');
if(ispc), o='/openmp"'; optsOmp={'OPTIMFLAGS="$OPTIMFLAGS' o}; else
  o='-fopenmp"'; optsOmp={'CFLAGS="\$CFLAGS' o 'LDFLAGS="\$LDFLAGS' o};
end

%% compile individual source files into object files
sDir = [fileparts(mfilename('fullpath')) '/mprivate'];
if(ispc()), oExt='.obj'; else oExt='.o'; end
fs = {'Savable','Clf','ClfBoost','ClfMpl','ClfTree','ClfFtr','ClfFern',...
  'DataSetImg'};
opts=['-O' '-outdir' sDir optsOmp]; n=length(fs); os=fs;
for i=1:n, fs{i}=[sDir '/' fs{i} '.cpp']; os{i}=[sDir '/' os{i} oExt]; end
for i=1:n, disp(fs{i}); mex('-c',fs{i},opts{:}); end

%% compile boost from object files and cleanup object files
mex([sDir '/boost.cpp'],os{:},opts{:});
if(1), for i=1:n, delete(os{i}); end; end
