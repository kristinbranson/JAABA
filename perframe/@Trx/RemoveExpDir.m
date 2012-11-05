function RemoveExpDir(obj,expdir)

% which experiment
n = find(strcmp(expdir,obj.expdirs),1);
if isempty(n),
  error('Experiment %s not loaded',expdir);
end

% clear file names
obj.expdirs(n) = [];
obj.outexpdirs(n) = [];

% clear video info
if numel(obj.nrs) >= n,
  obj.nrs(n) = [];
  obj.ncs(n) = [];
  obj.ncolors(n) = [];
  obj.movie_nframes(n) = [];
end

flies = obj.exp2flies{n};

% clear trajectory frame info
obj.firstframes(flies) = [];
obj.endframes(flies) = [];
obj.nframes(flies) = [];

% clear arena parameters
if numel(obj.tl_x) >= max(flies),
  obj.tl_x(flies) = [];
  obj.tl_y(flies) = []; 
  obj.tr_x(flies) = []; 
  obj.tr_y(flies) = []; 
  obj.bl_x(flies) = []; 
  obj.bl_y(flies) = []; 
  obj.br_x(flies) = []; 
  obj.br_y(flies) = []; 
end

% clear fps
if numel(obj.fps) >= n,
  obj.fps(n) = [];
end

% clear roi
if numel(obj.roi) >= max(flies),
  obj.roi(flies) = [];
end

% clear trajectory files
obj.trxfiles(n) = [];

% update number of flies
obj.nflies = obj.nflies - obj.nfliespermovie(n);
obj.nfliespermovie(n) = [];

% clear indexing info
obj.fly2exp(obj.exp2flies{n}) = [];
obj.exp2flies(n) = [];

% clear conversion from pixels to mm
if numel(obj.pxpermm) >= n,
  obj.pxpermm(n) = [];
end

% clear movie size in mm
if numel(obj.width_mms) >= n,
  obj.width_mms(n) = [];
end
if numel(obj.height_mms) >= n,
  obj.height_mms(n) = [];
end

% clear cache
if numel(obj.datacached) >= n,
  obj.datacached(n) = [];
  obj.fnscached(n) = [];
  obj.nfnscached(n) = [];
end
if numel(obj.ndatacachedperexp) >= n,
  obj.ndatacached = obj.ndatacached - obj.ndatacachedperexp(n);
  obj.ndatacachedperexp(n) = [];
end

% update number of experiments
obj.nexpdirs = obj.nexpdirs - 1;
