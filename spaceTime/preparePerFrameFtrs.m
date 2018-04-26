function preparePerFrameFtrs(moviefilename,trackfilename,stationary,usedeep, params)

% function preparePerFrameFtrs(moviefilename,trackfilename,stationary)
% The perframe features are stored in perframedir.
% 
setup;
% method = 'hs-brightness';
% flowname = 'hs_ff';
if nargin<4,
  usedeep = false;
end

if nargin < 5,
  params = getParams;
end

if usedeep,
  method = 'deep-sup';
else
  method = 'hs-sup';
end
mndx = find(strcmp(params.methods,method));
flowname = params.flownames{mndx};
% method can be LK, hs-brightness, hs-sup. 

% method = 'LK';
% flowname = 'ff';
ftrs = computeFeaturesParallel(moviefilename,trackfilename,stationary,method, params);
expdir = fileparts(moviefilename);
extractPerframeFtrs(fullfile(expdir,'perframe'),ftrs,stationary,flowname);
  