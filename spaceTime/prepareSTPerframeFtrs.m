function prepareSTPerframeFtrs(moviefilename,trackfilename,stationary,usedeep, params,perframedir)

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
ftrs = computeSTFeaturesParallel(moviefilename,trackfilename,stationary,method, params);
expdir = fileparts(moviefilename);
extractSTPerframeFtrs(fullfile(expdir,perframedir),ftrs,stationary,flowname,params);
  