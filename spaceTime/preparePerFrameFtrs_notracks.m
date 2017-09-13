function preparePerFrameFtrs_notracks(moviefilename,params)

%function preparePerFrameFtrs_notracks(moviefilename,params)
% The perframe features are stored in perframedir.
% 

setup;
method = 'hs-brightness';
flowname = 'ff';
ftrs = computeFeaturesParallel_notracks(moviefilename,method);
expdir = fileparts(moviefilename);
extractPerframeFtrs(fullfile(expdir,'perframe'),ftrs,false,flowname);
  