function [data,units] = compute_spacetime_hoghof(trx,n,m_spacetime)

params = trx.stInfo;
method = params.cur_method;
mndx = find(strcmp(params.methods,method));
flowname = params.flownames{mndx};
expdir = trx.expdirs{n};
movstr= trx.moviefilestr;
trxstr = trx.trxfilestr;
perframedir = trx.perframedir;
moviefilename = fullfile(expdir, movstr);
trackfilename = fullfile(expdir, trxstr);
ftrs = computeSTFeaturesParallel(moviefilename,trackfilename,params.is_stationary,method, params);
extractSTPerframeFtrs(fullfile(expdir,perframedir),ftrs,params.is_stationary,flowname,params);

aa = load(fullfile(expdir,perframedir,sprintf('st_%s.mat',m_spacetime)));
data = aa.data;
units = aa.units;