% test reading mouse data

%% data locations

% path to root directory containing the input data
indatadir = '/groups/egnor/egnorlab/for kristin';
% base name of the input experiment
inexpname = 'b6_popcage_16_110405_09.58.30.268.216000_324000';
% input .seq file name; this will be used for all intervals of frames
inexpname_seq = 'b6_popcage_16_110405_09.58.30.268.seq';

% path to root directory to hold the output data
outdatadir = '/groups/branson/home/bransonk/behavioranalysis/data/roian';

%% parameters

% scale
mperpx = 0.00098387;
pxpermm = 1/(1000*mperpx);

% i assigned this arbitrarily for now
sex = {'F','F','M','M'};

% whether to make links or to copy
makelinks = true;

%% main function call

ConvertMouseTrx2Trx(indatadir,outdatadir,inexpname,inexpname_seq,pxpermm,sex,'makelinks',makelinks);
