% test reading mouse data

%% data location

indatadir = '/groups/egnor/egnorlab/for kristin';
expname = 'b6_popcage_16_110405_09.58.30.268';

inmoviefilestr = [expname,'.seq'];
inmoviefilestr2 = [expname,'.mat'];
intrxfilestr = [expname,'t.mat'];

outrootdatadir = '/groups/branson/home/bransonk/behavioranalysis/data/roian';
outtrxfilestr = 'trx.mat';
outmoviefilestr = 'movie.seq';
outmoviefilestr2 = 'movie.mat';

%% create the experiment directory

outexpdir = fullfile(outrootdatadir,expname);
if ~exist(outexpdir,'dir'),
  mkdir(outexpdir);
end

% soft link to .seq file
inmoviename = fullfile(indatadir,inmoviefilestr);
outmoviename = fullfile(outexpdir,outmoviefilestr);
if ~exist(outmoviename,'file'),
  unix(sprintf('ln -s "%s" %s',inmoviename,outmoviename));
end

inmoviename2 = fullfile(indatadir,inmoviefilestr2);
outmoviename2 = fullfile(outexpdir,outmoviefilestr2);
if ~exist(outmoviename2,'file'),
  unix(sprintf('ln -s "%s" %s',inmoviename2,outmoviename2));
end


%% parameters
mperpx = 0.00098387;
pxpermm = 1/(1000*mperpx);
% i assigned this arbitrarily for now, change this
sex = {'F','F','M','M'};

%% convert trx file
intrxname = fullfile(indatadir,intrxfilestr);
outtrxname = fullfile(outexpdir,outtrxfilestr);
ConvertMouseTrx2Trx(intrxname,outmoviename,outtrxname,pxpermm,sex);

