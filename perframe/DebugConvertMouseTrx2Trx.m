% test reading mouse data

%% data location

indatadir = '/groups/egnor/egnorlab/for kristin';
%inexpname = 'b6_popcage_16_110405_09.58.30.268';
%inexpname = 'b6_popcage_16_110405_09.58.30.268.108000_216000';
inexpname = 'b6_popcage_16_110405_09.58.30.268.216000_324000';
%inexpname_seq = inexpname;
%inexpname_seq = 'b6_popcage_16_110405_09.58.30.268 - frames 108000-216000';
inexpname_seq = 'b6_popcage_16_110405_09.58.30.268 - frames 216000-324000';
outexpname = regexprep(inexpname,'[^\w.]+','_');
if ~strcmp(inexpname,outexpname),
  fprintf('Converting %s to %s\n',inexpname,outexpname);
end

inmoviefilestr = [inexpname_seq,'.seq'];
inmoviefilestr2 = [inexpname,'.mat'];
intrxfilestr = [inexpname,'t.mat'];

outrootdatadir = '/groups/branson/home/bransonk/behavioranalysis/data/roian';
outtrxfilestr = 'trx.mat';
outmoviefilestr = 'movie.seq';
outmoviefilestr2 = 'movie.mat';

%% create the experiment directory

outexpdir = fullfile(outrootdatadir,outexpname);
if ~exist(outexpdir,'dir'),
  mkdir(outexpdir);
end

% soft link to .seq file
inmoviename = fullfile(indatadir,inmoviefilestr);
outmoviename = fullfile(outexpdir,outmoviefilestr);
if ~exist(inmoviename,'file'),
  error('File %s does not exist',inmoviename);
end
if ~exist(outmoviename,'file'),
  unix(sprintf('ln -s "%s" "%s"',inmoviename,outmoviename));
end

inmoviename2 = fullfile(indatadir,inmoviefilestr2);
outmoviename2 = fullfile(outexpdir,outmoviefilestr2);
if ~exist(inmoviename2,'file'),
  error('File %s does not exist',inmoviename2);
end
if ~exist(outmoviename2,'file'),
  unix(sprintf('ln -s "%s" "%s"',inmoviename2,outmoviename2));
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

