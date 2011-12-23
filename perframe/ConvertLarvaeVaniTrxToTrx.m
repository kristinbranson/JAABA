function ConvertLarvaeVaniTrxToTrx(expdir,varargin)

% parse optional parameters
trxfilestr = 'JLabel_trx.mat';
outmoviefilestr = 'JLabel_movie.avi';
perframedirstr = 'perframe';

[trxfilestr,outmoviefilestr,perframedirstr] = ...
  myparse(varargin,...
  'trxfilestr',trxfilestr,...
  'outmoviefilestr',outmoviefilestr,...
  'perframedirstr',perframedirstr);

% get experiment_name from full path
[~,experiment_name] = myfileparts(expdir);

%% input file names

inmoviename = fullfile(expdir,sprintf('%s_video.avi',experiment_name));
indataname = fullfile(expdir,sprintf('%s_data.txt',experiment_name));

%% output file names

outtrxname = fullfile(expdir,trxfilestr);
outmoviename = fullfile(expdir,outmoviefilestr);
perframedir = fullfile(expdir,perframedirstr);

%% convert trx

%indata = importdata(indataname);
indata.data = [];
fid = fopen(indataname,'r');
if fid < 0,
  error('Could not open file %s for reading',indataname);
end
while true,
  s = fgetl(fid);
  if ~ischar(s),
    break;
  end
  ss = strsplit(s,',');
  ss = str2num(char(ss(1:17))); %#ok<ST2NM>
  indata.data(:,end+1) = ss;
end
%indata.data = [(0:size(indata.data,1)-1)',indata.data];

timestamps = indata.data(2,:);

trx = struct;
trx.x = indata.data(10,:);
trx.y = indata.data(11,:);
trx.skeletonx = indata.data([3,5,7],:);
trx.skeletony = indata.data([4,6,8],:);
trx.a = indata.data(9,:) / 4;
trx.b = trx.a;
trx.theta = zeros(size(trx.a));
trx.dt = diff(timestamps);
trx.id = 1; % only one larva
trx.moviename = outmoviename;
trx.firstframe = 1;
trx.endframe = numel(trx.x);
trx.nframes = numel(trx.x);
trx.off = 0;
trx.timestamps = timestamps;

save(outtrxname,'trx','timestamps');

%% copy/soft link movie

% this will only work on linux and mac
if ispc,
  copyfile(inmoviename,outmoviename);
else
  system(sprintf('ln -s %s %s',inmoviename,outmoviename));
end

%% per-frame statistics

% create the directory
if ~exist(perframedir,'dir'),
  mkdir(perframedir);
end

% create the per frame statistics

% change in head-to-body angle
fn = 'dhead2body';
dangle = modrange(diff(indata.data(12,:)),-180,180);
data = {dangle ./ trx.dt }; %#ok<NASGU>
units = struct('num','deg','den','s'); %#ok<NASGU>
save(fullfile(perframedir,[fn,'.mat']),'data','units');

% ADD MORE HERE
