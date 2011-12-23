%% set up path

[~,computername] = system('hostname');
computername = strtrim(computername);

switch computername,
  
  case 'bransonk-desktop',

    JCtrax_path = '/groups/branson/home/bransonk/tracking/code/JCtrax';
    configfilename = 'params/JLabelParams_bransonk-desktop.xml';
    rootdatadir = '/groups/branson/home/bransonk/behavioranalysis/data/larvae_vani';
    
  case 'bransonk-lw2',

    JCtrax_path = 'C:\Code\JCtrax';
    configfilename = 'params\JLabelParamsKristinChase.xml';
    rootdatadir = 'C:\Data\larvae_vani';

  otherwise
    
    warning('Unknown computer %s. Paths may not be setup correctly',computername);
    JCtrax_path = '/groups/branson/home/bransonk/tracking/code/JCtrax';
    configfilename = 'params/JLabelParams_bransonk-desktop.xml';
    rootdatadir = '/groups/branson/home/bransonk/behavioranalysis/data/larvae_vani';

end

addpath(genpath(fullfile(JCtrax_path,'pdollar_toolbox')));
addpath(fullfile(JCtrax_path,'misc'));
addpath(fullfile(JCtrax_path,'filehandling'));
jlabelpath = fileparts(which('JLabel'));
addpath(fullfile(jlabelpath,'compute_perframe_features'));


%% parameters

experiment_name = '20111213-144536';
trxfilestr = 'JLabel_trx.mat';
outmoviefilestr = 'JLabel_movie.avi';
perframedirstr = 'perframe';


%% input file names

inmoviename = fullfile(rootdatadir,experiment_name,sprintf('%s_video.avi',experiment_name));
indataname = fullfile(rootdatadir,experiment_name,sprintf('%s_data.txt',experiment_name));

%% output file names

outtrxname = fullfile(rootdatadir,experiment_name,trxfilestr);
outmoviename = fullfile(rootdatadir,experiment_name,outmoviefilestr);
perframedir = fullfile(rootdatadir,experiment_name,perframedirstr);

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
  ss = str2num(char(ss(1:17)));
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
data = {dangle ./ trx.dt };
units = struct('num','deg','den','s');
save(fullfile(perframedir,[fn,'.mat']),'data','units');

% ADD MORE HERE

%% DEBUG: plot the annotated movie

for i = 1:1000,
set(him,'CData',flipdim(readframe(i),1));
set(hskel,'XData',trx.skeletonx(:,i),'YData',trx.skeletony(:,i));
set(hcentroid,'XData',trx.x(i),'ydata',trx.y(i));
pause(.1);
end