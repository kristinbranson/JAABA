% test reading mouse data

%% data locations

% path to root directory containing the input data
inmoviedir = 'Y:\popcage_enriched\mousetrack_18';
intrxdir = 'Y:\adam\mousetrack_18\Results\Tracks';
% base name of the input experiment
inmovieexpname = 'b6_popcage_18_09.15.11_10.56.24.135';
intrxexpname = 'b6_popcage_18_2011.09.15_10.56.24.135';

% path to root directory to hold the output data
outdatadir = 'C:\Workspace\data';

%% parameters

% scale
mperpx = 0.00098387;
pxpermm = 1/(1000*mperpx);

% i assigned this arbitrarily for now
sex = {'F','M','M','F'};

% whether to make links or to copy
makelinks = true;

% frame interval length
frameintervallength = 108000;

% read in length of movie
inseqfile = fullfile(inmoviedir,[inmovieexpname,'.seq']);
headerinfo = r_readseqinfo(inseqfile);
nframes = headerinfo.m_iNumFrames;


%% main function call

for intervalstart = 1:frameintervallength:nframes,
  intervalend = intervalstart + frameintervallength - 1;
  if intervalend > nframes,
    break;
  end
  frameinterval = [intervalstart,intervalend];
  fprintf('Creating experiment for frames %d to %d...\n',intervalstart,intervalend);
  ConvertMouseTrx2Trx(inmoviedir,intrxdir,outdatadir,inmovieexpname,intrxexpname,...
    pxpermm,sex,'frameinterval',frameinterval,'makelinks',makelinks);
end
