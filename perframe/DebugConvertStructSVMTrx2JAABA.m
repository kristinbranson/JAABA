%% set up path

addpath ../misc;
addpath ../filehandling;

%% parameters

rootdatadir = 'C:\Code\Jdetect\fly_data.sharpTurns\TrainingData';
expdir = 'test';

%% grab files

trxfiles = dir(fullfile(rootdatadir,'*.trx'));
trxfiles = {trxfiles.name};
trxm = regexp(trxfiles,'^fly(?<id>\d+)_(?<tag>.*)\.trx$','names');
ism = ~cellfun(@isempty,trxm);
trxfiles = trxfiles(ism);
trxm = trxm(ism);
trxm = cell2mat(trxm);
[tags,~,tagidx] = unique({trxm.tag});
counts = hist(tagidx,1:numel(tags));
[~,tagi] = max(counts);
tag = tags{tagi};
trxm(tagidx~=tagi) = [];
trxfiles(tagidx~=tagi) = [];
ids = str2double({trxm.id});
[ids,order] = sort(ids);
trxfiles = trxfiles(order);

labelfiles = dir(fullfile(rootdatadir,'*.label'));
labelfiles = {labelfiles.name};
labelm = regexp(labelfiles,'^fly(?<id>\d+)_(?<tag>.*)\.trx\.label$','names');
ism = ~cellfun(@isempty,labelm);
labelfiles = labelfiles(ism);
labelm = labelm(ism);
labelm = cell2mat(labelm);
idx = strcmp({labelm.tag},tag);
labelm = labelm(idx);
labelfiles = labelfiles(idx);
labelids = str2double({labelm.id});
[ism,idx] = ismember(ids,labelids);
trxfiles(~ism) = [];
labelfiles = labelfiles(idx(ism));

trxfiles = cellfun(@(x) fullfile(rootdatadir,x),trxfiles,'UniformOutput',false);
labelfiles = cellfun(@(x) fullfile(rootdatadir,x),labelfiles,'UniformOutput',false);

%%

ConvertStructSVMTrx2JAABA(trxfiles,labelfiles,expdir);
%%

rootdatadir = 'C:\Data\JAABA\eyrun\midres_flies_cleanedmovies';
expdir = 'test_movie1seq1';

%% grab files

trxfiles = {fullfile(rootdatadir,'movie1_seq1.trx')};
labelfiles = {fullfile(rootdatadir,'movie1_seq1.label')};

%%
ConvertStructSVMTrx2JAABA(trxfiles,labelfiles,expdir);
