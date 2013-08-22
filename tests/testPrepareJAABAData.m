function success = testPrepareJAABAData(datatype,rootdir)

success = false;
if nargin < 2,
  rootdir = '/groups/branson/bransonlab/projects/JAABA/sampledata';
end

%% create the JAABA GUI
hfig = PrepareJAABAData;

% get handles
handles = guidata(hfig);

%% parameters

params = struct;
params.datatype = datatype;

switch datatype,
  
  case 'Ctrax'

    params.inputdir = fullfile(rootdir,'flies','GMR_14C07_AE_01_TrpA_Rig1Plate15BowlA_20120404T141155');
    params.inputmoviefilestr = 'movie.ufmf';
    params.inputmoviefile = fullfile(params.inputdir,params.inputmoviefilestr);
    params.inputfiledict = {'intrxfile','ctrax_results.mat',fullfile(params.inputdir,'ctrax_results.mat')
      'annfile','movie.ufmf.ann',fullfile(params.inputdir,'movie.ufmf.ann')};
    [~,params.experiment_name] = fileparts(params.inputdir);
    params.experiment_name = [params.experiment_name,'__Ctrax'];
    params.outputdir = '../tests/PrepareJAABADataTests';
    params.outputmoviefilestr = params.inputmoviefilestr;
    params.outputtrxfilestr = 'registered_trx.mat';
    params.outputperframedir = 'perframe';
    
    % sample jab file
    params.jabfile = fullfile(rootdir,'flies','ChaseAR_new_v10p0.jab');
    
    % options
    params.options = struct;
    params.options.SoftLinkFiles = true;
    params.options.fliplr = false;
    params.options.flipud = false;
    params.options.dotransposeimage = false;
    params.options.fps = 30.344385106161738;
    params.options.pxpermm = 7.888175824175824;
    params.options.OverRideFPS = true;
    params.options.OverRideArena = false;
    params.options.CropFirstFrame = 1;
    params.options.CropEndFrame = inf;
    params.options.ArenaType = 'circle';
    
    % read arena or fps?
    params.readarena = true;
    params.readfps = true;
    
  case 'CtraxPlusWings',

    params.inputdir = fullfile(rootdir,'flies','GMR_14C07_AE_01_TrpA_Rig1Plate15BowlA_20120404T141155');
    params.inputmoviefilestr = 'movie.ufmf';
    params.inputmoviefile = fullfile(params.inputdir,params.inputmoviefilestr);
    params.inputfiledict = {'intrxfile','wingtracking_results.mat',fullfile(params.inputdir,'wingtracking_results.mat')
      'annfile','movie.ufmf.ann',fullfile(params.inputdir,'movie.ufmf.ann')
      'inperframedir','perframe',fullfile(params.inputdir,'perframe')};
    [~,params.experiment_name] = fileparts(params.inputdir);
    params.experiment_name = [params.experiment_name,'__CtraxPlusWings'];
    params.outputdir = '../tests/PrepareJAABADataTests';
    params.outputmoviefilestr = params.inputmoviefilestr;
    params.outputtrxfilestr = 'wingtracking_results.mat';
    params.outputperframedir = 'perframe';
    
    % sample jab file
    params.jabfile = fullfile(rootdir,'flies','WingExtension_AR_v1p0.jab');
    
    % options
    params.options = struct;
    params.options.SoftLinkFiles = true;
    params.options.fliplr = false;
    params.options.flipud = false;
    params.options.dotransposeimage = false;
    params.options.fps = 30.344385106161738;
    params.options.pxpermm = 7.888175824175824;
    params.options.OverRideFPS = true;
    params.options.OverRideArena = false;
    params.options.CropFirstFrame = 1;
    params.options.CropEndFrame = inf;
    params.options.ArenaType = 'circle';

    % read arena or fps?
    params.readarena = true;
    params.readfps = true;
 
  case 'SimpleTwoFlies',

    params.inputdir = fullfile(rootdir,'courtship_bowls','FCF_pBDPGAL4U_1500437_TrpA_Rig2Plate17BowlD_20121121T152832');
    params.inputmoviefilestr = 'movie.ufmf';
    params.inputmoviefile = fullfile(params.inputdir,params.inputmoviefilestr);
    params.inputfiledict = {'intrxfile','registered_trx.mat',fullfile(params.inputdir,'registered_trx.mat')
      'inperframedir','perframe',fullfile(params.inputdir,'perframe')};
    [~,params.experiment_name] = fileparts(params.inputdir);
    params.experiment_name = [params.experiment_name,'__SimpleTwoFlies'];
    params.outputdir = '../tests/PrepareJAABADataTests';
    params.outputmoviefilestr = params.inputmoviefilestr;
    params.outputtrxfilestr = 'registered_trx.mat';
    params.outputperframedir = 'perframe';
    
    % sample jab file
    params.jabfile = fullfile(rootdir,'courtship_bowls','copulation_cb_ARv1p8.jab');
    
    % options
    params.options = struct;
    params.options.SoftLinkFiles = true;
    params.options.fliplr = false;
    params.options.flipud = false;
    params.options.dotransposeimage = false;
    params.options.OverRideFPS = true;
    params.options.OverRideArena = false;
    params.options.CropFirstFrame = 1;
    params.options.CropEndFrame = inf;
    params.options.ArenaType = 'None';

    % read arena or fps?
    params.readarena = false;
    params.readfps = true;

    
  case 'LarvaeRiveraAlba'
    
    params.inputdir = fullfile(rootdir,'larvae','OregonRA');
    params.inputmoviefilestr = 'movie.mov';
    params.inputmoviefile = fullfile(params.inputdir,params.inputmoviefilestr);
    params.inputfiledict = {'intrxfile','trx.mat',fullfile(params.inputdir,'trx.mat')
       'inperframedir','perframe',fullfile(params.inputdir,'perframe')};
    [~,params.experiment_name] = fileparts(params.inputdir);
    params.experiment_name = [params.experiment_name,'__LarvaeRiveraAlba'];
    params.outputdir = '../tests/PrepareJAABADataTests';
    params.outputmoviefilestr = params.inputmoviefilestr;
    params.outputtrxfilestr = 'trx.mat';
    params.outputperframedir = 'perframe';
    
    % sample jab file
    params.jabfile = fullfile(rootdir,'larvae','TODO.jab');
    
    
    % options
    params.options = struct;
    params.options.SoftLinkFiles = true;
    params.options.fliplr = false;
    params.options.flipud = false;
    params.options.dotransposeimage = false;
    params.options.OverRideFPS = false;
    params.options.OverRideArena = false;
    params.options.CropFirstFrame = 1;
    params.options.CropEndFrame = inf;
    params.options.ArenaType = 'circle';    
    
    % read arena or fps?
    params.readarena = false;
    params.readfps = false;

  case 'Qtrax',
    
    params.inputdir = fullfile(rootdir,'qtrax','EH091202_15A01_p1');
    params.inputmoviefilestr = 'EH091202_15A01_p1.mp4';
    params.inputmoviefile = fullfile(params.inputdir,params.inputmoviefilestr);
    params.inputfiledict = {'featfile','EH091202_15A01_p1_1_feat.mat',fullfile(params.inputdir,'EH091202_15A01_p1_1_feat.mat')
      'roifile','EH091202_15A01_p1_1_roi.mat',fullfile(params.inputdir,'EH091202_15A01_p1_1_roi.mat')};
    [~,params.experiment_name] = fileparts(params.inputdir);
    params.experiment_name = [params.experiment_name,'__Qtrax'];
    params.outputdir = '../tests/PrepareJAABADataTests';
    [~,~,ext] = fileparts(params.inputmoviefilestr);
    params.outputmoviefilestr = ['movie',ext];
    params.outputtrxfilestr = 'trx.mat';
    params.outputperframedir = 'perframe';
    
    % sample jab file
    params.jabfile = fullfile(rootdir,'qtrax','qtrax_dummy20130819.jab');
    
    % options
    params.options = struct;
    params.options.SoftLinkFiles = true;
    params.options.fliplr = false;
    params.options.flipud = false;
    params.options.dotransposeimage = false;
    params.options.OverRideFPS = false;
    params.options.OverRideArena = false;
    params.options.CropFirstFrame = 1;
    params.options.CropEndFrame = inf;
    params.options.ArenaType = 'Circle';
    
    % read arena or fps?
    params.readarena = false;
    params.readfps = false;
    
  case 'MAGATAnalyzer',
    
    params.inputdir = fullfile(rootdir,'larvae_samuel','20120206_100300');
    params.inputmoviefilestr = '20120206_100300@FCF_attp2_1500062@UAS_TNT_2_0003@t7@n#n#n#n@30@.mmf';
    params.inputmoviefile = fullfile(params.inputdir,params.inputmoviefilestr);
    params.inputfiledict = {'expfile','expfile.mat',fullfile(params.inputdir,'FCF_attp2_1500062@UAS_TNT_2_0003_experiment_20120206_100300.mat')};
    % experiment name is encoded in movie name
    [~,params.experiment_name,ext] = fileparts(params.inputmoviefile);
    params.experiment_name = [params.experiment_name,'__MAGATAnalyzer'];
    params.outputdir = '../tests/PrepareJAABADataTests';
    params.outputmoviefilestr = ['movie',ext];
    params.outputtrxfilestr = 'trx.mat';
    params.outputperframedir = 'perframe';
    
    % sample jab file
    params.jabfile = fullfile(rootdir,'larvae_samuel','MAGATAnalyzer_Test.jab');
    
    % options
    params.options = struct;
    params.options.SoftLinkFiles = true;
    params.options.fliplr = false;
    params.options.flipud = false;
    params.options.dotransposeimage = false;
    params.options.OverRideFPS = false;
    params.options.OverRideArena = false;
    params.options.CropFirstFrame = 1;
    params.options.CropEndFrame = inf;
    params.options.ArenaType = 'None';
    
    % read arena or fps?
    params.readarena = false;
    params.readfps = false;  
    
  case 'MWT',
    
    params.inputdir = fullfile(rootdir,'larvae_mwt','MZZ_ppk1d9GAL4_t2_TrpA_20120612_110553');
    params.inputmoviefilestr = '';
    params.inputmoviefile = '';
    tmp = dir(fullfile(params.inputdir,'*.blob*'));
    blobspaths = cellfun(@(x) fullfile(params.inputdir,x),{tmp.name},'UniformOutput',false);    
    tmp = dir(fullfile(params.inputdir,'*.dat'));
    datpaths = cellfun(@(x) fullfile(params.inputdir,x),{tmp.name},'UniformOutput',false);
    params.inputfiledict = {'blobsfile','trx.blobs',blobspaths
      'datfiles','chor.dat',datpaths};
    [~,params.experiment_name] = fileparts(params.inputdir);
    params.experiment_name = [params.experiment_name,'__MWT'];
    params.outputdir = '../tests/PrepareJAABADataTests';
    params.outputmoviefilestr = params.inputmoviefilestr;
    params.outputtrxfilestr = 'trx.mat';
    params.outputperframedir = 'perframe';
    
    % sample jab file
    params.jabfile = fullfile(rootdir,'larvae_mwt','LarvaeMWT_Roll_v18TO.jab');
    
    % options
    params.options = struct;
    params.options.SoftLinkFiles = true;
    params.options.fliplr = false;
    params.options.flipud = false;
    params.options.dotransposeimage = false;
    params.options.fps = 30.319466231015564;
    params.options.pxpermm = 17.857142387323400;
    params.options.OverRideFPS = false;
    params.options.OverRideArena = false;
    params.options.CropFirstFrame = 1;
    params.options.CropEndFrame = inf;
    params.options.ArenaType = 'None';
    
    % read arena or fps?
    params.readarena = true;
    params.readfps = false;
    
  case 'LarvaeReid',
    
    name1 = 'FCF_attP2_1500062_50';
    name2 = 'o_ethylbutyrate015m8ul_0s1x360s0s#n#n#n@20';
    name3 = '20120427_134403';
    inputdir1 = fullfile(rootdir,'larvae_reid',name1,name2);
    params.inputdir = fullfile(inputdir1,name3);
    params.inputmoviefilestr = '';
    params.inputmoviefile = '';
    tmp = dir(fullfile(params.inputdir,'*.blob*'));
    blobspaths = cellfun(@(x) fullfile(params.inputdir,x),{tmp.name},'UniformOutput',false);    
    params.inputfiledict = {'blobsfile','trx.blobs',blobspaths
      'kinmatfile','kinVariables.mat',fullfile(inputdir1,'kinVariables.mat')
      'eventmatfile','eventVariables.mat',fullfile(inputdir1,'eventVariables.mat')};
    params.experiment_name = sprintf('%s__%s__%s__LarvaeReid',name1,name2,name3);
    params.outputdir = '../tests/PrepareJAABADataTests';
    params.outputmoviefilestr = params.inputmoviefilestr;
    params.outputtrxfilestr = 'trx.mat';
    params.outputperframedir = 'perframe';
    
    % sample jab file
    params.jabfile = fullfile(rootdir,'larvae_reid','LarvaeReid_dummy.jab');
    
    % options
    params.options = struct;
    params.options.SoftLinkFiles = true;
    params.options.fliplr = false;
    params.options.flipud = false;
    params.options.dotransposeimage = false;
    params.options.fps = 30.319466231015564;
    params.options.pxpermm = 17.857142387323400;
    params.options.OverRideFPS = false;
    params.options.OverRideArena = false;
    params.options.CropFirstFrame = 1;
    params.options.CropEndFrame = inf;
    params.options.ArenaType = 'None';
    
    % read arena or fps?
    params.readarena = true;
    params.readfps = false;
    
  case 'LarvaeLouis',
    
    name1 = 'rawdata';
    name2 = '20111212-120514';
    inputdir1 = fullfile(rootdir,'larvae_louis',name1);
    params.inputdir = fullfile(inputdir1,name2);
    params.inputmoviefilestr = '20111212-120514_video.avi';
    params.inputmoviefile = fullfile(params.inputdir,params.inputmoviefilestr);
    params.inputfiledict = {'indatafile','20111212-120514_data.txt',fullfile(params.inputdir,'20111212-120514_data.txt')
      'inkinmatfile','kinVariables.mat',fullfile(inputdir1,'kinVariables.mat')};
    params.experiment_name = sprintf('%s__LarvaeLouis',name2);
    params.outputdir = '../tests/PrepareJAABADataTests';
    [~,~,ext] = fileparts(params.inputmoviefilestr);
    params.outputmoviefilestr = ['movie',ext];
    params.outputtrxfilestr = 'trx.mat';
    params.outputperframedir = 'perframe';
    
    % sample jab file
    params.jabfile = fullfile(rootdir,'larvae_louis','LarvaeLouis_dummy.jab');
    
    % options
    params.options = struct;
    params.options.SoftLinkFiles = true;
    params.options.fliplr = false;
    params.options.flipud = false;
    params.options.dotransposeimage = false;
    params.options.OverRideFPS = false;
    params.options.OverRideArena = false;
    params.options.CropFirstFrame = 1;
    params.options.CropEndFrame = inf;
    params.options.ArenaType = 'None';
    
    % read arena or fps?
    params.readarena = false;
    params.readfps = false;
    
    
end

%% set PrepareJAABAData fields to match

%% set data type
handles = guidata(hfig);
handles.InputDataTypeIndex = find(strcmp(params.datatype,handles.InputDataTypeNames));
handles.InputDataType = params.datatype;
handles = PrepareJAABAData('UpdateGUI',handles);
guidata(hfig,handles);

%% set input table entries
handles = guidata(hfig);

% movie
handles.inputmoviefilestr = params.inputmoviefilestr;
handles.InputVideoFile = params.inputmoviefile;

% rest of files
InputDataType = handles.InputDataTypes.(handles.InputDataType);
for i = 1:size(params.inputfiledict,1),
  j = find(strcmp({InputDataType.files.code},params.inputfiledict{i,1}));
  if numel(j) ~= 1,
    error('Could not match input file type %s',params.inputfiledict{i,1});
  end
  handles.inputfilestrs.(handles.InputDataType){j} = params.inputfiledict{i,2};
  handles.InputFiles.(handles.InputDataType){j} = params.inputfiledict{i,3};
end

handles = PrepareJAABAData('UpdateGUI',handles);
guidata(hfig,handles);

%% set output table entries

handles = guidata(hfig);

% experiment directory
handles.ExperimentDirectory = fullfile(params.outputdir,params.experiment_name);

% output file strings
handles.moviefilestr = params.outputmoviefilestr;
handles.trxfilestr = params.outputtrxfilestr;
handles.perframedirstr = params.outputperframedir;

handles = PrepareJAABAData('UpdateGUI',handles);
guidata(hfig,handles);

%% options

handles = guidata(hfig);

fns = fieldnames(params.options);
for i = 1:numel(fns),
  handles.(fns{i}) = params.options.(fns{i});
end

handles = PrepareJAABAData('UpdateGUI',handles);
guidata(hfig,handles);

%% read arena parameters

if params.readarena,

  PrepareJAABAData('pushbutton_ReadArenaParameters_Callback',hfig,[],guidata(hfig));
  handles = guidata(hfig);

end
  
%% read fps parameters

if params.readfps,

  PrepareJAABAData('pushbutton_ReadFPS_Callback',hfig,[],guidata(hfig));
  handles = guidata(hfig);

end

%% convert

outexpdir = fullfile(params.outputdir,params.experiment_name);
if exist(outexpdir,'dir'),
  rmdir(outexpdir,'s');
end

PrepareJAABAData('pushbutton_Convert_Callback',hfig,[],guidata(hfig));
handles = guidata(hfig);

%% check that this worked right

filesneeded = {
  fullfile(outexpdir,params.outputtrxfilestr)
  fullfile(outexpdir,params.outputperframedir)
  };
if ~isempty(params.outputmoviefilestr),
  filesneeded{end+1} = fullfile(outexpdir,params.outputmoviefilestr);
end

ismissing = false(size(filesneeded));
for i = 1:numel(filesneeded),
  ismissing(i) = ~exist(filesneeded{i},'file');
  if ismissing(i) && ispc,
    [~,didfind] = GetPCShortcutFileActualPath(filesneeded{i});
    ismissing(i) = didfind == 0;
  end
  fprintf('File %s exists? %d\n',filesneeded{i},~ismissing(i));
end
if any(ismissing),
  error(['Files missing:',sprintf(' %s',filesneeded{ismissing})]);
end

%% close the gui

PrepareJAABAData('figure1_CloseRequestFcn',hfig,[],guidata(hfig));

%% see if we can run JAABADetect on this experiment

JAABADetect(outexpdir,'jabfiles',{params.jabfile},'forcecompute',true);

success = true;