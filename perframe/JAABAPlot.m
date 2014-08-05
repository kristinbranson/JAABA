function varargout = JAABAPlot(varargin)
%JAABAPLOT: Plot the output of JAABA
%
% This program is part of JAABA.
%
% JAABA: The Janelia Automatic Animal Behavior Annotator
% Copyright 2012, Kristin Branson, HHMI Janelia Farm Resarch Campus
% http://jaaba.sourceforge.net/
% bransonk@janelia.hhmi.org
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License (version 3 pasted in LICENSE.txt) for 
% more details.

% Last Modified by GUIDE v2.5 11-Mar-2013 15:51:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @JAABAPlot_OpeningFcn, ...
                   'gui_OutputFcn',  @JAABAPlot_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% ---
function handles=initialize(handles)

handles.grouplist={};
handles.groupvalue=1;
handles.experimentlist={{}};
handles.experimentvalue={1};
handles.classifierlist={};
handles.classifiervalue=1;
%handles.configurations={};
handles.scorefiles={};
handles.analysis='';
handles.behaviorlist={};
handles.behaviornot=0;
handles.behaviorvalue=1;
handles.behaviorlogic=1;
handles.behaviorvalue2=1;
handles.behaviornormalizenot=0;
handles.behaviorvalue3=1;
handles.features={};
handles.featurelist={};
handles.featurevalue=1;
handles.individuals_behavior=[];
handles.individuals_feature=[];
handles.individuallist={'All'};
handles.individualvalue=1;
handles.individualidx='A';
handles.sexdata={};
handles.fps=nan;
handles.classify_forcecompute=false;
handles.behaviorbarchart_style=1;
handles.behaviortimeseries_style=1;
handles.featurehistogram_style=1;
handles.featurehistogram_style2=1;
handles.comparison=0;
handles.logbinsize=0;
handles.nbins=100;
handles.featuretimeseries_style=1;
handles.featuretimeseries_style2=1;
handles.subtractmean=0;
handles.windowradius=10;
handles.boutstats_style=1;
handles.boutstats_style2=1;
handles.omitnan=1;
handles.omitinf=1;
handles.absdprimezscore=1;
handles.comparison2=0;
handles.dump2csv=1;
handles.dump2mat=1;
handles.centraltendency=1;
handles.dispersion=1;
handles.xoffset=1;
handles.minimumtrajectorylength=1;
handles.convolutionwidth=10;
handles.pvalue=0.01;
handles.interestingfeaturehistograms_cache=[];
handles.interestingfeaturetimeseries_cache=[];
handles.defaultcolors=[1 0 0;  0 0.5 0;  0 0 1;  0 1 1;  1 0 1;  0.749 0.749 0;  0 0 0];
handles.colors=[];
handles.trx_file='registered_trx.mat';
handles.perframe_dir='perframe';

% fnssave = ConfigFileFeatures2Save()
% returns the field names that should get saved and loaded into the rc/configuration
% files

function fnssave = ConfigFileFeatures2Save()

fnssave = {
  'grouplist'
  'groupvalue'
  'experimentlist'
  'experimentvalue'
  'classifierlist'
  'classifiervalue'
  'analysis'
  'scorefiles'
  'behaviorlist'
  'behaviornot'
  'behaviorvalue'
  'behaviorlogic'
  'behaviorvalue2'
  'behaviornormalizenot'
  'behaviorvalue3'
  'features'
  'featurelist'
  'featurevalue'
  'individuals_behavior'
  'individuals_feature'
  'individuallist'
  'individualvalue'
  'individualidx'
  'sexdata'
  'fps'
  'classify_forcecompute'
  'behaviorbarchart_style'
  'behaviortimeseries_style'
  'featurehistogram_style'
  'featurehistogram_style2'
  'comparison'
  'logbinsize'
  'nbins'
  'featuretimeseries_style'
  'featuretimeseries_style2'
  'subtractmean'
  'windowradius'
  'boutstats_style'
  'boutstats_style2'
  'omitnan'
  'omitinf'
  'absdprimezscore'
  'comparison2'
  'dump2csv'
  'dump2mat'
  'centraltendency'
  'dispersion'
  'xoffset'
  'minimumtrajectorylength'
  'convolutionwidth'
  'pvalue'
  'interestingfeaturehistograms_cache'
  'interestingfeaturetimeseries_cache'
  'defaultcolors'
  'colors'
  'trx_file'
  'perframe_dir'
  };
  

% ---
function handles=load_configuration_file(filename,hObject,eventdata,handles)

% try to load from filename
if ~exist(filename,'file'),
  uiwait(warndlg('File %s does not exist.',filename));
  return;
end

handles_saved=load(filename);
handles_saved=handles_saved.handles;

fnssave = ConfigFileFeatures2Save();

iserror = false;
for i = 1:numel(fnssave),
  
  fn = fnssave{i};
  
  if isfield(handles_saved,fn),
    try
      handles.(fn) = handles_saved.(fn);
    catch ME,
      warning('Could not set handles.%s to handles_saved.%s: %s',fn,fn,getReport(ME));
      iserror = true;
      break;
    end
  else
    warning('Required field %s not saved in configuration file %s.',fn,filename);
    iserror = true;
  end
  
end

if iserror,
  
  handles=initialize(handles);
  handles.grouplist=handles_saved.grouplist;
  handles.groupvalue=handles_saved.groupvalue;
  handles.experimentlist=handles_saved.experimentlist;
  handles.experimentvalue=handles_saved.experimentvalue;
  handles.colors=zeros(length(handles.grouplist),3);
  handles=update_experiment_data(handles,true,true,true);
  uiwait(warndlg([filename ' is in an old format.  only group and experiment lists are salvageable.  save a new version']));  drawnow;
  
end

% handles.grouplist=handles_saved.grouplist;
% handles.groupvalue=handles_saved.groupvalue;
% handles.experimentlist=handles_saved.experimentlist;
% handles.experimentvalue=handles_saved.experimentvalue;
% try
%   handles.classifierlist=handles_saved.classifierlist;
%   handles.classifiervalue=handles_saved.classifiervalue;
% %  handles.configurations=handles_saved.configurations;
%   handles.analysis=handles_saved.analysis;
%   handles.scorefiles=handles_saved.scorefiles;
%   handles.behaviorlist=handles_saved.behaviorlist;
%   handles.behaviornot=handles_saved.behaviornot;
%   handles.behaviorvalue=handles_saved.behaviorvalue;
%   handles.behaviorlogic=handles_saved.behaviorlogic;
%   handles.behaviorvalue2=handles_saved.behaviorvalue2;
%   handles.behaviornormalizenot=handles_saved.behaviornormalizenot;
%   handles.behaviorvalue3=handles_saved.behaviorvalue3;
%   handles.features=handles_saved.features;
%   handles.featurelist=handles_saved.featurelist;
%   handles.featurevalue=handles_saved.featurevalue;
%   handles.individuals_behavior=handles_saved.individuals_behavior;
%   handles.individuals_feature=handles_saved.individuals_feature;
%   handles.individuallist=handles_saved.individuallist;
%   handles.individualvalue=handles_saved.individualvalue;
%   handles.individualidx=handles_saved.individualidx;
%   handles.sexdata=handles_saved.sexdata;
%   handles.fps=handles_saved.fps;
%   handles.classify_forcecompute=handles_saved.classify_forcecompute;
%   handles.behaviorbarchart_style=handles_saved.behaviorbarchart_style;
%   handles.behaviortimeseries_style=handles_saved.behaviortimeseries_style;
%   handles.featurehistogram_style=handles_saved.featurehistogram_style;
%   handles.featurehistogram_style2=handles_saved.featurehistogram_style2;
%   handles.comparison=handles_saved.comparison;
%   handles.logbinsize=handles_saved.logbinsize;
%   handles.nbins=handles_saved.nbins;
%   handles.featuretimeseries_style=handles_saved.featuretimeseries_style;
%   handles.featuretimeseries_style2=handles_saved.featuretimeseries_style2;
%   handles.subtractmean=handles_saved.subtractmean;
%   handles.windowradius=handles_saved.windowradius;
%   handles.boutstats_style=handles_saved.boutstats_style;
%   handles.boutstats_style2=handles_saved.boutstats_style2;
%   handles.omitnan=handles_saved.omitnan;
%   handles.omitinf=handles_saved.omitinf;
%   handles.absdprimezscore=handles_saved.absdprimezscore;
%   handles.comparison2=handles_saved.comparison2;
%   handles.dump2csv=handles_saved.dump2csv;
%   handles.dump2mat=handles_saved.dump2mat;
%   handles.centraltendency=handles_saved.centraltendency;
%   handles.dispersion=handles_saved.dispersion;
%   handles.xoffset=handles_saved.xoffset;
%   handles.minimumtrajectorylength=handles_saved.minimumtrajectorylength;
%   handles.convolutionwidth=handles_saved.convolutionwidth;
%   handles.pvalue=handles_saved.pvalue;
%   handles.interestingfeaturehistograms_cache=handles_saved.interestingfeaturehistograms_cache;
%   handles.interestingfeaturetimeseries_cache=handles_saved.interestingfeaturetimeseries_cache;
%   handles.defaultcolors=handles_saved.defaultcolors;
%   handles.colors=handles_saved.colors;
%   handles.trx_file=handles_saved.trx_file;
%   handles.perframe_dir=handles_saved.perframe_dir;
%   handles.table=[];
% catch me
%   handles=initialize(handles);
%   handles.grouplist=handles_saved.grouplist;
%   handles.groupvalue=handles_saved.groupvalue;
%   handles.experimentlist=handles_saved.experimentlist;
%   handles.experimentvalue=handles_saved.experimentvalue;
%   handles.colors=zeros(length(handles.grouplist),3);
%   handles=update_experiment_data(handles,true,true,true);
%   uiwait(warndlg([filename ' is in an old format.  only group and experiment lists are salvageable.  save a new version']));  drawnow;
% end


% ---
function behavior_data=update_t01s_from_postprocessed(behavior_data)

for i=1:length(behavior_data.allScores.postprocessed)
  [behavior_data.allScores.t0s{i} behavior_data.allScores.t1s{i}]=...
      get_interval_ends(behavior_data.allScores.postprocessed{i}==1);
end


% ---
function ret_val=get_nindividuals_behavior(experiment_path,score_file,classifier_timestamp)

try
  behavior_data=load(fullfile(experiment_path,score_file));
catch
  behavior_data.allScores.t0s=[];
  behavior_data.timestamp=nan;
end
if ((classifier_timestamp ~= behavior_data.timestamp) || (length(behavior_data.allScores.t0s)==0))
  ret_val=-1;
else
  ret_val=length(behavior_data.allScores.t0s);
end


% ---
function ret_val=get_nindividuals_feature(experiment_path,perframe_dir,feature_file)

% could do more error checking here by loading them all in

try
  feature_data=load(fullfile(experiment_path,perframe_dir,[feature_file{1} '.mat']));
catch
  feature_data.data={};
end
if (length(feature_data.data)==0)
  ret_val=-1;
else
  ret_val=length(feature_data.data);
end


% ---
function [fps,trxfile]=get_fps(expdir,trxfile)

filename=fullfile(expdir,trxfile);
if(~exist(filename) || isempty(trxfile))
  [trxfile,~,~]=uigetfile(fullfile(expdir,'*.mat'),...
      sprintf('Select trx file for %s',expdir));
  if(isnumeric(trxfile) && (trxfile==0))
    fps=nan;
    trxfile='';
    return;
  end
  filename=fullfile(expdir,trxfile);
  %[fps,trxfile]=get_fps(expdir,trxfile);
  %return;
end

t=load(filename);
if isfield(t.trx(1),'fps'),
  fps=t.trx(1).fps;
elseif isfield(t.trx(1),'dt'),
  fps = 1/mean(t.trx(1).dt(1:10));
  if isnan(fps)
    uiwait(warndlg('Trx file does not have recording fps (frames per second). Assuming fps as 30'));
    drawnow;
    fps = 30;
  end
else
  uiwait(warndlg('Trx file does not have recording fps (frames per second). Assuming fps as 30'));
  drawnow;
  fps = 30;
end


% ---
function handles=update_experiment_data(handles,features,sexdata,individuals)

handlesexperimentlist=[handles.experimentlist{:}];

handlesfeatures=cell(1,length(handlesexperimentlist));
handlessexdata=cell(1,length(handlesexperimentlist));
handlesindividualsbehavior=zeros(length(handlesexperimentlist),length(handles.scorefiles));
handlesindividualsfeature=zeros(1,length(handlesexperimentlist));

parfor ge=1:length(handlesexperimentlist)
  if(features)
    tmp=dir(fullfile(handlesexperimentlist{ge},handles.perframe_dir,'*.mat'));
    [handlesfeatures{ge}{1:length(tmp)}]=deal(tmp.name);
    handlesfeatures{ge}=cellfun(@(x) x(1:(end-4)),handlesfeatures{ge},'uniformoutput',false);
  end

  if(sexdata)
    fullfile(handlesexperimentlist{ge},handles.perframe_dir,'sex.mat');
    if(exist(ans,'file'))
      tmp=load(ans);
      cellfun(@(x) strcmp(x,'M'),tmp.data,'uniformoutput',false);
      handlessexdata(ge)={ans};
    else
      tmp=dir(fullfile(handlesexperimentlist{ge},handles.perframe_dir,'*.mat'));
      tmp=load(fullfile(handlesexperimentlist{ge},handles.perframe_dir,tmp(1).name));
      cellfun(@(x) nan(1,length(x)),tmp.data,'uniformoutput',false);
      handlessexdata(ge)={ans};
    end
  end

  if(individuals)
    behavior_data=[];
    parfor_tmp=zeros(1,length(handles.scorefiles));
    for s=1:length(handles.scorefiles)
      %classifier=load(handles.classifierlist{s});
      classifier=load(handles.classifierlist{s},'-mat');
      %classifier=x;
      parfor_tmp(s)=get_nindividuals_behavior(handlesexperimentlist{ge},handles.scorefiles{s},...
          classifier.x.classifierStuff.timeStamp);
    end
    handlesindividualsbehavior(ge,:)=parfor_tmp;

    handlesindividualsfeature(ge)=...
        get_nindividuals_feature(handlesexperimentlist{ge},handles.perframe_dir,handlesfeatures{ge});
  end
end

if(features)
  handles.features={handlesfeatures{:}};
  handles.featurelist=check_for_diff_and_return_intersection(handles.features);
end

if(sexdata)
  handles.sexdata={handlessexdata{:}};
end

if(individuals)
  handles.individuals_behavior=handlesindividualsbehavior;
  handles.individuals_feature=handlesindividualsfeature;
  handles=fillin_individuallist(handles);
end

%classifier=load(handles.classifierlist{1});
%handles.fps=get_fps(fullfile(handlesexperiments{1},classifier.trxfilename));
[handles.fps,handles.trx_file]=get_fps(handlesexperimentlist{1},handles.trx_file);

handles.interestingfeaturehistograms_cache=[];
handles.interestingfeaturetimeseries_cache=[];


% ---
function update_figure(handles)

set(handles.DumpToCSV,'enable','off');
set(handles.DumpToMAT,'enable','off');

if(isempty(handles.grouplist) || length(handles.experimentlist{handles.groupvalue})==0)
  set(handles.ExperimentList,'enable','off');
else
  set(handles.ExperimentList,'enable','on');
end
if(length(handles.grouplist)<2)
  set(handles.ExperimentMove,'enable','off');
else
  set(handles.ExperimentMove,'enable','on');
end
if(isempty(handles.grouplist))
  set(handles.GroupList,'enable','off');
  set(handles.GroupChange,'enable','off');
  set(handles.ExperimentAdd,'enable','off');
  set(handles.ExperimentDelete,'enable','off');
else
  set(handles.GroupList,'enable','on');
  set(handles.GroupChange,'enable','on');
  set(handles.ExperimentAdd,'enable','on');
  set(handles.ExperimentDelete,'enable','on');
end
if(sum(cellfun(@length,handles.experimentlist))==0)
  set(handles.ClassifierAuto,'enable','off');
  set(handles.ClassifierCheck,'enable','off');
else
  set(handles.ClassifierAuto,'enable','on');
  set(handles.ClassifierCheck,'enable','on');
  set(handles.DumpToCSV,'enable','on');
  set(handles.DumpToMAT,'enable','on');
end
if(isempty(handles.classifierlist))
  set(handles.ClassifierList,'enable','off');
  set(handles.ClassifierDelete,'enable','off');
else
  set(handles.ClassifierList,'enable','on');
  set(handles.ClassifierDelete,'enable','on');
end
if((sum(cellfun(@length,handles.experimentlist))==0) || (isempty(handles.classifierlist)))
  set(handles.ClassifierClassify,'enable','off');
else
  set(handles.ClassifierClassify,'enable','on');
end

if(sum(cellfun(@length,handles.experimentlist))==0)
  set(handles.FeatureHistogram,'enable','off');
  set(handles.FeatureTimeSeries,'enable','off');
  set(handles.InterestingFeatureHistograms,'enable','off');
else
  set(handles.FeatureHistogram,'enable','on');
  set(handles.FeatureTimeSeries,'enable','on');
  set(handles.InterestingFeatureHistograms,'enable','on');
end
if((sum(cellfun(@length,handles.experimentlist))==0) || ...
   (isempty(handles.classifierlist)))
%   (sum(sum(handles.individuals_behavior==-1))>0) || ...
%   (sum(sum(diff(handles.individuals_behavior,[],2)~=0))>0))
  set(handles.BehaviorBarChart,'enable','off');
  set(handles.BehaviorTimeSeries,'enable','off');
  set(handles.BoutStats,'enable','off');
else
  set(handles.BehaviorBarChart,'enable','on');
  set(handles.BehaviorTimeSeries,'enable','on');
  set(handles.BoutStats,'enable','on');
end

set(handles.StyleList,'enable','off');              set(handles.StyleList2,'enable','off');
set(handles.BehaviorNot,'enable','off');            set(handles.BehaviorList,'enable','off');
set(handles.BehaviorLogic,'enable','off');          set(handles.BehaviorList3,'enable','off');
set(handles.BehaviorNormalizeNot,'enable','off');   set(handles.BehaviorList2,'enable','off');
set(handles.FeatureList,'enable','off');            set(handles.IndividualList,'enable','off');
set(handles.ConvolutionWidth,'enable','off');       set(handles.PValue,'enable','off');
set(handles.NBins,'enable','off');                  set(handles.WindowRadius,'enable','off');
set(handles.XOffset,'enable','off');                set(handles.MinimumTrajectoryLength,'enable','off');
set(handles.OmitInf,'enable','off');                set(handles.OmitNaN,'enable','off');
set(handles.AbsDPrimeZScore,'enable','off');        set(handles.SubtractMean,'enable','off');
set(handles.LogBinSize,'enable','off');             set(handles.AllFrames,'enable','off');
set(handles.NotDuring,'enable','off');
set(handles.CentralTendency,'enable','off');        set(handles.Dispersion,'enable','off');
set(handles.DPrime,'enable','off');                 set(handles.ZScore,'enable','off');
set(handles.Plot,'enable','off');

bg=get(handles.GroupNew,'backgroundcolor');
fg=get(handles.GroupNew,'foregroundcolor');
set(handles.FeatureHistogram,'backgroundcolor',bg);
set(handles.FeatureHistogram,'foregroundcolor',fg);
set(handles.FeatureTimeSeries,'backgroundcolor',bg);
set(handles.FeatureTimeSeries,'foregroundcolor',fg);
set(handles.BehaviorBarChart,'backgroundcolor',bg);
set(handles.BehaviorBarChart,'foregroundcolor',fg);
set(handles.BehaviorTimeSeries,'backgroundcolor',bg);
set(handles.BehaviorTimeSeries,'foregroundcolor',fg);
set(handles.BoutStats,'backgroundcolor',bg);
set(handles.BoutStats,'foregroundcolor',fg);
set(handles.InterestingFeatureHistograms,'backgroundcolor',bg);
set(handles.InterestingFeatureHistograms,'foregroundcolor',fg);
analysis2=handles.analysis;

if(isempty(handles.grouplist))
  set(handles.GroupList,'String',{''},'Value',1);
else
  tmp=length(handles.grouplist);
  cellstr(strcat(repmat('<html><font color=#"',tmp,1),...
      reshape(dec2hex(round(handles.colors'*255),2)',6,size(handles.colors,1))',...
      repmat('">',tmp,1),handles.grouplist',repmat('</font></html>',tmp,1)));
  set(handles.GroupList,'String',ans,'Value',handles.groupvalue);
end
if(isempty(handles.experimentlist))
  set(handles.ExperimentList,'String',{''},'Value',1);
else
  set(handles.ExperimentList,'String',handles.experimentlist{handles.groupvalue});
  set(handles.ExperimentList,'Value',handles.experimentvalue{handles.groupvalue});
end
if(isempty(handles.classifierlist))
  set(handles.ClassifierList,'String',{''},'Value',1);
else
  set(handles.ClassifierList,'String',handles.classifierlist);
  set(handles.ClassifierList,'Value',handles.classifiervalue);
end
if(isempty(handles.behaviorlist))
  set(handles.BehaviorList,'String',{''},'Value',1);
  set(handles.BehaviorList2,'String',{''},'Value',1);
  set(handles.BehaviorList3,'String',{''},'Value',1);
else
  set(handles.BehaviorList,'String',{handles.behaviorlist{:} 'All behaviors'},'Value',handles.behaviorvalue);
  set(handles.BehaviorList2,'String',handles.behaviorlist,'Value',handles.behaviorvalue2);
  set(handles.BehaviorList3,'String',{'All frames' handles.behaviorlist{:}},'Value',handles.behaviorvalue3);
end
set(handles.BehaviorNot,'Value',handles.behaviornot);
set(handles.BehaviorLogic,'Value',handles.behaviorlogic);
set(handles.BehaviorNormalizeNot,'Value',handles.behaviornormalizenot);
if(isempty(handles.featurelist))
  set(handles.FeatureList,'String',{''},'Value',1);
else
  set(handles.FeatureList,'String',handles.featurelist,'Value',handles.featurevalue);
end
%set(handles.Table,'Data',[]);

switch(handles.analysis)
  case 'feature_histogram'
    set(handles.FeatureHistogram,'backgroundcolor',fg);
    set(handles.FeatureHistogram,'foregroundcolor',bg);
    set(handles.StyleList,'string',handles.featurehistogram_stylelist,...
        'value',handles.featurehistogram_style);
    set(handles.StyleList2,'string',handles.featurehistogram_stylelist2,...
        'value',handles.featurehistogram_style2);
    if(~isempty(handles.behaviorlist))
      get(handles.BehaviorList,'String');
      set(handles.BehaviorList,'String',{ans{:} 'All frames'});
    end
    if(strcmp(get(handles.FeatureHistogram,'enable'),'off'))
      analysis2='';
    else
      set(handles.FeatureList,'enable','on');
      set(handles.PValue,'enable','on');
      set(handles.CentralTendency,'enable','on');
      set(handles.Dispersion,'enable','on');
      set(handles.LogBinSize,'enable','on');
      if(~isempty(handles.classifierlist))
        set(handles.AllFrames,'enable','on');
        set(handles.NotDuring,'enable','on');
      end
      set(handles.NBins,'enable','on');
    end
  case 'feature_timeseries'
    set(handles.FeatureTimeSeries,'backgroundcolor',fg);
    set(handles.FeatureTimeSeries,'foregroundcolor',bg);
    set(handles.StyleList,'string',handles.featuretimeseries_stylelist,...
        'value',handles.featuretimeseries_style);
    if(length(handles.classifierlist)>0)
      set(handles.StyleList2,'string',handles.featuretimeseries_stylelist2,...
          'value',handles.featuretimeseries_style2);
    else
      set(handles.StyleList2,'string',{''},'value',1);
    end
    if(strcmp(get(handles.FeatureTimeSeries,'enable'),'off'))
      analysis2='';
    else
      set(handles.FeatureList,'enable','on');
      set(handles.ConvolutionWidth,'enable','on');
      set(handles.CentralTendency,'enable','on');
      set(handles.Dispersion,'enable','on');
      if(handles.featuretimeseries_style2==1)
        set(handles.XOffset,'enable','on');
      else
        set(handles.SubtractMean,'enable','on');
        set(handles.WindowRadius,'enable','on');
      end
    end
  case 'behavior_barchart'
    set(handles.BehaviorBarChart,'backgroundcolor',fg);
    set(handles.BehaviorBarChart,'foregroundcolor',bg);
    set(handles.StyleList,'string',handles.behaviorbarchart_stylelist,...
        'value',handles.behaviorbarchart_style);
    set(handles.StyleList2,'string',{''},'value',1);
    if(strcmp(get(handles.BehaviorBarChart,'enable'),'off'))
      analysis2='';
    else
      set(handles.PValue,'enable','on');
      set(handles.CentralTendency,'enable','on');
      set(handles.Dispersion,'enable','on');
    end
  case 'behavior_timeseries'
    set(handles.BehaviorTimeSeries,'backgroundcolor',fg);
    set(handles.BehaviorTimeSeries,'foregroundcolor',bg);
    set(handles.StyleList,'string',handles.behaviortimeseries_stylelist,...
        'value',handles.behaviortimeseries_style);
    set(handles.StyleList2,'string',{''},'value',1);
    if(strcmp(get(handles.BehaviorTimeSeries,'enable'),'off'))
      analysis2='';
    else
      set(handles.ConvolutionWidth,'enable','on');
      set(handles.XOffset,'enable','on');
      set(handles.CentralTendency,'enable','on');
      set(handles.Dispersion,'enable','on');
    end
  case 'bout_stats'
    set(handles.BoutStats,'backgroundcolor',fg);
    set(handles.BoutStats,'foregroundcolor',bg);
    set(handles.StyleList,'string',handles.boutstats_stylelist,...
        'value',handles.boutstats_style);
    set(handles.StyleList2,'string',handles.boutstats_stylelist2,...
        'value',handles.boutstats_style2);
    if(strcmp(get(handles.BoutStats,'enable'),'off'))
      analysis2='';
    else
      set(handles.PValue,'enable','on');
      set(handles.CentralTendency,'enable','on');
      set(handles.Dispersion,'enable','on');
    end
  case 'interesting_feature_histograms'
    set(handles.InterestingFeatureHistograms,'backgroundcolor',fg);
    set(handles.InterestingFeatureHistograms,'foregroundcolor',bg);
    set(handles.StyleList,'string',{''},'value',1);
%     set(handles.StyleList2,'string',{''},'value',1);
    set(handles.StyleList2,'string',handles.featurehistogram_stylelist2,...
        'value',handles.featurehistogram_style2);
    if(strcmp(get(handles.InterestingFeatureHistograms,'enable'),'off'))
      analysis2='';
    else
      set(handles.OmitInf,'enable','on');
      set(handles.OmitNaN,'enable','on');
      set(handles.AbsDPrimeZScore,'enable','on');
      set(handles.DPrime,'enable','on');
      set(handles.ZScore,'enable','on');
    end
  otherwise
    set(handles.StyleList,'string',{''},'value',1);
    set(handles.StyleList2,'string',{''},'value',1);
end

if(~isempty(analysis2))
  set(handles.IndividualList,'enable','on');
  set(handles.Plot,'enable','on');
  set(handles.MinimumTrajectoryLength,'enable','on');
  set(handles.DumpToCSV,'enable','on');
  set(handles.DumpToMAT,'enable','on');
  if((ismember(analysis2,{'feature_histogram','behavior_barchart','behavior_timeseries','bout_stats'}) || ...
     (strcmp(analysis2,'feature_timeseries')&&(handles.featuretimeseries_style2~=1))) && ...
      (length(handles.behaviorlist)>0))
    set(handles.BehaviorNot,'enable','on');
    set(handles.BehaviorList,'enable','on');
    if(length(handles.behaviorlist)>1)
      set(handles.BehaviorLogic,'enable','on');
    end
    if(handles.behaviorlogic>1)
      set(handles.BehaviorList2,'enable','on');
    end
  end
  if((ismember(analysis2,{'behavior_barchart','behavior_timeseries'})) && ...
      (length(handles.behaviorlist)>1))
    set(handles.BehaviorList3,'enable','on');
  end
  if(~isempty(handles.classifierlist) && (handles.behaviorvalue3>1))
    set(handles.BehaviorNormalizeNot,'enable','on');
  end
  if(~strcmp(analysis2,'interesting_feature_histograms'))
    set(handles.StyleList,'enable','on');
    if(isempty(handles.individuallist))
      set(handles.IndividualList,'String',{''},'Value',1);
    else
      set(handles.IndividualList,'String',handles.individuallist,'Value',handles.individualvalue);
    end
  else
    set(handles.IndividualList,'String',handles.individuallist(1:3),'Value',handles.individualvalue);
  end
  if(length(get(handles.StyleList2,'string'))>1)
    set(handles.StyleList2,'enable','on');
  end
end

menu_classify_forcecompute_set(handles);
button_comparison_set(handles);
button_comparison2_set(handles);

set(handles.MinimumTrajectoryLength,'string',handles.minimumtrajectorylength);
set(handles.ConvolutionWidth,'string',handles.convolutionwidth);
set(handles.PValue,'string',handles.pvalue);
set(handles.NBins,'string',handles.nbins);
set(handles.WindowRadius,'string',handles.windowradius);

set(handles.LogBinSize,'value',handles.logbinsize);
set(handles.SubtractMean,'value',handles.subtractmean);
set(handles.AbsDPrimeZScore,'value',handles.absdprimezscore);
set(handles.DumpToCSV,'value',handles.dump2csv);
set(handles.DumpToMAT,'value',handles.dump2mat);
set(handles.OmitNaN,'value',handles.omitnan);
set(handles.OmitInf,'value',handles.omitinf);
set(handles.XOffset,'value',handles.xoffset);
set(handles.CentralTendency,'value',handles.centraltendency);
set(handles.Dispersion,'value',handles.dispersion);


% --- Executes just before JAABAPlot is made visible.
function JAABAPlot_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

if ~isdeployed,
  SetUpJAABAPath;
end


handles.computation_threads = SetUpMatlabPoolforJAABAPlot;


handles.featurehistogram_stylelist=...
    {'Central Tendency','Central Tendency & Dispersion','Overlayed per-Exp Means', 'Box Plot'};
handles.featurehistogram_stylelist2=...
    {'per Frame','Mean per Bout','Median per Bout','Max per Bout','Min per Bout','Std. Dev. per Bout'};
handles.featuretimeseries_stylelist=...
    {'Central Tendency','Central Tendency & Dispersion','Overlayed per-Exp Means'};
handles.featuretimeseries_stylelist2=...
    {'Entire Recording','Onset Triggered','Offset Triggered'};
handles.behaviorbarchart_stylelist=...
    {'per Group','per Experiment, CT & D','per Experiment, box',...
     'per Fly, grouped','per Fly, scatter','per Fly, trajectory length'};
handles.behaviortimeseries_stylelist=...
    {'Central Tendency','Central Tendency & Dispersion','Overlayed per-Exp Means'};
handles.boutstats_stylelist=...
    {'per Experiment, CT & D', 'per Fly, grouped'};
handles.boutstats_stylelist2=...
    {'Bout Length','Inter-bout Length'};

if isdeployed,
  handles.rcfilename = deployedRelative2Global('most_recent_config.mat');
else
  handles.rcfilename = 'most_recent_config.mat';
end

if(exist(handles.rcfilename)==2)
  try
    handles=load_configuration_file(handles.rcfilename,hObject,eventdata,handles);
  catch ME
    %errordlg({'Error loading last used configuration. Using default values.',getReport(ME)},'Error loading last used configuration');
    handles=initialize(handles);
  end
else
  handles=initialize(handles);
end
update_figure(handles);

addedsystemmonitor = false;
try
  if ~isdeployed
    javaaddpath(fullfile(baseDir,'misc','javasysmon-0.3.4.jar'));
    addedsystemmonitor = true;
  else
    javaaddpath('javasysmon-0.3.4.jar');
    addedsystemmonitor = true;
  end
  catch ME,
    fprintf('Could not add javasysmon to path, disabling system monitor:\n%s\n',getReport(ME));
end

if addedsystemmonitor,
  import com.jezhumble.javasysmon.JavaSysMon.*
  handles.system_monitor.object=com.jezhumble.javasysmon.JavaSysMon();
  handles.system_monitor.timer=timer('Name','system_monitor','Period',1,'ExecutionMode','fixedRate',...
    'TimerFcn',@(hObject,eventdata)system_monitor_callback(hObject,eventdata,handles));
  warning('off','MATLAB:Java:ConvertFromOpaque');
  start(handles.system_monitor.timer);
else
  handles.system_monitor = [];
end


% Choose default command line output for JAABAPlot
handles.output = hObject;

set(hObject,'CloseRequestFcn',@figure_CloseRequestFcn);
set(handles.ExperimentList,'Callback',@ExperimentListCallback);
set(handles.ClassifierList,'Callback',@ClassifierListCallback);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes JAABAPlot wait for user response (see UIRESUME)
% uiwait(handles.figure1);


function system_monitor_callback(obj,src,handles)

persistent last_cpu

next_cpu=handles.system_monitor.object.cpuTimes();
mem=handles.system_monitor.object.physical();
if(~isempty(last_cpu))
  %disp(['cpu: ' num2str(last_cpu.getCpuUsage(next_cpu)) ', mem=' num2str(mem.getFreeBytes()/mem.getTotalBytes())]);
  set(handles.SystemMonitor,'string',[num2str(round(100*last_cpu.getCpuUsage(next_cpu))) '% cpu, ' ...
      num2str(round(100*(mem.getTotalBytes()-mem.getFreeBytes())/mem.getTotalBytes())) '% mem']);
  drawnow
end
last_cpu=next_cpu;

function [success,msg] = SaveConfiguration(handles,filename)

fnssave = ConfigFileFeatures2Save();
handles_save = struct;

for i = 1:numel(fnssave),
  
  fn = fnssave{i};
  
  if isfield(handles,fn),
    try
      handles_save.(fn) = handles.(fn);
    catch ME,
      warning('Could not set handles_save.%s to handles.%s: %s',fn,fn,getReport(ME));
    end
  end
  
end

handles = handles_save;

try
  save(filename,'handles');
  success = true;
  msg = '';
catch ME
  success = false;
  msg = getReport(ME);
end


% ---
function figure_CloseRequestFcn(hObject, eventdata)

handles=guidata(hObject);

if ~isempty(handles.system_monitor),
  try
    stop(handles.system_monitor.timer);
    delete(handles.system_monitor.timer);
  catch ME,
    fprintf('Error stopping system monitor:\n%s\n',getReport(ME));
  end
end
handles.system_monitor=[];

filename = handles.rcfilename;

try
  [success,msg] = SaveConfiguration(handles,filename);
  if ~success,
    error(msg);
  end
catch ME,
  warning(getReport(ME));
  warndlg('Could not save the last configuration. State will not be saved.');
end
delete(hObject);


% --- Outputs from this function are returned to the command line.
function varargout = JAABAPlot_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% ---
function ExperimentListCallback(hObject,eventdata)

handles=guidata(hObject);
handles.experimentvalue{handles.groupvalue}=get(handles.ExperimentList,'Value');
handles.interestingfeaturehistograms_cache=[];
handles.interestingfeaturetimeseries_cache=[];
guidata(hObject,handles);


% ---
function ClassifierListCallback(hObject,eventdata)

handles=guidata(hObject);
handles.classifiervalue=get(handles.ClassifierList,'Value');
handles.interestingfeaturehistograms_cache=[];
handles.interestingfeaturetimeseries_cache=[];
guidata(hObject,handles);


% --- Executes on button press in UpdatePlot.
function UpdatePlot_Callback(hObject, eventdata, handles)
% hObject    handle to UpdatePlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)


% --- Executes on selection change in IndividualList.
function IndividualList_Callback(hObject, eventdata, handles)
% hObject    handle to IndividualList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns IndividualList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from IndividualList

handles.individualvalue=get(handles.IndividualList,'Value');
tmp=cellfun(@(x) x(1),get(handles.IndividualList,'String'));
handles.individualidx=tmp(handles.individualvalue);
if handles.individualidx=='G'
  handles.individualidx=handles.individualvalue-find(tmp=='G',1,'first')+1;
end
if ((handles.individualidx=='M') && length(handles.individuallist{handles.individualvalue})>5)
  handles.individualidx={'M' 'F'};
end
handles.interestingfeaturehistograms_cache=[];
handles.interestingfeaturetimeseries_cache=[];
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function IndividualList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IndividualList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in FeatureList.
function FeatureList_Callback(hObject, eventdata, handles)
% hObject    handle to FeatureList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns FeatureList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from FeatureList

handles.featurevalue=get(handles.FeatureList,'Value');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function FeatureList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FeatureList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in BehaviorNot.
function BehaviorNot_Callback(hObject, eventdata, handles)
% hObject    handle to BehaviorNot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of BehaviorNot

handles.behaviornot=get(handles.BehaviorNot,'Value');
guidata(hObject,handles);


% --- Executes on selection change in BehaviorList.
function BehaviorList_Callback(hObject, eventdata, handles)
% hObject    handle to BehaviorList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns BehaviorList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from BehaviorList

handles.behaviorvalue=get(handles.BehaviorList,'Value');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function BehaviorList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BehaviorList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in BehaviorLogic.
function BehaviorLogic_Callback(hObject, eventdata, handles)
% hObject    handle to BehaviorLogic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns BehaviorLogic contents as cell array
%        contents{get(hObject,'Value')} returns selected item from BehaviorLogic

handles.behaviorlogic=get(handles.BehaviorLogic,'Value');
if(handles.behaviorlogic==1)
  set(handles.BehaviorList2,'enable','off');
else
  set(handles.BehaviorList2,'enable','on');
end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function BehaviorLogic_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BehaviorLogic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in BehaviorList2.
function BehaviorList2_Callback(hObject, eventdata, handles)
% hObject    handle to BehaviorList2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns BehaviorList2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from BehaviorList2

handles.behaviorvalue2=get(handles.BehaviorList2,'Value');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function BehaviorList2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BehaviorList2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in BehaviorNormalizeNot.
function BehaviorNormalizeNot_Callback(hObject, eventdata, handles)
% hObject    handle to BehaviorNormalizeNot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of BehaviorNormalizeNot

handles.behaviornormalizenot=get(handles.BehaviorNormalizeNot,'Value');
guidata(hObject,handles);


% --- Executes on selection change in BehaviorList3.
function BehaviorList3_Callback(hObject, eventdata, handles)
% hObject    handle to BehaviorList3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns BehaviorList3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from BehaviorList3

handles.behaviorvalue3=get(handles.BehaviorList3,'Value');
if(handles.behaviorvalue3==1)
  handles.behaviornormalizenot=0;
  set(handles.BehaviorNormalizeNot,'enable','off');
  update_figure(handles);
else
  set(handles.BehaviorNormalizeNot,'enable','on');
end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function BehaviorList3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BehaviorList3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ExperimentList.
function ExperimentList_Callback(hObject, eventdata, handles)
% hObject    handle to ExperimentList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ExperimentList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ExperimentList



% --- Executes during object creation, after setting all properties.
function ExperimentList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ExperimentList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% ---
function feature_intersection=check_for_diff_and_return_intersection(arg)

if(length(arg)==0)  feature_intersection=[];  return;  end
if(length(arg)==1)  feature_intersection=arg{1};  return;  end

feature_union=unique([arg{:}]);
feature_intersection=arg{1};
for i=2:length(arg)
  feature_intersection=intersect(feature_intersection,arg{i});
end
if(numel(feature_intersection)<numel(feature_union))
  setdiff(feature_union,feature_intersection);
  cellfun(@(x) [char(x) ', '],ans,'uniformoutput',false);
  uiwait(warndlg([{'the following features are not in all experiments and so will be ignored.' '' [ans{:}]}]));
  drawnow;
end


% ---
function handles=fillin_individuallist(handles)

%if((numel(handles.individuals_behavior)==0) || ...
%    (sum(sum(diff(handles.individuals_behavior,[],2)~=0))>0))  return;  end

cumsum_num_exp_per_group=[0 cumsum(cellfun(@length,handles.experimentlist))];

tmp(1)={'All'};
if(sum(cellfun(@(x) islogical([x{:}]),handles.sexdata)==0)==0)
  tmp(2:4)={'Male' 'Female' 'Male vs. Female'};
end
k=1;
for e=1:length(handles.individuals_feature)
  g=find(cumsum_num_exp_per_group<e,1,'last');
  for i=1:handles.individuals_feature(e)
    tmp{end+1}=['Grp ' handles.grouplist{g} ', Exp #' num2str(e-cumsum_num_exp_per_group(g)) ', Indi #' num2str(i)];
    k=k+1;
  end
end
handles.individuallist=tmp;
if(handles.individualvalue>length(handles.individuallist))
  handles.individualvalue=1;
  handles.individualidx='A';
end
set(handles.IndividualList,'String',handles.individuallist,'Value',handles.individualvalue);


% ---
function handles=experiment_add(handles,newexperiments,newgroups,newcolors)

newexperiments=[cellfun(@(x) regexprep(x,'/$',''),newexperiments,'uniformoutput',false)];
tmp=ismember(newexperiments,[handles.experimentlist{:}]);
if(sum(tmp)>0)
  msg{1}='The following experiments have already been added:';
  msg{2}='';
  msg(3:(2+sum(tmp)))=newexperiments(tmp);
  uiwait(errordlg(msg));
  newexperiments(tmp)=[];
  if(~isempty(newgroups))  newgroups(tmp)=[];  end
  if(~isempty(newcolors))  newcolors(tmp)=[];  end
end
if isempty(newexperiments),
  return;
end

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

if(~exist(fullfile(newexperiments{1},handles.perframe_dir),'dir'))

  [~,name] = fileparts(newexperiments{1});
  res = questdlg(sprintf('Cannot find per-frame directory in default location %s, Generate per-frame data, Locate per-frame directory, or Cancel?',fullfile(name,'perframe')),'Per-Frame Directory Not Found','Generate Per-Frame Data','Locate Per-Frame Directory','Cancel','Cancel');
  if strcmp(res,'Cancel'),
    set(handles.Status,'string','Ready.','foregroundcolor','g');
    set(handles.figure1,'pointer','arrow');
    drawnow;
    return;
  elseif strcmp(res,'Generate Per-Frame Data'),
    if ~(isfield(handles,'perframe_dir') && ischar(handles.perframe_dir) && ~isempty(handles.perframe_dir)),
      handles.perframe_dir = 'perframe';
    end
    if ~(isfield(handles,'trx_file') && ischar(handles.trx_file) && ~isempty(handles.trx_file)) || ...
        ~exist(fullfile(newexperiments{1},handles.trx_file),'file'),
      filecurr = uigetfile(fullfile(newexperiments{1},'*.mat'),sprintf('Select trx file for %s',newexperiments{i}));
      if ~ischar(filecurr),
        set(handles.Status,'string','Ready.','foregroundcolor','g');
        set(handles.figure1,'pointer','arrow');
        drawnow;
        return;
      end
      handles.trx_file = filecurr;
    end 
    perframetrx = Trx('trxfilestr',handles.trx_file,...
      'perframedir',handles.perframe_dir);
    perframetrx.AddExpDir(newexperiments{1},'dooverwrite',false,'openmovie',false);
  else
    res = uigetdir2(newexperiments{1},...
      'Select per-frame directory');
    if isempty(res),
      set(handles.Status,'string','Ready.','foregroundcolor','g');
      set(handles.figure1,'pointer','arrow');
      drawnow;
      return;
    end
    [~,handles.perframe_dir,~]=fileparts(res);
  end
end

if(isnan(handles.fps))
%if((isnan(handles.fps))&&(length(handles.classifierlist)>0))
%  classifier=load(handles.classifierlist{1});
%  handles.fps=get_fps(fullfile(newexperiments{1},classifier.trxfilename));
  [handles.fps,handles.trx_file]=get_fps(newexperiments{1},handles.trx_file);
  if(isempty(handles.trx_file))
    set(handles.Status,'string','Ready.','foregroundcolor','g');
    set(handles.figure1,'pointer','arrow');
    drawnow;
    return;
  end
end

handlesfeatures=cell(1,length(newexperiments));
handlessexdata=cell(1,length(newexperiments));
handlesindividualsbehavior=zeros(length(newexperiments),length(handles.scorefiles));
handlesindividualsfeature=zeros(1,length(newexperiments));
idx_error=false(1,length(newexperiments));
parfor n=1:length(newexperiments)
%for n=1:length(newexperiments)
  tmp=dir(fullfile(newexperiments{n},handles.perframe_dir,'*.mat'));
  if isempty(tmp),
    handlesfeatures{n} = {};
  else
    [handlesfeatures{n}{1:length(tmp)}]=deal(tmp.name);
    handlesfeatures{n}=cellfun(@(x) x(1:(end-4)),handlesfeatures{n},'uniformoutput',false);
  end

  fullfile(newexperiments{n},handles.perframe_dir,'sex.mat');
  if(exist(ans,'file'))
    tmp=load(ans);
    cellfun(@(x) strcmp(x,'M'),tmp.data,'uniformoutput',false);
    handlessexdata(n)={ans};
  else
    tmp=dir(fullfile(newexperiments{n},handles.perframe_dir,'*.mat'));
    if(~isempty(tmp))
      tmp=load(fullfile(newexperiments{n},handles.perframe_dir,tmp(1).name));
      cellfun(@(x) nan(1,length(x)),tmp.data,'uniformoutput',false);
      handlessexdata(n)={ans};
    else
      idx_error(n)=true;
      continue;
    end
  end

  behavior_data=[];
  parfor_tmp=zeros(1,length(handles.scorefiles));
  for s=1:length(handles.scorefiles)
    %classifier=load();
    classifier=load(handles.classifierlist{s},'-mat');
    %classifier=x;
    parfor_tmp(s)=get_nindividuals_behavior(newexperiments{n},handles.scorefiles{s},...
        classifier.x.classifierStuff.timeStamp);
  end
  handlesindividualsbehavior(n,:)=parfor_tmp;

  handlesindividualsfeature(n)=...
      get_nindividuals_feature(newexperiments{n},handles.perframe_dir,handlesfeatures{n});
end

if(sum(idx_error)>0)
  uiwait(errordlg([{'skipping the following experiments because they have not been tracked:  ' '' ...
      newexperiments{idx_error}}],''));  drawnow;
  newexperiments=newexperiments(~idx_error);
  handlesfeatures=handlesfeatures(~idx_error);
  handlessexdata=handlessexdata(~idx_error);
  if(~isempty(handlesindividualsbehavior))  handlesindividualsbehavior=handlesindividualsbehavior(~idx_error);  end
  handlesindividualsfeature=handlesindividualsfeature(~idx_error);
  if(~isempty(newgroups))  newgroups=newgroups(~idx_error);  end
  if(~isempty(newcolors))  newcolors=newcolors(~idx_error);  end
end
if isempty(newexperiments),
  set(handles.Status,'string','Ready.','foregroundcolor','g');
  set(handles.figure1,'pointer','arrow');
  drawnow;
  return;
end

if(isempty(newgroups))
  handles.experimentlist{handles.groupvalue}={handles.experimentlist{handles.groupvalue}{:} newexperiments{:}};
  handles.experimentvalue{handles.groupvalue}=1:length(handles.experimentlist{handles.groupvalue});
else
%  [newgroups,i]=sort(newgroups);
%  newexperiments=newexperiments(i);
%  if(~isempty(newcolors))  newcolors=newcolors(i);  end
  ridx=[];
  [unique_groups,ia,ic]=unique(newgroups,'stable');
  for ugi=1:length(unique_groups)
    k=length(handles.grouplist);
    handles.grouplist{k+1}=unique_groups{ugi};
    if(~isempty(newcolors))
      handles.colors(k+1,1)=hex2dec(newcolors{ia(ugi)}(1:2))/255;
      handles.colors(k+1,2)=hex2dec(newcolors{ia(ugi)}(3:4))/255;
      handles.colors(k+1,3)=hex2dec(newcolors{ia(ugi)}(5:6))/255;
    else
      handles.colors(k+1,:)=[0 0 0];
    end
    find(cellfun(@(x) strcmp(x,unique_groups{ugi}),newgroups));
    handles.experimentlist{k+1}={newexperiments{ans}};
    handles.experimentvalue{k+1}=1:length(handles.experimentlist{k+1});
    ridx=[ridx; ans];
  end
  handlesfeatures=handlesfeatures(ridx);
  handlessexdata=handlessexdata(ridx);
  handlesindividualsbehavior=handlesindividualsbehavior(ridx,:);
  handlesindividualsfeature=handlesindividualsfeature(ridx);
end

handles.features={handles.features{:} handlesfeatures{:}};
handles.sexdata={handles.sexdata{:} handlessexdata{:}};
handles.individuals_behavior=[handles.individuals_behavior; handlesindividualsbehavior];
handles.individuals_feature=[handles.individuals_feature handlesindividualsfeature];

if(isempty(newgroups))
  tmp=length(handles.features);
  idx=[1 : (sum(cellfun(@length,handles.experimentlist(1:handles.groupvalue)))-length(newexperiments)) ...
      ((tmp-length(newexperiments)+1) : tmp) ...
      ((tmp-length(newexperiments)-sum(cellfun(@length,handles.experimentlist((handles.groupvalue+1):end)))+1) : ...
          (tmp-length(newexperiments)))];
  handles.features=handles.features(idx);
  handles.sexdata=handles.sexdata(idx);
  handles.individuals_behavior=handles.individuals_behavior(idx,:);
  handles.individuals_feature=handles.individuals_feature(idx);
end

handles.featurelist=check_for_diff_and_return_intersection(handles.features);
handles=fillin_individuallist(handles);

handles.interestingfeaturehistograms_cache=[];
handles.interestingfeaturetimeseries_cache=[];

update_figure(handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% --- Executes on button press in ExperimentBatch.
function ExperimentBatch_Callback(hObject, eventdata, handles)
% hObject    handle to ExperimentBatch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

persistent directory
if(isempty(directory))  directory=pwd;  end

tmp=directory;
[newexperiments directory]=uigetfile([fullfile(directory,'*.txt') ';*.csv'],'Select batch file');
if(isnumeric(newexperiments)&&(newexperiments==0))  directory=tmp; return;  end
fid=fopen(fullfile(directory,newexperiments),'r');
textscan(fgetl(fid),'%s','delimiter',',');
ncol=numel(ans{1});
fclose(fid);
newexperiments=textread(fullfile(directory,newexperiments),'%s','delimiter',',');
%sum(cellfun(@(x) exist(x,'dir'),newexperiments)~=0);
%if(ans==(length(newexperiments)/2))
switch(ncol)
  case 1
    newgroups=[];
    newcolors=[];
  case 2
    newgroups=newexperiments(2:2:end);
    newcolors=[];
    newexperiments=newexperiments(1:2:end);
  case 3
    newgroups=newexperiments(2:3:end);
    newcolors=newexperiments(3:3:end);
    newexperiments=newexperiments(1:3:end);
end
find(cellfun(@(x) ~exist(x,'dir'),newexperiments));
if(~isempty(ans))
  uiwait(errordlg({'These experiments were not found:','',newexperiments{ans}}));
  return;
end
if((length(handles.grouplist)==0)&&(isempty(newgroups)))
  uiwait(errordlg('Add a new group before adding ungrouped experiments'));
  return;
end
uniqueexperiments=unique(newexperiments);
if(length(uniqueexperiments)~=length(newexperiments))
  tmp=logical(cellfun(@(x) sum(strcmp(x,newexperiments)),uniqueexperiments)>1);
  msg{1}='The following experiments are duplicated:';
  msg{2}='';
  msg(3:(2+sum(tmp)))=uniqueexperiments(tmp);
  uiwait(errordlg(msg));
  for i=uniqueexperiments(tmp)
    foo=find(strcmp(i,newexperiments));
    foo=foo(1:end-1);
    newexperiments(foo)=[];
    if(~isempty(newgroups))  newgroups(foo)=[];  end
    if(~isempty(newcolors))  newcolors(foo)=[];  end
  end
end

handles=experiment_add(handles,newexperiments,newgroups,newcolors);

guidata(hObject,handles);


% --- Executes on button press in ExperimentAdd.
function ExperimentAdd_Callback(hObject, eventdata, handles)
% hObject    handle to ExperimentAdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

persistent directory
if(isempty(directory))  directory=pwd;  end

newexperiments=uipickfiles('prompt','Select experiment directory','filterspec',directory);
if(~iscell(newexperiments) || (length(newexperiments)==0))  return;  end
if(sum(~cellfun(@(x) exist(x,'dir'),newexperiments))>0)
  uiwait(errordlg('please select only experimental directories'));  drawnow;
  return;
end
[directory,~,~]=fileparts(newexperiments{1});

handles=experiment_add(handles,newexperiments,[],[]);

guidata(hObject,handles);


% --- Executes on button press in ExperimentDelete.
function ExperimentDelete_Callback(hObject, eventdata, handles)
% hObject    handle to ExperimentDelete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

if(length(handles.experimentlist)==0)  return;  end

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

idx=handles.experimentvalue{handles.groupvalue};
handles.experimentlist{handles.groupvalue}(idx)=[];
handles.experimentvalue{handles.groupvalue}=1:length(handles.experimentlist{handles.groupvalue});

idx=idx+sum(cellfun(@length,handles.experimentlist(1:(handles.groupvalue-1))));
handles.features(idx)=[];
handles.sexdata(idx)=[];
handles.individuals_behavior(idx,:)=[];
handles.individuals_feature(idx)=[];

if(isempty(handles.experimentlist{handles.groupvalue}))
  handles.colors=handles.colors(setdiff(1:size(handles.colors,1),handles.groupvalue),:);
  handles.experimentlist(handles.groupvalue)=[];
  handles.experimentvalue(handles.groupvalue)=[];
  handles.grouplist(handles.groupvalue)=[];
  handles.groupvalue=max(1,min([handles.groupvalue length(handles.grouplist)]));
end

if(isempty(handles.experimentlist))
  handles.grouplist={};
  handles.groupvalue=1;
  handles.experimentlist={{}};
  handles.experimentvalue={1};
end

handles.featurelist=check_for_diff_and_return_intersection(handles.features);
handles=fillin_individuallist(handles);

update_figure(handles);

handles.interestingfeaturehistograms_cache=[];
handles.interestingfeaturetimeseries_cache=[];

guidata(hObject,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% --- Executes on button press in ExperimentMove.
function ExperimentMove_Callback(hObject, eventdata, handles)
% hObject    handle to ExperimentMove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(isempty(handles.experimentlist{handles.groupvalue}))  return;  end
[to_group,ok]=listdlg('liststring',...
    {handles.grouplist{setdiff(1:length(handles.grouplist),handles.groupvalue)}},...
    'selectionmode','single');
if(~ok)  return;  end
if(to_group>=handles.groupvalue)  to_group=to_group+1;  end

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

from_group=handles.groupvalue;
idx=handles.experimentvalue{from_group};
idxF=idx+sum(cellfun(@length,handles.experimentlist(1:(from_group-1))));
idxT=sum(cellfun(@length,handles.experimentlist(1:to_group)));

handles.experimentlist{to_group}=...
    {handles.experimentlist{to_group}{:} handles.experimentlist{from_group}{idx}};
handles.experimentvalue{to_group}=1:length(handles.experimentlist{to_group});
handles.experimentlist{from_group}(idx)=[];
handles.experimentvalue{from_group}=1:length(handles.experimentlist{from_group});
%if(isempty(handles.experimentlist{from_group}))
%  set(handles.ExperimentList,'String',{''},'Value',1);
%else
%  set(handles.ExperimentList,'String',handles.experimentlist{from_group});
%  set(handles.ExperimentList,'Value',handles.experimentvalue{from_group});
%end
if(idxF(1)<=idxT)  idxT=idxT-length(idxF);  end
tmp=setdiff(1:length(handles.features),idxF);
tmp=[tmp(1:idxT) idxF tmp((idxT+1):end)];
%handles.behaviors=handles.behaviors(tmp);
handles.features=handles.features(tmp);
handles.sexdata=handles.sexdata(tmp);
handles.individuals_behavior=handles.individuals_behavior(tmp,:);
handles.individuals_feature=handles.individuals_feature(tmp);

handles=fillin_individuallist(handles);

update_figure(handles);

handles.interestingfeaturehistograms_cache=[];
handles.interestingfeaturetimeseries_cache=[];

guidata(hObject,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% --- Executes on selection change in GroupList.
function GroupList_Callback(hObject, eventdata, handles)
% hObject    handle to GroupList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns GroupList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from GroupList

handles.groupvalue=get(handles.GroupList,'Value');
set(handles.ExperimentList,'String',handles.experimentlist{handles.groupvalue},...
    'Value',handles.experimentvalue{handles.groupvalue});
update_figure(handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function GroupList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GroupList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in GroupNew.
function GroupNew_Callback(hObject, eventdata, handles)
% hObject    handle to GroupNew (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

inputdlg({'Name:'},'Create new experiment group');
if((isempty(ans))||strcmp(ans,''))  return;  end
handles.colors(end+1,:)=...
    uisetcolor(handles.defaultcolors(1+mod(length(handles.grouplist),size(handles.defaultcolors,1)),:));
handles.grouplist{end+1}=char(ans);
handles.groupvalue=length(handles.grouplist);
handles.experimentlist{handles.groupvalue}={};
handles.experimentvalue{handles.groupvalue}=[];
update_figure(handles);
guidata(hObject,handles);


% --- Executes on button press in GroupChange.
function GroupChange_Callback(hObject, eventdata, handles)
% hObject    handle to GroupChange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

inputdlg({'Name:'},'New name for this group',1,handles.grouplist(handles.groupvalue));
if((isempty(ans))||strcmp(ans,''))  return;  end
handles.colors(handles.groupvalue,:)=uisetcolor(handles.colors(handles.groupvalue,:));
handles.grouplist{handles.groupvalue}=char(ans);
update_figure(handles);
guidata(hObject,handles);


% ---
function handles=classifier_add(handles,newclassifiers)

%handlesconfigurations=cell(1,length(newclassifiers));
handlesbehaviorlist=cell(1,length(newclassifiers));
handlesscorefiles=cell(1,length(newclassifiers));
handlesindividualsbehavior=zeros(sum(cellfun(@length,handles.experimentlist)),length(newclassifiers));
timestamps=nan(1,length(newclassifiers));
parfor c=1:length(newclassifiers)
%for c=1:length(newclassifiers)
  classifier=load(newclassifiers{c},'-mat'); %#ok<PFTUS>

  timestamps(c)=classifier.x.classifierStuff.timeStamp;

  handlesbehaviorlist{c}=classifier.x.behaviors.names{1};
  %handlesscorefiles{c}=params.file.scorefilename;
  handlesscorefiles{c}=classifier.x.file.scorefilename;

  ee=0;  behavior_data=[];  parfor_tmp=zeros(sum(cellfun(@length,handles.experimentlist)),1);
  for g=1:length(handles.grouplist)
    for e=1:length(handles.experimentlist{g})
      parfor_tmp(ee+e)=get_nindividuals_behavior(handles.experimentlist{g}{e},handlesscorefiles{c},...
          classifier.x.classifierStuff.timeStamp);
    end
    ee=ee+length(handles.experimentlist{g});
  end
  handlesindividualsbehavior(:,c)=parfor_tmp;
end

handlesexperimentlist=[handles.experimentlist{:}];
if((isnan(handles.fps))&&(length(handlesexperimentlist)>0))
  handles.trx_file=classifier.x.file.trxfilename;
  [handles.fps,handles.trx_file]=get_fps(handlesexperimentlist{1},handles.trx_file);
end

idx=find(~cellfun(@isempty,newclassifiers));
newclassifiers=newclassifiers(idx);
%handlesconfigurations=handlesconfigurations(idx);
handlesbehaviorlist=handlesbehaviorlist(idx);
handlesscorefiles=handlesscorefiles(idx);
handlesindividualsbehavior=handlesindividualsbehavior(:,idx);
timestamps=timestamps(idx);

i=1;
while i<length(handlesscorefiles)
  idx=find(strcmp(handlesscorefiles{i},handlesscorefiles((i+1):end)) & ...
      (timestamps(i)==timestamps((i+1):end)));
  if(~isempty(idx))
    uiwait(warndlg({'the following classifiers have the same scores filename and timestamp.  keeping only the first.' '' newclassifiers{[i i+idx]}},''));  drawnow;
    newclassifiers(i+idx)=[];
%    handlesconfigurations(i+idx)=[];
    handlesbehaviorlist(i+idx)=[];
    handlesscorefiles(i+idx)=[];
    handlesindividualsbehavior(:,i+idx)=[];
    timestamps(i+idx)=[];
  end
  i=i+1;
end

handles.classifierlist={handles.classifierlist{:} newclassifiers{:}};
%handles.configurations={handles.configurations{:} handlesconfigurations{:}};
handles.behaviorlist={handles.behaviorlist{:} handlesbehaviorlist{:}};
handles.scorefiles={handles.scorefiles{:} handlesscorefiles{:}};
handles.individuals_behavior=[handles.individuals_behavior handlesindividualsbehavior];

handles.classifiervalue=1:length(handles.classifierlist);

handles=fillin_individuallist(handles);

handles.interestingfeaturehistograms_cache=[];
handles.interestingfeaturetimeseries_cache=[];


% --- Executes on button press in ClassifierBatch.
function ClassifierBatch_Callback(hObject, eventdata, handles)
% hObject    handle to ClassifierBatch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

persistent directory
if(isempty(directory))  directory=pwd;  end

tmp=directory;
[newclassifiers directory]=uigetfile(fullfile(directory,'*.csv'),'Select classifier files');
if(isnumeric(newclassifiers)&&(newclassifiers==0))  directory=tmp; return;  end
newclassifiers=textread(fullfile(directory,newclassifiers),'%s','delimiter',',');

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

handles=classifier_add(handles,newclassifiers);
update_figure(handles);

guidata(hObject,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% ---
function delete_classifer_check_table()

hh=findobj('type','figure');
for h=1:length(hh)
  if(strcmp(get(hh(h),'name'),'classifier check'))
    close(hh(h));
  end
end


% --- Executes on button press in ClassifierAdd.
function ClassifierAdd_Callback(hObject, eventdata, handles)
% hObject    handle to ClassifierAdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

persistent directory
if(isempty(directory))  directory=pwd;  end

tmp=directory;
[newclassifiers directory]=uigetfile(fullfile(directory,'*.jab'),...
    'Select classifier files','multiselect','on');
if(isnumeric(newclassifiers)&&(newclassifiers==0))  directory=tmp; return;  end
if(~iscell(newclassifiers))  newclassifiers={newclassifiers};  end
newclassifiers=cellfun(@(x) fullfile(directory,x),newclassifiers,'uniformoutput',false);

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

handles=classifier_add(handles,newclassifiers);
update_figure(handles);

delete_classifer_check_table();

guidata(hObject,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% --- Executes on button press in ClassifierDelete.
function ClassifierDelete_Callback(hObject, eventdata, handles)
% hObject    handle to ClassifierDelete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

idx=handles.classifiervalue;
handles.classifierlist(idx)=[];
handles.classifiervalue=1:max(1,length(handles.classifierlist));
%handles.configurations(idx)=[];
handles.behaviorlist(idx)=[];
handles.behaviorvalue=1;
handles.behaviorvalue2=1;
handles.behaviorvalue3=1;
handles.scorefiles(idx)=[];
handles.individuals_behavior(:,idx)=[];

%if(length(handles.classifierlist)==0)
%  handles.fps=nan;
%end

handles=fillin_individuallist(handles);
update_figure(handles);

handles.interestingfeaturehistograms_cache=[];
handles.interestingfeaturetimeseries_cache=[];

delete_classifer_check_table();

guidata(hObject,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% --- Executes on button press in ClassifierAuto.
function ClassifierAuto_Callback(hObject, eventdata, handles)
% hObject    handle to ClassifierAuto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

handlesexperimentlist=[handles.experimentlist{:}];

classifiers_found=cell(1,length(handlesexperimentlist));
classifiers_notfound=cell(1,length(handlesexperimentlist));
parfor ge=1:length(handlesexperimentlist)
%for ge=1:length(handlesexperimentlist)
  tmp=dir(fullfile(handlesexperimentlist{ge},'*.mat'));
  possiblescorefiles=setdiff({tmp.name},handles.scorefiles);
  classifiers_found{ge}={};
  classifiers_notfound{ge}={};
  for p=1:length(possiblescorefiles)
    %scoresfilename = fullfile(handlesexperimentlist{ge},possiblescorefiles{p});
    %if ~exist(scoresfilename,'file'),
    %  continue;
    %end
    tmp=load(fullfile(handlesexperimentlist{ge},possiblescorefiles{p}));
    if(~isfield(tmp,'jabFileNameAbs'))
      continue;
    end
    if(exist(tmp.jabFileNameAbs)==2)
      classifiers_found{ge}{end+1}=tmp.jabFileNameAbs;
    else
      classifiers_notfound{ge}{end+1}=tmp.jabFileNameAbs;
    end
  end
end
classifiers_found=unique([classifiers_found{:}]);
classifiers_notfound=unique([classifiers_notfound{:}]);

classifiers_found=setdiff(classifiers_found,handles.classifierlist);
handles=classifier_add(handles,classifiers_found);

if(length(classifiers_notfound)>2)
  uiwait(errordlg({'Could not find these classifiers:' '' [classifiers_notfound{:}]}));
end

update_figure(handles);

delete_classifer_check_table();

guidata(hObject,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% --- Executes on button press in ClassifierCheck.
function ClassifierCheck_Callback(hObject, eventdata, handles)
% hObject    handle to ClassifierCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

handlesexperimentlist=[handles.experimentlist{:}];

table=cell(length(handlesexperimentlist),length(handles.scorefiles));

parfor ge=1:length(handlesexperimentlist)
%for ge=1:length(handlesexperimentlist)
  behavior_data=[];
  parfor_tmp=cell(1,length(handles.scorefiles));
  for c=1:length(handles.classifierlist)
    classifier=load(handles.classifierlist{c},'-mat');
    %classifier=x;
    try
      behavior_data=load(fullfile(handlesexperimentlist{ge},handles.scorefiles{c}));
    catch ME,
%      warning(getReport(ME));
      parfor_tmp{c}='missing';
      continue;
    end
    if (classifier.x.classifierStuff.timeStamp ~= behavior_data.timestamp)
      parfor_tmp{c}=[datestr(classifier.x.classifierStuff.timeStamp) ' ~= ' datestr(behavior_data.timestamp)];
    else
      parfor_tmp{c}='up-to-date';
    end
  end
  table(ge,:)=parfor_tmp;
  tmp=dir(fullfile(handlesexperimentlist{ge},'*.mat'));
  extrascorefiles=setdiff({tmp.name},handles.scorefiles);
  table2{ge}={};
  for s=1:length(extrascorefiles)
    tmp=load(fullfile(handlesexperimentlist{ge},extrascorefiles{s}));
    if(~isfield(tmp,'jabFileNameAbs'))
      continue;
    end
    table2{ge}{end+1}=extrascorefiles{s};
  end
end

[~,tmp,~]=cellfun(@(x) fileparts(x),handlesexperimentlist,'uniformoutput',false);
table=[tmp' table];
[~,tmp,~]=cellfun(@(x) fileparts(x),handles.classifierlist,'uniformoutput',false);
table=[{'' tmp{:}}; table];

tmp=unique([table2{:}]);
if(length(tmp)>0)
  table{1,end+1}='   ';
  table(1,(end+1):(end+length(tmp)))=tmp;
end
for i=1:length(table2)
  for j=1:length(table2{i})
    idx=find(cellfun(@(x) strcmp(x,table2{i}{j}),table(1,:)));
    table{i+1,idx}='extra';
  end
end

if(length(handles.experimentlist)>1)
  tmp=zeros(1,1+length([handles.experimentlist{:}]));
  cumsum(cellfun(@length,handles.experimentlist));
  tmp(2+ans(1:end-1))=1;
  cumsum(tmp);
  tmp=ans+(1:(1+length([handles.experimentlist{:}])));
  tmp2(tmp,:)=[table(:,:)];
  table=tmp2;
end

hf=figure('menubar','none','toolbar','none','numbertitle','off',...
    'name','classifier check');
ht=uitable('data',table,'columnwidth',num2cell(8*max(cellfun(@length,table))),...
    'rowname',[],'columnname',{''},'rowstriping','off');
extT=get(ht,'extent');
posF=get(hf,'position');
set(hf,'position',[posF(1) posF(2) min(posF(3),16*(posF(4)<extT(4))+extT(3)) min(posF(4),extT(4))]);
set(ht,'units','normalized','position',[0 0 1 1]);

if(handles.dump2csv)
  fid=fopen('most_recent_table.csv','w');
  transpose(table);
  fprintf(fid,[repmat('%s, ',1,size(table,2)) '\n'],ans{:});
  fclose(fid);
end

if(handles.dump2mat)
  save('most_recent_table.mat','table');
end

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% --- Executes on button press in ClassifierClassify.
function ClassifierClassify_Callback(hObject, eventdata, handles)
% hObject    handle to ClassifierClassify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

%h=waitbar(0,'This will likely take awhile...',...
%    'CreateCancelBtn','fid=fopen(fullfile(tempdir,''cancel.txt''),''w''); fclose(fid);');
h=waitbar(0,'This will likely take awhile...');
fid=fopen(fullfile(tempdir,'progressbar.txt'),'w');
fwrite(fid,length(handles.behaviorlist)*sum(cellfun(@length,handles.experimentvalue)),'uint32');
fclose(fid);
t = timer('TimerFcn',{@progress_bar,h}, 'Period', 3, 'ExecutionMode', 'fixedRate');
start(t);

handlesexperimentlist=[handles.experimentlist{:}];

for ge=1:length(handlesexperimentlist)
%for ge=1:length(handlesexperimentlist)
%   JAABADetect(handlesexperimentlist{ge},...
%       'classifierfiles',handles.classifierlist(handles.classifiervalue),...
%       'forcecompute',handles.classify_forcecompute);
  JAABADetect(handlesexperimentlist{ge},...
      'jabfiles',handles.classifierlist(handles.classifiervalue),...
      'forcecompute',handles.classify_forcecompute);
%      'configfiles',handles.configurations(handles.classifiervalue),...
  fid=fopen(fullfile(tempdir,'progressbar.txt'),'a');  fwrite(fid,1,'uint32');  fclose(fid);
end

stop(t);
delete(t);
delete(h);
delete(fullfile(tempdir,'progressbar.txt'));

%if(exist(fullfile(tempdir,'cancel.txt')))
%  delete(fullfile(tempdir,'cancel.txt'));
%  set(handles.Status,'string','Ready.','foregroundcolor','g');
%  drawnow;
%  set(handles.figure1,'pointer','arrow');
%  return;
%end

handles=update_experiment_data(handles,true,true,true);
update_figure(handles);

handles.interestingfeaturehistograms_cache=[];

guidata(hObject,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% --- Executes on selection change in ClassifierList.
function ClassifierList_Callback(hObject, eventdata, handles)
% hObject    handle to ClassifierList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ClassifierList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ClassifierList


% --- Executes during object creation, after setting all properties.
function ClassifierList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ClassifierList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in FeatureHistogram.
function FeatureHistogram_Callback(hObject, eventdata, handles)
% hObject    handle to FeatureHistogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

handles.analysis='feature_histogram';
update_figure(handles);
guidata(hObject,handles);


% --- Executes on button press in FeatureTimeSeries.
function FeatureTimeSeries_Callback(hObject, eventdata, handles)
% hObject    handle to FeatureTimeSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

handles.analysis='feature_timeseries';
update_figure(handles);
guidata(hObject,handles);


% --- Executes on button press in BehaviorBarChart.
function BehaviorBarChart_Callback(hObject, eventdata, handles)
% hObject    handle to BehaviorBarChart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.analysis='behavior_barchart';
handles.behaviorvalue = min(handles.behaviorvalue, length(handles.behaviorlist)+1);
update_figure(handles);
guidata(hObject,handles);


% --- Executes on button press in BehaviorTimeSeries.
function BehaviorTimeSeries_Callback(hObject, eventdata, handles)
% hObject    handle to BehaviorTimeSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.analysis='behavior_timeseries';
handles.behaviorvalue = min(handles.behaviorvalue, length(handles.behaviorlist)+1);
update_figure(handles);
guidata(hObject,handles);

%
% --- Executes on button press in BoutStats.
function BoutStats_Callback(hObject, eventdata, handles)
% hObject    handle to BoutStats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.analysis='bout_stats';
handles.behaviorvalue = min(handles.behaviorvalue, length(handles.behaviorlist)+1);
update_figure(handles);
guidata(hObject,handles);


% --- Executes on button press in InterestingFeatureHistograms.
function InterestingFeatureHistograms_Callback(hObject, eventdata, handles)
% hObject    handle to InterestingFeatureHistograms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

handles.analysis='interesting_feature_histograms';
handles.individualvalue=min(3,handles.individualvalue);
tmp=cellfun(@(x) x(1),get(handles.IndividualList,'String'));
handles.individualidx=tmp(handles.individualvalue);
handles.behaviorvalue = min(handles.behaviorvalue, length(handles.behaviorlist)+1);
update_figure(handles);
guidata(hObject,handles);


% --- Executes on button press in Plot.
function Plot_Callback(hObject, eventdata, handles)
% hObject    handle to Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(handles.behaviorvalue>length(handles.behaviorlist))
  if((sum(sum(handles.individuals_behavior==-1))>0) || ...
     (sum(sum(diff(handles.individuals_behavior,[],2)~=0))>0))
    uiwait(errordlg('some experiments either have no individiuals or differ in the number of individuals across classifiers'));  drawnow;
    return;
  end
else
  if(sum(handles.individuals_behavior(:,handles.behaviorvalue)==-1)>0)
    uiwait(errordlg('some experiments have no individiuals for this classifier'));  drawnow;
    return;
  end
end

switch(handles.analysis)
  case 'feature_histogram'
    handles=feature_histogram_plot(handles);
  case 'feature_timeseries'
    handles=feature_timeseries_plot(handles);
  case 'behavior_barchart'
    handles=behavior_barchart_plot(handles);
  case 'behavior_timeseries'
    handles=behavior_timeseries_plot(handles);
  case 'bout_stats'
    handles=bout_stats_plot(handles);
  case 'interesting_feature_histograms'
    handles=interesting_feature_histograms_plot(handles);
end

guidata(hObject,handles);


% ---
function ret_val=get_label(feature_name,feature_units)

ret_val=[strrep(char(feature_name),'_','-') ' ('];
if(~isempty(feature_units.num))
  if(length(feature_units.num)>1)
    ret_val=[ret_val '['];
  end
  tmp=char(feature_units.num);
  tmp=[tmp repmat('-',size(tmp,1),1)];
  tmp=reshape(tmp',1,numel(tmp));
  tmp=tmp(1:end-1);
  ret_val=[ret_val tmp];
  if(length(feature_units.num)>1)
    ret_val=[ret_val ']'];
  end
  if(~isempty(feature_units.den))
    ret_val=[ret_val ' / ' ];
  end
end
if(~isempty(feature_units.den))
  if(length(feature_units.den)>1)
    ret_val=[ret_val '['];
  end
  tmp=char(feature_units.den);
  tmp=[tmp repmat('-',size(tmp,1),1)];
  tmp=reshape(tmp',1,numel(tmp));
  tmp=tmp(1:end-1);
  ret_val=[ret_val tmp];
  if(length(feature_units.den)>1)
    ret_val=[ret_val ']'];
  end
  if(isempty(feature_units.num))
    ret_val=[ret_val ' ^ -1 ' ];
  end
end
ret_val=[ret_val ')'];


% ---
function print_csv_help(fid,type,tstr,xstr,ystr)

tstr=strrep(tstr,'%','%%');
xstr=strrep(xstr,'%','%%');
ystr=strrep(ystr,'%','%%');

fprintf(fid,['%% type=' type '\n%% title=' tstr '\n%% xlabel=' xstr '\n%% ylabel=' ystr '\n\n']);


% ---
function print_csv_data(fid,data)

idx=1;
limit=2^14-1;
while idx<length(data)
  fprintf(fid,'%g, ',data(idx:min(idx+limit-1,end)));
  fprintf(fid,'\n');
  idx=idx+limit;
end


% ---
function h=plot_it(ha,xdata,ydata,style,centraltendency,dispersion,color,linewidth,linestyle,...
    fid,experimentlist)

if(~isnan(fid))
  fprintf(fid,'%% xdata\n');
  print_csv_data(fid,xdata);
end
if(style~=3)
  switch(centraltendency)
    case 1
      data_ct=nanmean(ydata,1);
      str_ct='mean';
    case 2
      data_ct=nanmedian(ydata,1);
      str_ct='median';
    case 3
      data_ct=mode(ydata,1);
      str_ct='mode';
  end
end
switch(style)
  case 1
    h=plot(ha,xdata,data_ct,'color',color,'linewidth',linewidth,'linestyle',linestyle);
    if(~isnan(fid))
      fprintf(fid,['%% ydata, ' str_ct '\n']);
      print_csv_data(fid,data_ct);
    end
  case 2
    switch(dispersion)
      case 1
        tmp=nanstd(ydata,[],1);
        data_dp=data_ct+tmp;
        data_dn=data_ct-tmp;
        str_dp=[str_ct ' + std dev'];
        str_dn=[str_ct ' - std dev'];
      case 2
        tmp=nanstd(ydata,[],1)./sqrt(sum(~isnan(ydata),1));
        data_dp=data_ct+tmp;
        data_dn=data_ct-tmp;
        str_dp=[str_ct ' + std err'];
        str_dn=[str_ct ' - std err'];
      case 3
        data_dp=prctile(ydata,95);
        data_dn=prctile(ydata,5);
        str_dp='95%%';
        str_dn='5%%';
      case 4
        data_dp=prctile(ydata,75);
        data_dn=prctile(ydata,25);
        str_dp='75%%';
        str_dn='25%%';
    end
    h=plot(ha,xdata,data_ct,'color',color,'linewidth',3*linewidth,'linestyle',linestyle);
    idx=isnan(data_dp) | isnan(data_dn);
    xdata=xdata(~idx);  data_dp=data_dp(~idx);  data_dn=data_dn(~idx);
    color2=(get(h,'color')+[4 4 4])/5;
    k=1;  m=0;  step=10000;
    hfront = get(ha,'Children');
    while(k<=length(xdata))
      idx=k:min(k+step,length(xdata));
      patch([xdata(idx) fliplr(xdata(idx))],[data_dp(idx) fliplr(data_dn(idx))],color2,'edgecolor','none','parent',ha);
      k=k+step+1;  m=m+1;
    end
    hback = setdiff(get(ha,'Children'),hfront);
    set(ha,'children',[hfront;hback]);
    %get(ha,'children');  set(ha,'children',circshift(ans,-m));  % send to back
    if(~isnan(fid))
      fprintf(fid,['%% ydata, ' str_dp '\n']);
      print_csv_data(fid,data_dp);
      fprintf(fid,['%% ydata, ' str_dn '\n']);
      print_csv_data(fid,data_dn);
      fprintf(fid,['%% ydata, ' str_ct '\n']);
      print_csv_data(fid,data_ct);
    end
  case 3
    h=plot(ha,xdata,ydata','color',color,'linewidth',linewidth,'linestyle',linestyle);
    h=h(1);
    for e=1:size(ydata,1)
      if(~isnan(fid))
        fprintf(fid,['%% ydata, experiment ' experimentlist{e} '\n']);
        print_csv_data(fid,ydata(e,:));
      end
    end
end
if(~isnan(fid)) fprintf(fid,'\n');  end


% --- 
function [during not_during]=calculate_feature_histogram(behavior_data,behavior_logic,behavior_data2,...
    feature_data,sexdata,individual,perwhat,behaviornot)

%if(isempty(behavior_data.allScores.scores))
if(isempty(feature_data.data))
  during=[];
  not_during=[];
  return;
end

if(iscell(feature_data.data{1}))
  vals=unique([feature_data.data{:}]);
  if(length(vals)>2)  error('uhoh');  end
  feature_data.data=cellfun(@(x) strcmp(x,vals{1}),feature_data.data,'uniformoutput',false);
end

if(isempty(behavior_data))
  idx=1:length(feature_data.data);
  if(~isnan(individual))  idx=individual;  end
  [cellfun(@(x,y) x(true & y(1:length(x))), feature_data.data(idx), sexdata(idx), 'uniformoutput',false)];
  during=[ans{:}];
  not_during=[];
  return;
end

during={};  not_during={};
for i=1:length(behavior_data.allScores.t0s)  % individual
  if((~isnan(individual)) && (i~=individual))  continue;  end

  tmp1 = compute_behavior_logic(behavior_data.allScores, i);
  tstart = behavior_data.allScores.tStart(i);
  tmp1 = tmp1(tstart : tstart+length(feature_data.data{i})-1);
  
  if(behavior_logic>1)
    tmp2 = compute_behavior_logic(behavior_data2.allScores, i);
    tstart = behavior_data2.allScores.tStart(i);
    tmp2 = tmp2(tstart : tstart+length(feature_data.data{i})-1);
  end

  if(behaviornot)  tmp1=~tmp1;  end

  switch(behavior_logic)
    case(1)
      partition_idx=tmp1;
    case(2)
      partition_idx=tmp1 & tmp2;
    case(3)
      partition_idx=tmp1 & ~tmp2;
    case(4)
      partition_idx=tmp1 | tmp2;
    case(5)
      partition_idx=tmp1 | ~tmp2;
  end
  if length(partition_idx)==(length(feature_data.data{i})+1)  % n-1 features
    partition_idx=partition_idx(1:end-1);
  end

  if(perwhat==1)  % per frame
    during{i}=feature_data.data{i}(partition_idx & sexdata{i}(1:length(partition_idx)));
    not_during{i}=feature_data.data{i}((~partition_idx) & sexdata{i}(1:length(partition_idx)));
  else  % per bout
    partition_idx=[0 partition_idx 0];
    start=find(~partition_idx(1:(end-1)) &  partition_idx(2:end));
    stop =find( partition_idx(1:(end-1)) & ~partition_idx(2:end));
    during{i}=[];  not_during{i}=[];
    if(length(start)>0)
      for j=1:length(start)
        if(sum(sexdata{i}(start(j):(stop(j)-1))) < ((stop(j)-start(j))/2))  continue;  end
        switch(perwhat)
          case 2
            during{i}(j)=mean(feature_data.data{i}(start(j):(stop(j)-1)));
          case 3
            during{i}(j)=median(feature_data.data{i}(start(j):(stop(j)-1)));
          case 4
            during{i}(j)=max(feature_data.data{i}(start(j):(stop(j)-1)));
          case 5
            during{i}(j)=min(feature_data.data{i}(start(j):(stop(j)-1)));
          case 6
            during{i}(j)=std(feature_data.data{i}(start(j):(stop(j)-1)));
        end
      end
      if(length(start)>1)
        for j=1:(length(start)-1)
          if(sum(sexdata{i}(stop(j):(start(j+1)-1))) < ((start(j+1)-stop(j))/2))  continue;  end
          switch(perwhat)
            case 2
              not_during{i}(j)=mean(feature_data.data{i}(stop(j):start(j+1)));
            case 3
              not_during{i}(j)=median(feature_data.data{i}(stop(j):start(j+1)));
            case 4
              not_during{i}(j)=max(feature_data.data{i}(stop(j):start(j+1)));
            case 5
              not_during{i}(j)=min(feature_data.data{i}(stop(j):start(j+1)));
            case 6
              not_during{i}(j)=std(feature_data.data{i}(stop(j):start(j+1)));
          end
        end
      end
    end
  end
end

during=[during{:}];
not_during=[not_during{:}];


% ---
function tmp=calculate_statistics2(table_data,behaviorlist,grouplist,comparison,fid,crit)

nrows=10+4*(comparison>0);
if(((size(table_data{1},1)==2)&&(comparison==0)) || ((size(table_data{1},1)==1)&&(comparison>0)))
  nrows=8+4*(comparison>0);
end

switch(comparison)
  case 0, ctxt='';
  case 1, ctxt=' allframes';
  case 2, ctxt=' notduring';
end

tmp={'behavior'};
for b=1:length(table_data)

  grouplistcopy=grouplist;
  %g=1;
  %while g<=length(table_data{b})
  %  if(isempty(table_data{b}{g}) || (sum(~isnan(table_data{b}{g}))==0))
  %    table_data{b}(g)=[];
  %    grouplistcopy(g)=[];
  %  else
  %    g=g+1;
  %  end
  %end

  tmp{nrows*b+1,1}=behaviorlist{b};
  for g=1:size(table_data{b},1)
    for c=1:(1+(comparison>0))
      idx=2:5;  if(comparison>0)  idx=c+(1:2:7);  end
      p=nan;
      if(sum(~isnan(table_data{b}{g,c}))>0)
        [~,p,~,~]=kstest(table_data{b}{g,c});
      end
      tmp{nrows*b+idx(1),g}=color_red(length(table_data{b}{g,c}),0);
      tmp{nrows*b+idx(2),g}=color_red(p,p<crit);
      tmp{nrows*b+idx(3),g}=color_red(nanmean(table_data{b}{g,c}),0);
      tmp{nrows*b+idx(4),g}=color_red(nanstd(table_data{b}{g,c}),0);
      ctxt2='';
      if(comparison>0)
        ctxt2=ctxt;  if(c==1)  ctxt2=' during';  end
      end
      if(b==1)
        tmp(idx,g)=repmat({[grouplistcopy{g} ctxt2]},4,1);
        if(g==1)
          idx=2:5;  if(comparison>0)  idx=[2:3; 4:5; 6:7; 8:9]';  end
          tmp(idx(:,1),size(table_data{b},1)+1)=repmat({'(sample size)'},size(idx,1),1);
          tmp(idx(:,2),size(table_data{b},1)+1)=repmat({'(K-S normal)'},size(idx,1),1);
          tmp(idx(:,3),size(table_data{b},1)+1)=repmat({'(mean)'},size(idx,1),1);
          tmp(idx(:,4),size(table_data{b},1)+1)=repmat({'(std. dev.)'},size(idx,1),1);
        end
      end
    end
  end

  if(((size(table_data{1},1)==2)&&(comparison==0)) || ((size(table_data{1},1)==1)&&(comparison>0)))
  
    if(b==1)  tmp{6+4*(comparison>0),1}='t-test';  end
    if(comparison>0)
      [~,p]=ttest2(table_data{b}{1,1},table_data{b}{1,2});
    else
      [~,p]=ttest2(table_data{b}{1,1},table_data{b}{2,1});
    end
    tmp{nrows*b+6+4*(comparison>0),1}=color_red(p,p<crit);

    if(b==1)  tmp{7+4*(comparison>0),1}='Wilcoxen';  end
    if(comparison>0)
      p=ranksum(table_data{b}{1,1},table_data{b}{1,2});
    else
      p=ranksum(table_data{b}{1,1},table_data{b}{2,1});
    end
    tmp{nrows*b+7+4*(comparison>0),1}=color_red(p,p<crit);

  else

    f1=cellfun(@(x,y) y*ones(1,length(x)),table_data{b},...
         num2cell(repmat((1:size(table_data{b},1))',1,size(table_data{b},2))),...
        'uniformoutput',false);
    f2=cellfun(@(x,y) y*ones(1,length(x)),table_data{b},...
         num2cell(repmat((1:size(table_data{b},2)) ,size(table_data{b},1),1)),...
        'uniformoutput',false);
    if(comparison>0)
      [p,table,stats]=anovan([table_data{b}{:}],{[f1{:}] [f2{:}]},'model','interaction','display','off');
    else
      [p,table,stats]=anova1([table_data{b}{:}],[f1{:}],'off');
    end
    c=multcompare(stats,'alpha',crit,'display','off');
    for(pp=1:length(p))
      tmp{nrows*b+6+4*(comparison>0),pp}=color_red(p(pp),p(pp)<crit);
    end
    if(b==1)
      if(comparison>0)
        tmp(10,1:4)={'group','comparison','interaction','(ANOVA)'};
      else
        tmp(6,1)={'(ANOVA)'};
      end
    end
    for g2=1:size(c,1)
      if(b==1)
        tmp((7:9)+4*(comparison>0),g2)=repmat({strcat(grouplistcopy{c(g2,1)},'-',grouplistcopy{c(g2,2)})},3,1);
      end
      foo=(sign(c(g2,3))==sign(c(g2,5)));
      tmp{nrows*b+7+4*(comparison>0),g2}=color_red(c(g2,3),foo);
      tmp{nrows*b+8+4*(comparison>0),g2}=color_red(c(g2,4),foo);
      tmp{nrows*b+9+4*(comparison>0),g2}=color_red(c(g2,5),foo);
    end
    pooh=size(c,1);
    if(comparison>0)
      [p,table,stats]=anovan([table_data{b}{:}],{[f2{:}] [f1{:}]},'model','interaction','display','off');
      c=multcompare(stats,'alpha',crit,'display','off');
      if(b==1)
        tmp((7:9)+4*(comparison>0),pooh+1)=repmat({strcat('during-',ctxt)},3,1);
      end
      foo=(sign(c(1,3))==sign(c(1,5)));
      tmp{nrows*b+7+4*(comparison>0),pooh+1}=color_red(c(1,3),foo);
      tmp{nrows*b+8+4*(comparison>0),pooh+1}=color_red(c(1,4),foo);
      tmp{nrows*b+9+4*(comparison>0),pooh+1}=color_red(c(1,5),foo);
    end
    if(b==1)
      tmp((7:9)+4*(comparison>0),pooh+1+(comparison>0))={'(5%, Tukey post-hoc)'; '(mean)'; '(95%)'};
    end
  end
end
tmp{end+1,1}='';

if(~isnan(fid))
  fprintf(fid,'%% statistics\n');
  for j=1:(size(tmp,1)/nrows)
    for i=1:nrows
      if(j==1)  fprintf(fid,'%% ');  end
      tmp(nrows*(j-1)+i,:);
      ans(cellfun(@isempty,ans))='';
      cellfun(@(x) regexprep(x,'<[^<]*>','*'),ans,'uniformoutput',false);
      fprintf(fid,'%s, ',ans{:});
      fprintf(fid,'\n');
    end
  end
end


% ---
function [behavior_data, behavior_data2, behavior_data3, feature_data, sex_data]=...
    cull_short_trajectories(handles, behavior_data, behavior_data2, behavior_data3, feature_data, sex_data)

if(~isempty(behavior_data))
  idx=find((behavior_data.allScores.tEnd-behavior_data.allScores.tStart+1)./handles.fps < ...
      handles.minimumtrajectorylength);
else
  idx=find((cellfun(@length,feature_data.data)./handles.fps) < handles.minimumtrajectorylength);
end

if(~isempty(behavior_data))
  behavior_data.allScores.scores(idx)=[];
  behavior_data.allScores.tStart(idx)=[];
  behavior_data.allScores.tEnd(idx)=[];
  behavior_data.allScores.t0s(idx)=[];
  behavior_data.allScores.t1s(idx)=[];
end
if(~isempty(behavior_data2))
  behavior_data2.allScores.scores(idx)=[];
  behavior_data2.allScores.tStart(idx)=[];
  behavior_data2.allScores.tEnd(idx)=[];
  behavior_data2.allScores.t0s(idx)=[];
  behavior_data2.allScores.t1s(idx)=[];
end
if(~isempty(behavior_data3))
  behavior_data3.allScores.scores(idx)=[];
  behavior_data3.allScores.tStart(idx)=[];
  behavior_data3.allScores.tEnd(idx)=[];
  behavior_data3.allScores.t0s(idx)=[];
  behavior_data3.allScores.t1s(idx)=[];
end
if(~isempty(feature_data))
  feature_data.data(idx)=[];
end
if(~isempty(sex_data))
  sex_data(idx)=[];
end


% ---
function handles=feature_histogram_plot(handles)

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

handlesexperimentlist=[handles.experimentlist{:}];

cumsum_num_indi_per_exp=[0 cumsum(handles.individuals_feature)];
cumsum_num_exp_per_group=[0 cumsum(cellfun(@length,handles.experimentlist))];
mat2cell(cumsum_num_exp_per_group(1:end-1),1,ones(1,length(cumsum_num_exp_per_group)-1));
cellfun(@(x,y) x+y,handles.experimentvalue,ans,'uniformoutput',false);
selected_exp=[ans{:}];
cumsum_num_selexp_per_group=[0 cumsum(cellfun(@length,handles.experimentvalue))];

%ggee=1:length(handlesexperimentlist);
ggee=selected_exp;
individual=handles.individualidx;
if(isnumeric(individual))
  ggee=find(cumsum_num_indi_per_exp<individual,1,'last');
  individual=individual-cumsum_num_indi_per_exp(ggee);
end

fid=nan;  if(handles.dump2csv)  fid=fopen('most_recent_figure.csv','w');  end

handles.type='feature histogram';

hf=figure('toolbar','figure');

bb=handles.behaviorvalue;
if (bb==(length(handles.behaviorlist)+1))
  bb=1:(bb-1);
elseif((bb==(length(handles.behaviorlist)+2)) || (strcmp(get(handles.BehaviorList,'enable'),'off')))
  bb=0;
end

behavior_logic=handles.behaviorlogic;
score_file2=[];
if((length(bb)>1) || (bb>0))  score_file2=handles.scorefiles{handles.behaviorvalue2};  end
feature_value=handles.featurevalue;
feature_list=handles.featurelist;
comparison=handles.comparison;
if((length(bb)==1) && (bb==0))  comparison=0;  end
nbins=handles.nbins;
style=handles.featurehistogram_style;
centraltendency=handles.centraltendency;
dispersion=handles.dispersion;
behaviornot=handles.behaviornot;

h=[];
table_data={};
for b=bb

  tstr='all frames';
  if(b>0)
    tstr='';  if(behaviornot)  tstr='NOT ';  end
    tstr=[tstr char(strrep(handles.behaviorlist(b),'_','-'))];
  end
  switch(handles.behaviorlogic)
    case 2
      tstr=[tstr ' AND '];
    case 3
      tstr=[tstr ' AND NOT '];
    case 4
      tstr=[tstr ' OR '];
    case 5
      tstr=[tstr ' OR NOT '];
  end
  if(handles.behaviorlogic>1)
    tstr=[tstr char(strrep(handles.behaviorlist(handles.behaviorvalue2),'_','-'))];
  end
  units=load(fullfile(handlesexperimentlist{ggee(1)},handles.perframe_dir,...
      [feature_list{feature_value} '.mat']),'units');
  if(style<4)
    xstr=get_label(feature_list(feature_value),units.units);
    ystr='normalized';
  else
    xstr='group';
    ystr=get_label(feature_list(feature_value),units.units);
  end

  if(handles.dump2csv)  print_csv_help(fid,handles.type,tstr,xstr,ystr);  end

  if(length(bb)>1)
    ceil(sqrt(length(bb)));
    ha=subplot(ceil(length(bb)/ans),ans,b,'parent',hf);
  else
    ha=subplot(1,1,1,'parent',hf);
  end
  hold(ha,'on');

  score_file=[];
  if(b>0)  score_file=handles.scorefiles{b};  end

  num_indi=0;
  during_data=cell(length(ggee),length(individual));
  not_during_data=cell(length(ggee),length(individual));
  parfor gei=1:numel(ggee)
  %for gei=1:numel(ggee)
    ge = ggee(gei);

    if(b>0)
      behavior_data=load(fullfile(handlesexperimentlist{ge},score_file));
      behavior_data=update_t01s_from_postprocessed(behavior_data);
      if(behavior_logic>1)
        behavior_data2=load(fullfile(handlesexperimentlist{ge},score_file2));
        behavior_data2=update_t01s_from_postprocessed(behavior_data2);
      else
        behavior_data2=[];
      end
    else
      behavior_data=[];
      behavior_data2=[];
    end
    feature_data=load(fullfile(handlesexperimentlist{ge},handles.perframe_dir,...
          [feature_list{feature_value} '.mat']));

    [behavior_data,behavior_data2,~,feature_data,sex_data]=...
        cull_short_trajectories(handles,behavior_data,behavior_data2,[],feature_data,handles.sexdata{ge});
    num_indi=num_indi+length(feature_data.data);

    ii=0;
    parfor_tmp=cell(1,length(individual));
    not_parfor_tmp=cell(1,length(individual));
    for i = individual
      ii=ii+1;
      if(iscell(i))  i=char(i);  end
      tmp2=[];
      switch(i)
        case('M')
          tmp2=sex_data;
        case('F')
          tmp2=cellfun(@not,sex_data,'uniformoutput',false);
        otherwise
          tmp2=cellfun(@(x) ones(1,length(x)),sex_data,'uniformoutput',false);
      end
      tmploop=nan;  if isnumeric(i)  tmploop=i;  end

      [parfor_tmp{ii} not_parfor_tmp{ii}]=calculate_feature_histogram(...
          behavior_data,behavior_logic,behavior_data2,feature_data,tmp2,tmploop,...
          handles.featurehistogram_style2,handles.behaviornot);

      if(comparison==1)
        not_parfor_tmp{ii}=[parfor_tmp{ii} not_parfor_tmp{ii}];
      end
    end
    during_data(gei,:)=parfor_tmp;
    not_during_data(gei,:)=not_parfor_tmp;
  end
  during_data=reshape(during_data,1,prod(size(during_data)));
  not_during_data=reshape(not_during_data,1,prod(size(not_during_data)));

  if(num_indi==0)
    delete(hf);
    set(handles.Status,'string','Ready.','foregroundcolor','g');
    set(handles.figure1,'pointer','arrow');
    uiwait(errordlg('no valid data.  check minimum trajectory length.'));  drawnow;
    return;
  end

  max(cellfun(@(x) size(x,2),during_data));
  cellfun(@(x) [x nan(size(x,1),ans-size(x,2))],during_data,'uniformoutput',false);
  during_data=cat(1,ans{:});
  max(cellfun(@(x) size(x,2),not_during_data));
  cellfun(@(x) [x nan(size(x,1),ans-size(x,2))],not_during_data,'uniformoutput',false);
  not_during_data=cat(1,ans{:});

  tmp=[];
  if(~isempty(during_data))
    tmp=reshape(during_data,1,numel(during_data));
  end
  if((comparison>0) && ~isempty(not_during_data))
    tmp=[tmp reshape(not_during_data,1,numel(not_during_data))];
  end

  %low=[];  high=[];  nearzero=[];
  %if(~isempty(during_data))
  %  low=min(min(during_data));
  %  high=max(max(during_data));
  %  unique(reshape(abs(during_data),1,prod(size(during_data))));
  %  nearzero=ans(1);  if(ans(1)==0)  nearzero=ans(2);  end
  %end
  %if((comparison>0) && ~isempty(not_during_data))
  %  low=min([low min(not_during_data)]);
  %  high=max([high max(not_during_data)]);
  %  unique(reshape(abs(not_during_data),1,prod(size(not_during_data))));
  %  tmp=ans(1);  if(ans(1)==0)  tmp=ans(2);  end
  %  nearzero=min(tmp,nearzero);
  %end

  foo=prctile(tmp,[1 99]);
  low=foo(1);  high=foo(2);
  unique(abs(tmp));
  nearzero=ans(1);  if(ans(1)==0)  nearzero=ans(2);  end
  if(~isempty(low) && ~isempty(high) && ~isempty(nearzero))
    if(handles.logbinsize)
      if((low>=0) && (high>0))
        bins=logspace(log10(max(low,nearzero)),log10(high),nbins);
      elseif((low<0) && (high<=0))
        bins=fliplr(-logspace(log10(max(abs(high),nearzero)),log10(abs(low)),nbins));
      elseif((low<0) && (high>0))
        bins=[fliplr(-logspace(log10(nearzero),log10(abs(low)),nbins)) ...
            logspace(log10(nearzero),log10(abs(high)),nbins)];
      end
      plot(ha,bins,zeros(1,length(bins)),'k.');
    else
      bins=linspace(low,high,nbins);
    end
  end

  ii=0;
  for i = individual
    ii=ii+1;
    if(iscell(i))  i=char(i);  end
    for g=1:length(handles.grouplist)
      color=handles.colors(g,:);

      if(handles.dump2csv)  fprintf(fid,['%% group ' handles.grouplist{g} '\n']);  end

      if ~isnumeric(i)
        idx=(cumsum_num_selexp_per_group(g)+1):(cumsum_num_selexp_per_group(g+1));
      else
        find(cumsum_num_exp_per_group<ggee,1,'last');
        if(ans~=g)  continue;  end
        idx=1;
      end
      idx2=idx+(ii-1)*numel(ggee);
      linestyle='-';  if(ii>1)  linestyle='--';  end

      if(style<4)
        if(comparison>0)
          if(handles.dump2csv)
            if(comparison==1)
              fprintf(fid,['%% all frames\n']);
            else
              fprintf(fid,['%% not during\n']);
            end
          end
          if(~isempty(not_during_data))
            hist_not_during=histc(not_during_data(idx2,:)',bins);
            if(size(hist_not_during,1)==1)  hist_not_during=hist_not_during';  end
            hist_not_during.*repmat(([0 diff(bins)]+[diff(bins) 0])'/2,1,size(hist_not_during,2));
            hist_not_during=hist_not_during./repmat(sum(ans,1),size(hist_not_during,1),1);
            plot_it(ha,bins,hist_not_during',style,centraltendency,dispersion,color,1,linestyle,...
                fid,handlesexperimentlist(idx));
          end
        end
        if(handles.dump2csv)
          if(comparison==1)
            fprintf(fid,['%% all frames\n']);
          else
            fprintf(fid,['%% not during\n']);
          end
        end
        if(~isempty(during_data))
          hist_during=histc(during_data(idx2,:)',bins);
          if(size(hist_during,1)==1)  hist_during=hist_during';  end
          hist_during.*repmat(([0 diff(bins)]+[diff(bins) 0])'/2,1,size(hist_during,2));
          hist_during=hist_during./repmat(sum(ans,1),size(hist_during,1),1);
          linewidth=1;  if(comparison>0)  linewidth=2;  end
          plot_it(ha,bins,hist_during',style,centraltendency,dispersion,color,linewidth,linestyle,...
              fid,handlesexperimentlist(idx));
          if(ii==1)  h(g)=ans;  end
        end
      else
        if(comparison>0)
          xticklabels{g+(ii-1)*length(handles.grouplist)}=handles.grouplist{g};
          tmp=nanmean(not_during_data(idx2,:),2);
          if(handles.dump2csv)
            if(comparison==1)
              fprintf(fid,['%% all frames\n']);
            else
              fprintf(fid,['%% not during\n']);
            end
            fprintf(fid,['%% percentiles 1, 5, 25, 50, 75, 95, 99\n']);
            fprintf(fid,'%g, ',prctile(tmp,[1 5 25 50 75 95 99]));
            fprintf(fid,'\n');
          end
          tmp=boxplot(ha,tmp,'positions',g+(ii-1)*length(handles.grouplist)+0.5,...
              'widths',0.5/((comparison>0)+1),'colors',color);
  %        if(ii==1)  h{g}=tmp;  end
          findobj(tmp,'type','line');
          if(ii==1)
            set(ans,'linestyle','-');
          else
            set(ans,'linestyle','--');
          end
          findobj(tmp,'tag','Outliers');
          set(ans,'markeredgecolor',color);
        end
        xticklabels{g+(ii-1)*length(handles.grouplist)}=handles.grouplist{g};
        tmp=nanmean(during_data(idx2,:),2);
        if(handles.dump2csv)
          if(comparison==1)
            fprintf(fid,['%% all frames\n']);
          else
            fprintf(fid,['%% not during\n']);
          end
          fprintf(fid,['%% percentiles 1, 5, 25, 50, 75, 95, 99\n']);
          fprintf(fid,'%g, ',prctile(tmp,[1 5 25 50 75 95 99]));
          fprintf(fid,'\n');
        end
        tmp=boxplot(ha,tmp,'positions',g+(ii-1)*length(handles.grouplist),...
            'widths',0.5/((comparison>0)+1),'colors',color);
%        if(ii==1)  h{g}=tmp;  end
        findobj(tmp,'type','line');
        if(ii==1)
          set(ans,'linestyle','-');
        else
          set(ans,'linestyle','--');
        end
        if(comparison>0)
          set(ans,'linewidth',2);
        end
        findobj(tmp,'tag','Outliers');
        set(ans,'markeredgecolor',color);
      end
    end

    table_data{end+1}={};
    if(handles.dump2csv)  fprintf(fid,'\n%% raw data\n');  end
    for g=1:length(handles.grouplist)
      if ischar(i)
        idx=(cumsum_num_selexp_per_group(g)+1):(cumsum_num_selexp_per_group(g+1));
      else
        find(cumsum_num_exp_per_group<ggee,1,'last');
        if(ans~=g)  continue;  end
        idx=1;
      end
      if(handles.dump2csv)
        fprintf(fid,['%% group ' handles.grouplist{g} '\n']);
        for e=1:length(idx)
          fprintf(fid,'%% experiment %s\n',handlesexperimentlist{selected_exp(idx(e))});
          if(comparison>0)
            if(comparison==1)
              fprintf(fid,'%% all frames\n');
            else
              fprintf(fid,'%% not during\n');
            end
            if(~isempty(not_during_data))
              print_csv_data(fid,not_during_data(idx(e),:));
            end
          end
          if(comparison==1)
            fprintf(fid,'%% all frames\n');
          else
            fprintf(fid,'%% during\n');
          end
          if(~isempty(during_data))
            print_csv_data(fid,during_data(idx(e),:));
          end
          fprintf(fid,'\n');
        end
      end

      if(~isempty(during_data))  table_data{end}{g,1}=nanmean(during_data(idx,:),2)';  end
      if(comparison>0)
        if(~isempty(not_during_data))  table_data{end}{g,2}=nanmean(not_during_data(idx,:),2)';  end
      end
    end
  end
  
  if(handles.dump2mat)
    find(b==bb);
    raw_during_data{ans}=during_data;
    raw_not_during_data{ans}=not_during_data;
  end

  title(ha,tstr,'interpreter','none');
  xlabel(ha,xstr,'interpreter','none');
  ylabel(ha,ystr,'interpreter','none');
  if(style<4)
    axis(ha,'tight');
  else
    axisalmosttight([],ha);
    set(ha,'xtick',(1:length(xticklabels))+0.25*(comparison>0),'xticklabel',xticklabels);
  end
  zoom(ha,'reset');
end

if(handles.dump2mat)
  save('most_recent_figure.mat','handles','raw_during_data','raw_not_during_data');
end

h2=[];  hh2={};
if(~isnumeric(individual) && (style<4))
  idx=find(h>0);
  h2=[h2 h(idx)];
  hh2=[hh2 handles.grouplist];
  %legend(ha,h(idx),[cellfun(@(x) [x ' ' handles.individuallist{handles.individualvalue}],...
  %    handles.grouplist,'uniformoutput',false)],'interpreter','none');
%else
%  legend(ha,h(idx),handles.individuallist(handles.individualvalue),'interpreter','none');
end
if(iscell(individual))
  h2=[h2 plot(0,0,'k-')];   hh2={hh2{:} 'males'};
  h2=[h2 plot(0,0,'k--')];  hh2={hh2{:} 'females'};
  %set(h((end-1):end),'visible','off');
end
if(comparison>0)
  h2=[h2 plot(0,0,'k-','linewidth',2)];  hh2={hh2{:} 'during'};
  if(comparison==1)
    h2=[h2 plot(0,0,'k-','linewidth',1)];  hh2={hh2{:} 'all frames'};
  else
    h2=[h2 plot(0,0,'k-','linewidth',1)];  hh2={hh2{:} 'not during'};
  end
  %set(h2,'visible','off');
  %legend(ha,h2,hh2,'interpreter','none');
end
if(~isempty(h2))
  legend(ha,h2,hh2,'interpreter','none');
end

uicontrol(hf,'style','pushbutton','string','Params','position',[5 5 60 20],...
    'callback',@figure_params_callback);
if(ischar(individual) && ((length(handles.grouplist)>1) || (comparison>0)) && ...
    (~isempty(during_data)) && (~isempty(not_during_data)))
  uicontrol(hf,'style','pushbutton','string','Stats','position',[70 5 50 20],...
      'callback',@figure_stats_callback);
  if((length(bb)==1) && (bb==0))
    {'all frames'};
  else
    handles.behaviorlist(bb);
  end
  handles.statistics=calculate_statistics2(table_data,ans,handles.grouplist,...
      comparison,fid,handles.pvalue);
end

if(handles.dump2csv)  fclose(fid);  end

guidata(hf,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% ---
function progress_bar(~,~,h)

fid=fopen(fullfile(tempdir,'progressbar.txt'),'r');
den=fread(fid,1,'uint32');
fseek(fid,0,1);
num=ftell(fid)/4;
waitbar(num/den,h);
fclose(fid);
drawnow;


% ---
function handles=interesting_feature_histograms_plot(handles)

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

cumsum_num_exp_per_group=[0 cumsum(cellfun(@length,handles.experimentlist))];
mat2cell(cumsum_num_exp_per_group(1:end-1),1,ones(1,length(cumsum_num_exp_per_group)-1));
cellfun(@(x,y) x+y,handles.experimentvalue,ans,'uniformoutput',false);
selected_exp=[ans{:}];
cumsum_num_selexp_per_group=[0 cumsum(cellfun(@length,handles.experimentvalue))];

handlesexperimentlist=[handles.experimentlist{:}];
handlesexperimentlist=handlesexperimentlist(selected_exp);
handlessexdata=handles.sexdata(selected_exp);

nexperiments=length(handlesexperimentlist);
nbehaviors=length(handles.behaviorlist);
nfeatures=length(handles.featurelist);

if(isempty(handles.interestingfeaturehistograms_cache))
  table_data={};

  h=waitbar(0,'This will likely take awhile...',...
      'CreateCancelBtn','fid=fopen(fullfile(tempdir,''cancel.txt''),''w''); fclose(fid);');
  fid=fopen(fullfile(tempdir,'progressbar.txt'),'w');
  fwrite(fid,nfeatures*sum(cellfun(@length,handles.experimentvalue)),'uint32');
  fclose(fid);
  t = timer('TimerFcn',{@progress_bar,h}, 'Period', 3, 'ExecutionMode', 'fixedRate');
  start(t);

  num_indi=0;
  table_data=zeros(nexperiments,max(1,nbehaviors),nfeatures,9);
  parfor ge=1:nexperiments
%   for ge=1:nexperiments
    bdata={};  behavior_data=[];
    for b=1:nbehaviors
      behavior_data=load(fullfile(handlesexperimentlist{ge},handles.scorefiles{b}));
      behavior_data=update_t01s_from_postprocessed(behavior_data);
      [bdata{b},~,~,~,~]=cull_short_trajectories(handles,behavior_data,[],[],[],[]);
      num_indi=num_indi+length(bdata{b}.allScores.scores);
    end
    [~,~,~,~,sex_data]=cull_short_trajectories(handles,behavior_data,[],[],[],handlessexdata{ge});
    sexdata=[];
    switch(handles.individualidx)
      case('A')
        sexdata=cellfun(@(x) ones(1,length(x)),sex_data,'uniformoutput',false);
      case('M')
        sexdata=sex_data;
      case('F')
        sexdata=cellfun(@not,sex_data,'uniformoutput',false);
    end

    bad{ge}={};
    parfor_tmp=zeros(nbehaviors,nfeatures,9);
    for f=1:nfeatures
      if(exist(fullfile(tempdir,'cancel.txt')))  break;  end
      feature_data=load(fullfile(handlesexperimentlist{ge},handles.perframe_dir,...
          [handles.featurelist{f} '.mat']));
      [~,~,~,fdata,~]=cull_short_trajectories(handles,[],[],[],feature_data,[]);

      %sexdata={};
      %for s=1:length(fdata.data)
      %  sexdata{s}=ones(1,length(fdata.data{s}));
      %end
      for b=1:nbehaviors
        if(exist(fullfile(tempdir,'cancel.txt')))  break;  end

        [during not_during]=calculate_feature_histogram(bdata{b},1,[],...
            fdata,sexdata,nan,handles.featurehistogram_style2,0);
        parfor_tmp(b,f,:)=[mean(during) mean(not_during) mean([during not_during]) ...
            std(during) std(not_during) std([during not_during]) ...
            length(during) length(not_during) length([during not_during])];
      end
      if(nbehaviors==0)
        if(exist(fullfile(tempdir,'cancel.txt')))  break;  end

        [during not_during]=calculate_feature_histogram([],1,[],...
            fdata,sexdata,nan,handles.featurehistogram_style2,0);
        parfor_tmp(1,f,:)=[mean(during) mean(not_during) mean([during not_during]) ...
            std(during) std(not_during) std([during not_during]) ...
            length(during) length(not_during) length([during not_during])];
      end
      fid=fopen(fullfile(tempdir,'progressbar.txt'),'a');  fwrite(fid,1,'uint32');  fclose(fid);
    end
    table_data(ge,:,:,:)=parfor_tmp;
  end

  stop(t);
  delete(t);
  delete(h);
  delete(fullfile(tempdir,'progressbar.txt'));

  if(exist(fullfile(tempdir,'cancel.txt')))
    delete(fullfile(tempdir,'cancel.txt'));
    set(handles.Status,'string','Ready.','foregroundcolor','g');
    set(handles.figure1,'pointer','arrow');
    drawnow;
    return;
  end

  if (num_indi==0) && (nbehaviors>0)
    handles.interestingfeaturehistograms_cache=nan;
  else
    handles.interestingfeaturehistograms_cache=table_data;
  end
else
  table_data=handles.interestingfeaturehistograms_cache;
end

if(isnan(handles.interestingfeaturehistograms_cache))
  set(handles.Status,'string','Ready.','foregroundcolor','g');
  set(handles.figure1,'pointer','arrow');
  uiwait(errordlg('no valid data.  check minimum trajectory length.'));  drawnow;
  return;
end

if(handles.comparison2==0)
  tmp2=[];
  for g=1:length(handles.grouplist)
    gg=cumsum_num_selexp_per_group(g)+(1:length(handles.experimentlist{g}));
    if(nbehaviors>0)
      tmp2=[tmp2; ...
          repmat(g,nbehaviors*nfeatures,1) ...
          reshape(squeeze(sum(table_data(gg,:,:,7))),nbehaviors*nfeatures,1) ...
          repmat(-2,nbehaviors*nfeatures,1) ...
          reshape(squeeze(sum(table_data(gg,:,:,8))),nbehaviors*nfeatures,1) ...
          reshape(repmat((1:nbehaviors)',1,nfeatures),nbehaviors*nfeatures,1) ...
          reshape(repmat(1:nfeatures,nbehaviors,1),nbehaviors*nfeatures,1) ...
          reshape(squeeze((mean(table_data(gg,:,:,1))-mean(table_data(gg,:,:,2)))./ ...
            sqrt((std(table_data(gg,:,:,1)).^2+std(table_data(gg,:,:,2)).^2)/2)), ...
            nbehaviors*nfeatures,1) ...
          nan(nbehaviors*nfeatures,1)];
%             sqrt((mean(table_data(gg,:,:,4).^2)+mean(table_data(gg,:,:,5).^2))/2)), ...
      tmp2=[tmp2; ...
          repmat(g,nbehaviors*nfeatures,1) ...
          reshape(squeeze(sum(table_data(gg,:,:,7))),nbehaviors*nfeatures,1) ...
          repmat(-1,nbehaviors*nfeatures,1) ...
          reshape(squeeze(sum(table_data(gg,:,:,9))),nbehaviors*nfeatures,1) ...
          reshape(repmat((1:nbehaviors)',1,nfeatures),nbehaviors*nfeatures,1) ...
          reshape(repmat(1:nfeatures,nbehaviors,1),nbehaviors*nfeatures,1) ...
          reshape(squeeze((mean(table_data(gg,:,:,1))-mean(table_data(gg,:,:,3)))./ ...
            sqrt((std(table_data(gg,:,:,1)).^2+std(table_data(gg,:,:,3)).^2)/2)), ...
            nbehaviors*nfeatures,1) ...
          nan(nbehaviors*nfeatures,1)];
%             sqrt((mean(table_data(gg,:,:,4).^2)+mean(table_data(gg,:,:,6).^2))/2)), ...
    end

    if(g==length(handles.grouplist))  break;  end
    for g2=(g+1):length(handles.grouplist)
      gg2=cumsum_num_selexp_per_group(g2)+(1:length(handles.experimentlist{g2}));
      if(nbehaviors>0)
        tmp2=[tmp2; ...
            repmat(g,nbehaviors*nfeatures,1) ...
            reshape(squeeze(sum(table_data(gg,:,:,7))),nbehaviors*nfeatures,1) ...
            repmat(g2,nbehaviors*nfeatures,1) ...
            reshape(squeeze(sum(table_data(gg2,:,:,7))),nbehaviors*nfeatures,1) ...
            reshape(repmat((1:nbehaviors)',1,nfeatures),nbehaviors*nfeatures,1) ...
            reshape(repmat(1:nfeatures,nbehaviors,1),nbehaviors*nfeatures,1) ...
            reshape(squeeze((mean(table_data(gg,:,:,1))-mean(table_data(gg2,:,:,1)))./ ...
              sqrt((std(table_data(gg,:,:,1)).^2+std(table_data(gg2,:,:,1)).^2)/2)), ...
              nbehaviors*nfeatures,1) ...
            reshape(squeeze((mean(table_data(gg,:,:,3))-mean(table_data(gg2,:,:,3)))./ ...
              sqrt((std(table_data(gg,:,:,3)).^2+std(table_data(gg2,:,:,3)).^2)/2)), ...
              nbehaviors*nfeatures,1)];
%               sqrt((mean(table_data(gg,:,:,4).^2)+mean(table_data(gg2,:,:,4).^2))/2)), ...
%               sqrt((mean(table_data(gg,:,:,6).^2)+mean(table_data(gg2,:,:,6).^2))/2)), ...
        tmp2=[tmp2; ...
            repmat(g,nbehaviors*nfeatures,1) ...
            reshape(squeeze(sum(table_data(gg,:,:,8))),nbehaviors*nfeatures,1) ...
            repmat(g2,nbehaviors*nfeatures,1) ...
            reshape(squeeze(sum(table_data(gg2,:,:,8))),nbehaviors*nfeatures,1) ...
            reshape(repmat(-(1:nbehaviors)',1,nfeatures),nbehaviors*nfeatures,1) ...
            reshape(repmat(1:nfeatures,nbehaviors,1),nbehaviors*nfeatures,1) ...
            reshape(squeeze((mean(table_data(gg,:,:,2))-mean(table_data(gg2,:,:,2)))./ ...
              sqrt((std(table_data(gg,:,:,2)).^2+std(table_data(gg2,:,:,2)).^2)/2)), ...
              nbehaviors*nfeatures,1) ...
            reshape(squeeze((mean(table_data(gg,:,:,3))-mean(table_data(gg2,:,:,3)))./ ...
              sqrt((mean(table_data(gg,:,:,6)).^2+mean(table_data(gg2,:,:,6)).^2)/2)), ...
              nbehaviors*nfeatures,1)];
%               sqrt((mean(table_data(gg,:,:,5).^2)+mean(table_data(gg2,:,:,5).^2))/2)), ...
      end
      tmp2=[tmp2; ...
          repmat(g,nfeatures,1) ...
          squeeze(sum(table_data(gg,1,:,9))) ...
          repmat(g2,nfeatures,1) ...
          squeeze(sum(table_data(gg2,1,:,9))) ...
          zeros(nfeatures,1) ...
          (1:nfeatures)' ...
          squeeze((mean(table_data(gg,1,:,3))-mean(table_data(gg2,1,:,3)))./ ...
            sqrt((std(table_data(gg,1,:,3)).^2+std(table_data(gg2,1,:,3)).^2)/2)), ...
          nan(nfeatures,1)];
%             sqrt((mean(table_data(gg,1,:,6).^2)+mean(table_data(gg2,1,:,6).^2))/2)), ...
    end
  end

  if(handles.omitnan)
    idx=find(~isnan(tmp2(:,7)));
    tmp2=tmp2(idx,:);
  end
  if(handles.omitinf)
    idx=find(~isinf(tmp2(:,7)));
    tmp2=tmp2(idx,:);
  end
  if(handles.absdprimezscore)
    tmp2(:,7)=abs(tmp2(:,7));
    tmp2(:,8)=abs(tmp2(:,8));
  end
  tmp2=sortrows(tmp2,-7);

  if(isempty(tmp2))
    set(handles.Status,'string','Ready.','foregroundcolor','g');
    set(handles.figure1,'pointer','arrow');
    uiwait(errordlg('no valid data.  check individual sex.'));  drawnow;
    return;
  end

  tmp=cell(size(tmp2,1),8);
  tmp(:,1)=handles.grouplist(tmp2(:,1));
  tmp(:,2)=cellstr(num2str(tmp2(:,2),'%-d'));
  idx=(tmp2(:,3)>0);    tmp(idx,3)=handles.grouplist(tmp2(idx,3));
  idx=(tmp2(:,3)==-1);  tmp(idx,3)=cellstr('all frames');
  idx=(tmp2(:,3)==-2);  tmp(idx,3)=cellstr('not during');
  tmp(:,4)=cellstr(num2str(tmp2(:,4),'%-d'));
  idx=(tmp2(:,5)>0);   tmp(idx,5)=handles.behaviorlist(tmp2(idx,5));
  idx=(tmp2(:,5)==0);  tmp(idx,5)=cellstr('all frames');
  idx=(tmp2(:,5)<0);   tmp(idx,5)=cellstr([repmat('not ',sum(idx),1) char(handles.behaviorlist(-tmp2(idx,5)))]);
  tmp(:,6)=handles.featurelist(tmp2(:,6));
  tmp(:,7)=cellstr(num2str(tmp2(:,7)));
  idx=find(~isnan(tmp2(:,8)));
  tmp(idx,8)=cellstr(num2str(tmp2(idx,8)));

  handles2.figure1=handles.figure1;
  handles2.figure2=figure('menubar','none','toolbar','none','numbertitle','off',...
      'name','interesting feature histograms');
  handles2.Table=uitable('Data',tmp,...
      'ColumnName',{'Group1' 'n1' 'Group2' 'n2' 'Behavior' 'Feature' 'd''' 'd'' all frames'},...
      'ColumnWidth',num2cell(7*max(cellfun(@length,tmp))),...
      'RowName',[],'RowStriping','on','BackgroundColor',[1 1 1; 0.95 0.95 0.95],...
      'CellSelectionCallback',@CellSelectionCallback);
  extT=get(handles2.Table,'extent');
  posF=get(handles2.figure2,'position');
  %set(hf,'position',[posF(1) posF(2) 16*(posF(4)<extT(4))+extT(3) posF(4)]);
  set(handles2.figure2,'position',[posF(1) posF(2) min(posF(3),16*(posF(4)<extT(4))+extT(3)) min(posF(4),extT(4))]);
  %set(ht,'position',[0 0 16*(posF(4)<extT(4))+extT(3) posF(4)]);
  set(handles2.Table,'units','normalized','position',[0 0 1 1]);
  handles2.table_data=tmp2;
  handles2.table='dprime';
  guidata(handles2.figure2,handles2);

  if(handles.dump2csv)
    fid=fopen('most_recent_table.csv','w');
    fprintf(fid,'%% group, sample size, group 2, sample size 2, behavior, feature, d-prime, d-prime all-frames\n');
    transpose(tmp);
    fprintf(fid,'%s, %s, %s, %s, %s, %s, %s, %s\n',ans{:});
    fclose(fid);
  end
end

if(handles.comparison2==1)
%   if(nbehaviors>0)
%     mu2=nanmean(reshape(shiftdim(table_data(:,:,:,1:3),3),numel(table_data(:,:,:,1:3))/nfeatures,nfeatures),1);
%     sigma2=nanstd(reshape(shiftdim(table_data(:,:,:,1:3),3),numel(table_data(:,:,:,1:3))/nfeatures,nfeatures),1);
%   else
%     mu=nanmean(reshape(shiftdim(table_data(:,:,:,3),3),numel(table_data(:,:,:,3))/nfeatures,nfeatures),1);
%     sigma=nanstd(reshape(shiftdim(table_data(:,:,:,3),3),numel(table_data(:,:,:,3))/nfeatures,nfeatures),1);
%   end
  mu=nanmean(squeeze(table_data(:,1,:,3)));
%   sigma=nanmean(squeeze(table_data(:,1,:,6)));
  sigma=nanstd(squeeze(table_data(:,1,:,3)));
  tmp2=[];
  for g=1:length(handles.grouplist)
    gg=cumsum_num_selexp_per_group(g)+(1:length(handles.experimentlist{g}));
    if(nbehaviors>0)
      mu2=repmat(shiftdim(mu,-1),[length(gg) nbehaviors 1]);
      sigma2=repmat(shiftdim(sigma,-1),[length(gg) nbehaviors 1]);
      tmp2=[tmp2; ...
          repmat(g,nbehaviors*nfeatures,1) ...
          reshape(repmat((1:nbehaviors)',1,nfeatures),nbehaviors*nfeatures,1) ...
          reshape(repmat(1:nfeatures,nbehaviors,1),nbehaviors*nfeatures,1) ...
          reshape(sum(squeeze(table_data(gg,:,:,7))),nbehaviors*nfeatures,1) ...
          reshape(mean((squeeze(table_data(gg,:,:,1))-mu2)./sigma2),nbehaviors*nfeatures,1) ...
          reshape(std((squeeze(table_data(gg,:,:,1))-mu2)./sigma2),nbehaviors*nfeatures,1)];
      tmp2=[tmp2; ...
          repmat(g,nbehaviors*nfeatures,1) ...
          reshape(repmat(-(1:nbehaviors)',1,nfeatures),nbehaviors*nfeatures,1) ...
          reshape(repmat(1:nfeatures,nbehaviors,1),nbehaviors*nfeatures,1) ...
          reshape(sum(squeeze(table_data(gg,:,:,8))),nbehaviors*nfeatures,1) ...
          reshape(mean((squeeze(table_data(gg,:,:,2))-mu2)./sigma2),nbehaviors*nfeatures,1) ...
          reshape(std((squeeze(table_data(gg,:,:,2))-mu2)./sigma2),nbehaviors*nfeatures,1)];
    end
    mu2=repmat(mu,[length(gg) 1]);
    sigma2=repmat(sigma,[length(gg) 1]);
    tmp2=[tmp2; ...
        repmat(g,nfeatures,1) ...
        zeros(nfeatures,1) ...
        (1:nfeatures)' ...
        sum(squeeze(table_data(gg,1,:,9)))' ...
        mean((squeeze(table_data(gg,1,:,3))-mu2)./sigma2)' ...
        std((squeeze(table_data(gg,1,:,3))-mu2)./sigma2)'];
  end

  if(handles.omitnan)
    idx=find(~isnan(tmp2(:,5)));
    tmp2=tmp2(idx,:);
  end
  if(handles.omitinf)
    idx=find(~isinf(tmp2(:,5)));
    tmp2=tmp2(idx,:);
  end
  if(handles.absdprimezscore)
    tmp2(:,5)=abs(tmp2(:,5));
  end
  tmp2=sortrows(tmp2,-5);

  if(isempty(tmp2))
    set(handles.Status,'string','Ready.','foregroundcolor','g');
    set(handles.figure1,'pointer','arrow');
    uiwait(errordlg('no valid data.  check individual sex.'));  drawnow;
    return;
  end

  tmp=cell(size(tmp2,1),5);
  tmp(:,1)=handles.grouplist(tmp2(:,1));
  idx=(tmp2(:,2)>0);   tmp(idx,2)=handles.behaviorlist(tmp2(idx,2));
  idx=(tmp2(:,2)==0);  tmp(idx,2)=cellstr('all frames');
  idx=(tmp2(:,2)<0);   tmp(idx,2)=cellstr([repmat('not ',sum(idx),1) char(handles.behaviorlist(-tmp2(idx,2)))]);
  tmp(:,3)=handles.featurelist(tmp2(:,3));
  tmp(:,4)=cellstr(num2str(tmp2(:,4),'%-d'));
  tmp(:,5)=cellstr(num2str(tmp2(:,5),'%0.3g'));

  handles2.figure1=handles.figure1;
  handles2.figure2=figure('menubar','none','toolbar','none','numbertitle','off',...
      'name','interesting feature histograms');
  handles2.Table=uitable('Data',tmp,...
      'ColumnName',{'Group' 'Behavior' 'Feature' 'n' 'z'},...
      'ColumnWidth',num2cell(7*max(cellfun(@length,tmp))),...
      'RowName',[],'RowStriping','on','BackgroundColor',[1 1 1; 0.95 0.95 0.95],...
      'CellSelectionCallback',@CellSelectionCallback);
  extT=get(handles2.Table,'extent');
  posF=get(handles2.figure2,'position');
  %set(hf,'position',[posF(1) posF(2) 16*(posF(4)<extT(4))+extT(3) posF(4)]);
  set(handles2.figure2,'position',[posF(1) posF(2) min(posF(3),16*(posF(4)<extT(4))+extT(3)) min(posF(4),extT(4))]);
  %set(ht,'position',[0 0 16*(posF(4)<extT(4))+extT(3) posF(4)]);
  set(handles2.Table,'units','normalized','position',[0 0 1 1]);
  handles2.table_data=tmp2;
  handles2.table='zscore';
  guidata(handles2.figure2,handles2);

  if(handles.dump2csv)
    fid=fopen('most_recent_table.csv','w');
    fprintf(fid,'%% group, behavior, feature, sample size, z-score\n');
    transpose(tmp);
    fprintf(fid,'%s, %s, %s, %s, %s\n',ans{:});
    fclose(fid);
  end
end

if(handles.dump2mat)
  save('most_recent_table.mat','tmp');
end

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% --- Executes on button press in InterestingFeatureTimeSeries.
function InterestingFeatureTimeSeries_Callback(hObject, eventdata, handles)
% hObject    handle to InterestingFeatureTimeSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- 
function data=calculate_entiretimeseries(behavior_data,feature_data,sexdata,individual,xoffset);

if(isempty(feature_data.data))
  data=[];
  return;
end

if(iscell(feature_data.data{1}))
  vals=unique([feature_data.data{:}]);
  if(length(vals)>2)  error('uhoh');  end
  feature_data.data=cellfun(@(x) strcmp(x,vals{1}),feature_data.data,'uniformoutput',false);
end

if(isempty(behavior_data))
  behavior_data.allScores.tStart=1;
  behavior_data.allScores.tEnd=max(cellfun(@length,feature_data.data))+1;
end

k=1;
if(xoffset==1)
  data=nan(length(feature_data.data),max(behavior_data.allScores.tEnd));
else
  data=nan(length(feature_data.data),max(behavior_data.allScores.tEnd)-min(behavior_data.allScores.tStart));
end
for i=1:length(feature_data.data)  % individual
  if((~isnan(individual)) && (i~=individual))  continue;  end
  if(sum(sexdata{i}) < length(sexdata{i})/2)  continue;  end
  switch(xoffset)
    case(1)
      behavior_data.allScores.tStart(i);
    case(2)
      1;
    case(3)
      behavior_data.allScores.tStart(i)-min(behavior_data.allScores.tStart)+1;
  end
  data(k,ans:(ans+length(feature_data.data{i})-1))=feature_data.data{i};
  k=k+1;
end


% --- 
function triggered_data=calculate_triggeredtimeseries(behavior_data,behavior_logic,behavior_data2, ...
    feature_data,sexdata,individual,timing,windowradius,subtractmean,behaviornot)

timing=2+xor(timing-2,behaviornot);

if(iscell(feature_data.data{1}))
  vals=unique([feature_data.data{:}]);
  if(length(vals)>2)  error('uhoh');  end
  feature_data.data=cellfun(@(x) strcmp(x,vals{1}),feature_data.data,'uniformoutput',false);
end

k=1;
triggered_data=[];
for i=1:length(behavior_data.allScores.t0s)  % individual
  if((~isnan(individual)) && (i~=individual))  continue;  end
  feature_data_padded=[nan(1,windowradius) feature_data.data{i} nan(1,windowradius)];
  sexdata_padded=[nan(1,windowradius) sexdata{i} nan(1,windowradius)];
  if(behavior_logic>1)
    idx2=[behavior_data2.allScores.t0s{i}'-behavior_data2.allScores.tStart(i) ...
       behavior_data2.allScores.t1s{i}'-behavior_data2.allScores.tStart(i)];
  end
  for j=1:length(behavior_data.allScores.t0s{i})  % bout
    switch(timing)
      case(2)
        idx=behavior_data.allScores.t0s{i}(j)-behavior_data.allScores.tStart(i);
      case(3)
        idx=behavior_data.allScores.t1s{i}(j)-behavior_data.allScores.tStart(i);
    end
    if((behavior_logic==2) && (diff(sum(idx2<idx))==0))   continue;  end
    if((behavior_logic==3) && (diff(sum(idx2<idx))==-1))  continue;  end
    if(behavior_logic>3) error('or & or not are not implemented yet');  end
    if(sum(sexdata_padded(idx+(0:(2*windowradius))+(3-timing))) < windowradius)  continue;  end
    triggered_data(k,:)=feature_data_padded(idx+(0:(2*windowradius))+(3-timing));
    k=k+1;
  end
end
if(k==1)
  triggered_data(k,:)=nan(1,2*windowradius+1);
end
if(subtractmean)
  triggered_data=triggered_data-...
      repmat(nanmean(triggered_data(:,1:windowradius),2),1,size(triggered_data,2));
end


% ---
function handles=feature_timeseries_plot(handles)

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

handlesexperimentlist=[handles.experimentlist{:}];

cumsum_num_indi_per_exp=[0 cumsum(handles.individuals_feature)];
cumsum_num_exp_per_group=[0 cumsum(cellfun(@length,handles.experimentlist))];
mat2cell(cumsum_num_exp_per_group(1:end-1),1,ones(1,length(cumsum_num_exp_per_group)-1));
cellfun(@(x,y) x+y,handles.experimentvalue,ans,'uniformoutput',false);
selected_exp=[ans{:}];
cumsum_num_selexp_per_group=[0 cumsum(cellfun(@length,handles.experimentvalue))];

%ggee=1:length(handlesexperimentlist);
ggee=selected_exp;
individual=handles.individualidx;
if isnumeric(individual)
  ggee=find(cumsum_num_indi_per_exp<individual,1,'last');
  individual=individual-cumsum_num_indi_per_exp(ggee);
end

fid=nan;  if(handles.dump2csv)  fid=fopen('most_recent_figure.csv','w');  end

handles.type='feature time series';

hf=figure('toolbar','figure');

bb=handles.behaviorvalue;
if(handles.featuretimeseries_style2==1)  bb=1;  end
if(bb==(length(handles.behaviorlist)+1))  bb=1:(bb-1);  end
if(strcmp(get(handles.BehaviorList,'enable'),'off'))  bb=0;  end

behavior_logic=handles.behaviorlogic;
score_file2=[];
if((length(bb)>1) || (bb>0))  score_file2=handles.scorefiles{handles.behaviorvalue2};  end
feature_value=handles.featurevalue;
feature_list=handles.featurelist;
%sexdata=handles.sexdata;
timing=handles.featuretimeseries_style2;
xoffset=handles.xoffset;
if((length(bb)==1) && (bb==0))
  timing=1;
  xoffset=2;
end
style=handles.featuretimeseries_style;
centraltendency=handles.centraltendency;
dispersion=handles.dispersion;
convolutionwidth=round(handles.convolutionwidth*handles.fps);
subtractmean=handles.subtractmean;
windowradius=handles.windowradius;
behaviornot=handles.behaviornot;

h=[];
for b=bb

  if(length(bb)>1)
    ceil(sqrt(length(bb)));
    ha=subplot(ceil(length(bb)/ans),ans,b,'parent',hf);
  else
    ha=subplot(1,1,1,'parent',hf);
  end
  hold(ha,'on');

  score_file=[];
  if(b>0)  score_file=handles.scorefiles{b};  end

  num_indi=0;
  raw_data=cell(length(ggee),length(individual));
  data=cell(length(ggee),length(individual));
  parfor gei=1:length(ggee)
  %for gei=1:length(ggee)
    ge = ggee(gei);

    if(b>0)
      behavior_data=load(fullfile(handlesexperimentlist{ge},score_file));
      behavior_data=update_t01s_from_postprocessed(behavior_data);
      if(behavior_logic>1)
        behavior_data2=load(fullfile(handlesexperimentlist{ge},score_file2));
        behavior_data2=update_t01s_from_postprocessed(behavior_data2);
      else
        behavior_data2=[];
      end
    else
      behavior_data=[];
      behavior_data2=[];
    end
    feature_data=load(fullfile(handlesexperimentlist{ge},handles.perframe_dir,...
        [feature_list{feature_value} '.mat']));

    [behavior_data,behavior_data2,~,feature_data,sex_data]=...
        cull_short_trajectories(handles,behavior_data,behavior_data2,[],feature_data,handles.sexdata{ge});
    num_indi=num_indi+length(feature_data.data);

    ii=0;
    raw_parfor_tmp=cell(1,length(individual));
    parfor_tmp=cell(1,length(individual));
    for i = individual
      ii=ii+1;
      if(iscell(i))  i=char(i);  end
      tmp2=[];
      switch(i)
        case('M')
          tmp2=sex_data;
        case('F')
          tmp2=cellfun(@not,sex_data,'uniformoutput',false);
        otherwise
          tmp2=cellfun(@(x) ones(1,length(x)),sex_data,'uniformoutput',false);
      end
      tmp=nan;  if isnumeric(i)  tmp=i;  end

      if(timing==1)
        calculate_entiretimeseries(behavior_data,feature_data,tmp2,tmp,xoffset);
        raw_parfor_tmp{ii}=nanmean(ans,1);
        if(~isempty(raw_parfor_tmp{ii}))
          conv(raw_parfor_tmp{ii},ones(1,convolutionwidth),'valid');
          ans./conv(ones(1,length(raw_parfor_tmp{ii})),ones(1,convolutionwidth),'valid');
          parfor_tmp{ii}=[nan(1,floor((convolutionwidth-1)/2)) ans nan(1,ceil((convolutionwidth-1)/2))];
        else
          parfor_tmp{ii}=raw_parfor_tmp{ii};
        end
      else
        calculate_triggeredtimeseries(behavior_data,behavior_logic,behavior_data2,...
            feature_data,tmp2,tmp,timing,windowradius,subtractmean,behaviornot);
        raw_parfor_tmp{ii}=nanmean(ans,1);
        parfor_tmp{ii}=raw_parfor_tmp{ii};
      end
    end
    raw_data(gei,:)=raw_parfor_tmp;
    data(gei,:)=parfor_tmp;
  end
  raw_data=reshape(raw_data,1,prod(size(raw_data)));
  data=reshape(data,1,prod(size(data)));

  if(num_indi==0)
    delete(hf);
    set(handles.Status,'string','Ready.','foregroundcolor','g');
    set(handles.figure1,'pointer','arrow');
    uiwait(errordlg('no valid data.  check minimum trajectory length.'));  drawnow;
    return;
  end

  if(timing==1)
    max(cellfun(@(x) size(x,2),data));
    cellfun(@(x) [x nan(size(x,1),ans-size(x,2))],data,'uniformoutput',false);
    ydata=cat(1,ans{:});
    xdata=1:size(ydata,2);
  else
    ydata=cat(1,data{:});
    xdata=-windowradius:windowradius;
  end

  tstr='';
  if(timing>1)
    tstr='';  if(behaviornot)  tstr='NOT ';  end
    tstr=[tstr char(strrep(handles.behaviorlist(b),'_','-'))];
    switch(handles.behaviorlogic)
      case 2
        tstr=[tstr ' AND '];
      case 3
        tstr=[tstr ' AND NOT '];
%      case 4
%        tstr=[tstr ' OR '];
%      case 5
%        tstr=[tstr ' OR NOT '];
    end
    if(behavior_logic>1)
      tstr=[tstr char(strrep(handles.behaviorlist(handles.behaviorvalue2),'_','-'))];
    end
  end
  time_base=xdata./handles.fps;
  xstr='time (sec)';
  if(time_base(end)>300)
    time_base=time_base./60;
    xstr='time (min)';
  end
  if(time_base(end)>300)
    time_base=time_base./60;
    xstr='time (hr)';
  end
  if(time_base(end)>120)
    time_base=time_base./24;
    xstr='time (d)';
  end
  units=load(fullfile(handlesexperimentlist{ggee(1)},handles.perframe_dir,...
      [feature_list{feature_value} '.mat']),'units');
  ystr=get_label(feature_list(feature_value),units.units);

  if(handles.dump2csv)  print_csv_help(fid,handles.type,tstr,xstr,ystr);  end

  ii=0;
  for i = individual
    ii=ii+1;
    if(iscell(i))  i=char(i);  end
    for g=1:length(handles.grouplist)
      color=handles.colors(g,:);

      if ischar(i)
        idx=(cumsum_num_selexp_per_group(g)+1):(cumsum_num_selexp_per_group(g+1));
      else
        find(cumsum_num_exp_per_group<ggee,1,'last');
        if(ans~=g)  continue;  end
        idx=1;
      end
      idx2=idx+(ii-1)*numel(ggee);
      linestyle='-';  if(ii>1)  linestyle='--';  end

      if(handles.dump2csv)  fprintf(fid,['%% group ' handles.grouplist{g} '\n']);  end
      plot_it(ha,time_base,ydata(idx2,:),style,centraltendency,dispersion,color,1,linestyle,...
          fid,handlesexperimentlist(idx));
      if (ii==1)
        h(g)=ans;
        hh{g}=handles.grouplist{g};
      end
    end

    if(handles.dump2csv)  fprintf(fid,'\n%% raw data\n');  end
    for g=1:length(handles.grouplist)
      if ischar(i)
        idx=(cumsum_num_selexp_per_group(g)+1):(cumsum_num_selexp_per_group(g+1));
      else
        find(cumsum_num_exp_per_group<ggee,1,'last');
        if(ans~=g)  continue;  end
        idx=1;
      end
      idx2=idx+(ii-1)*numel(ggee);

      if(handles.dump2csv)
        fprintf(fid,['%% group ' handles.grouplist{g} '\n']);
        for e=1:length(idx2)
          fprintf(fid,'%% experiment %s\n',handlesexperimentlist{selected_exp(idx(e))});
          print_csv_data(fid,raw_data{idx2(e)});
          fprintf(fid,'\n');
        end
      end
    end
  end

  if(handles.dump2mat)
    find(b==bb);
    raw_data2{ans}=raw_data;
  end

  xlabel(ha,xstr,'interpreter','none');
  ylabel(ha,ystr,'interpreter','none');
  title(ha,tstr,'interpreter','none');
  axis(ha,'tight');  zoom(ha,'reset');
end

if(handles.dump2mat)
  raw_data=raw_data2;
  save('most_recent_figure.mat','handles','raw_data');
end

if(iscell(individual))
  h(end+1)=plot(0,0,'k-');   hh{end+1}='males';
  h(end+1)=plot(0,0,'k--');  hh{end+1}='females';
  set(h((end-1):end),'visible','off');
end

idx=find(h>0);
if ~isnumeric(individual)
  %legend(ha,h(idx),[cellfun(@(x) [x ' ' handles.individuallist{handles.individualvalue}],...
  %    handles.grouplist,'uniformoutput',false)],'interpreter','none');
  legend(ha,h(idx),hh(idx),'interpreter','none');
%else
%  legend(ha,h(idx),handles.individuallist(handles.individualvalue),'interpreter','none');
end

uicontrol(hf,'style','pushbutton','string','Params','position',[5 5 60 20],...
    'callback',@figure_params_callback);

if(handles.dump2csv)  fclose(fid);  end

guidata(hf,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% ---
function table_data=calculate_interesting_timeseries(experiment_value,experiment_list,...
    behavior_list,feature_list,perframe_dir,windowradius)

table_data=cell(length(behavior_list),length(feature_list),2);
parfor b=1:length(behavior_list)
%for b=1:length(behavior_list)
  parfor_tmp=cell(length(feature_list),2);
  for f=1:length(feature_list)
    for t=2:3  % timing
      parfor_tmp{f,t-1}=[];
      for e=1:length(experiment_value)
        behavior_data=load(fullfile(experiment_list{experiment_value(e)},[behavior_list{b} '.mat']));
        behavior_data=update_t01s_from_postprocessed(behavior_data);
        feature_data=load(fullfile(experiment_list{experiment_value(e)},perframe_dir,...
          [feature_list{f} '.mat']));
        sexdata={};
        for i=1:length(feature_data.data)
          sexdata{i}=ones(1,length(feature_data.data{i}));
        end

        tmp=calculate_triggeredtimeseries(behavior_data,1,[],feature_data,sexdata,nan,t,windowradius,1);
        parfor_tmp{f,t-1}=[parfor_tmp{f,t-1}; tmp];
      end
      parfor_tmp{f,t-1}=[...
          sqrt(nanmean(parfor_tmp{f,t-1}(:,1:windowradius).^2,2))...
          sqrt(nanmean(parfor_tmp{f,t-1}(:,(windowradius+1):end).^2,2))];
    end
  end
  table_data(b,:,:)=parfor_tmp;
  disp([num2str(b) ' of ' num2str(length(behavior_list))]);
end


% --- Executes on button press in InterestingTimeSeries.
function InterestingTimeSeries_Callback(hObject, eventdata, handles)
% hObject    handle to InterestingTimeSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

experiment_value=get(handles.ExperimentList,'Value');
experiment_list=get(handles.ExperimentList,'String');
experiment_value2=get(handles.ExperimentList2,'Value');
experiment_list2=get(handles.ExperimentList2,'String');
behavior_list=get(handles.BehaviorList,'String');
feature_list=get(handles.FeatureList,'String');
windowradius=handles.windowradius;

if(isempty(handles.interestingfeaturetimeseries_cache))
  if(length(experiment_value)>0)
    table_data=calculate_interesting_timeseries(experiment_value,experiment_list,...
        behavior_list,feature_list,handles.perframe_dir,windowradius);
  end
  if(length(experiment_value2)>0)
    table_data2=calculate_interesting_timeseries(experiment_value2,experiment_list2,...
        behavior_list,feature_list,handles.perframe_dir,windowradius);
  end

  if((length(experiment_value)>0) && (length(experiment_value2)>0))
    tmp=(cellfun(@(x) nanmean(x(:,2)),table_data2) - cellfun(@(x) nanmean(x(:,2)),table_data)) ./ ...
        sqrt(cellfun(@(x) nanstd(x(:,2)),table_data2).^2 + cellfun(@(x) nanstd(x(:,2)),table_data).^2);
  elseif(length(experiment_value)>0)
    tmp=cellfun(@(x) (nanmean(x(:,2))-nanmean(x(:,1)))/sqrt(nanstd(x(:,2)).^2+nanstd(x(:,1)).^2),table_data);
  elseif(length(experiment_value2)>0)
    tmp=cellfun(@(x) (nanmean(x(:,2))-nanmean(x(:,1)))/sqrt(nanstd(x(:,2)).^2+nanstd(x(:,1)).^2),table_data2);
  end
  [tmp2(:,1) tmp2(:,2) tmp2(:,3)]=ind2sub(size(tmp),1:prod(size(tmp)));
  tmp2(:,4)=tmp(1:prod(size(tmp)));
  tmp2=sortrows(tmp2,-4);

  handles.interestingfeaturetimeseries_cache=tmp2;
else
  tmp2=handles.interestingfeaturetimeseries_cache;
end

tmp=cell(size(tmp2));
tmp(:,1)=behavior_list(tmp2(:,1));
tmp(:,2)=feature_list(tmp2(:,2));
tmp(:,3)=repmat({'Onset'},size(tmp2,1),1);
  find(tmp2(:,3)==2);
  tmp(ans,3)=repmat({'Offset'},size(ans),1);
tmp(:,4)=num2cell(tmp2(:,4));

set(handles.Table,'Data',tmp);
set(handles.Table,'ColumnName',{'Behavior' 'Feature' 'Timing' 'dRMS'});
set(handles.Table,'ColumnWidth',{150 150 50 50 50});
set(handles.Table,'RowStriping','on','BackgroundColor',[1 1 1; 0.95 0.95 0.95]);

handles.table_data=tmp2;
handles.table='timeseries';
guidata(hObject,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% ---
function ret_val=print_stats(data,stat)

switch(stat)
  case(0)
    ret_val=num2str(data,3);
  case(1)
    ret_val=num2str(mean(data),3);
  case(2)
    ret_val=[num2str(mean(data),3) ' (' num2str(std(data),3) ')'];
  case(3)
    ret_val=[num2str(mean(data),3) ' (' num2str(std(data)./sqrt(length(data)),3) ')'];
  case(4)
    ret_val=num2str(median(data),3);
  case(5)
    ret_val=[num2str(median(data),3) ' (' num2str(prctile(data,5),3) '-' num2str(prctile(data,95),3) ')'];
  case(6)
    ret_val=[num2str(median(data),3) ' (' num2str(prctile(data,25),3) '-' num2str(prctile(data,75),3) ')'];
end


% ---
function [ct,dp,dn]=calculate_ct_d(data,centraltendency,dispersion)

switch(centraltendency)
  case 1
    ct=nanmean(data);
  case 2
    ct=nanmedian(data);
  case 3
    ct=mode(data);
end
switch(dispersion)
  case 1
    tmp=std(data);
    dp=ct+tmp;
    dn=ct-tmp;
  case 2
    tmp=std(data)./sqrt(length(data));
    dp=ct+tmp;
    dn=ct-tmp;
  case 3
    dp=prctile(data,95);
    dn=prctile(data,5);
  case 4
    dp=prctile(data,75);
    dn=prctile(data,25);
end


% ---
function [h,h2]=errorbarplot(ha,x,b,dn,dp,color)

h = bar(ha,x,b,'grouped');
set(h,'facecolor',color);
fxdata = get(h, 'XData');
fydata = get(h ,'YData');
if(~iscell(fxdata)) xdata{1,1} = fxdata; else xdata = fxdata; end
if(~iscell(fydata)) ydata{1,1} = fydata; else ydata = fydata; end

% Determine number of bars
sizz = size(b);
nb = sizz(1)*sizz(2);
xb = [];
yb = [];
for i = 1:length(xdata),
    xb = [xb xdata{i,1}];
    yb = [yb ydata{i,1}];
end;

%% To find the center of each bar - need to look at the output vectors xb, yb
%% find where yb is non-zero - for each bar there is a pair of non-zero yb values.
%% The center of these values is the middle of the bar
%
%nz = find(yb);
%for i = 1:nb,
%    center(i) = (xb(nz(i*2))-xb(nz((i*2)-1)))/2 + xb(nz((i*2)-1));
%end;

% To place the error bars - use the following:

%errdata=[.1 .2; .3 .4; .5 .6];

hold(ha,'on');
h2=errorbar(ha,xb, b, dn, dp);
set(h2(1),'linewidth',1);            % This changes the thickness of the errorbars
set(h2(1),'color','k');              % This changes the color of the errorbars
set(h2(1),'linestyle','none');       % This removes the connecting


% ---
function ret_val=color_red(arg,flag)

ret_val=num2str(arg);
if(flag)
  ret_val=['<html><font color="red">', ret_val, '</font></html>'];
end


% ---
function tmp=calculate_statistics(table_data,behaviorlist,grouplist,fid,crit)

nrows=14;  if(length(table_data{1})==2)  nrows=8;  end

tmp={'behavior'};
for b=1:length(table_data)

  grouplistcopy=grouplist;
  g=1;
  while g<=length(table_data{b})
    if(isempty(table_data{b}{g}) || (sum(~isnan(table_data{b}{g}))==0))
      table_data{b}(g)=[];
      grouplistcopy(g)=[];
    else
      g=g+1;
    end
  end

  tmp{nrows*b+1,1}=behaviorlist{b};
  for g=1:length(table_data{b})
    if(b==1)  tmp(2:5,g)=repmat({grouplistcopy{g}},4,1);  end
    p=nan;
    if(sum(~isnan(table_data{b}{g}))>0)
      [~,p,~,~]=kstest(table_data{b}{g});
    end
    tmp{nrows*b+2,g}=color_red(length(table_data{b}{g}),0);
    tmp{nrows*b+3,g}=color_red(p,p<crit);
    tmp{nrows*b+4,g}=color_red(nanmean(table_data{b}{g}),0);
    tmp{nrows*b+5,g}=color_red(nanstd(table_data{b}{g}),0);
  end
  if(b==1)
    tmp{2,length(table_data{b})+1}='(sample size)';
    tmp{3,length(table_data{b})+1}='(K-S normal)';
    tmp{4,length(table_data{b})+1}='(mean)';
    tmp{5,length(table_data{b})+1}='(std. dev.)';
  end

  if(length(table_data{b})==2)
  
    if(b==1)  tmp{6,1}='t-test';  end
    [~,p]=ttest2(table_data{b}{1},table_data{b}{2});
    tmp{nrows*b+6,1}=color_red(p,p<crit);

    if(b==1)  tmp{7,1}='Wilcoxen';  end
    p=ranksum(table_data{b}{1},table_data{b}{2});
    tmp{nrows*b+7,1}=color_red(p,p<crit);

  else

    foo=cellfun(@(x,y) y*ones(1,length(x)),table_data{b},num2cell(1:length(table_data{b})),...
        'uniformoutput',false);
    [p,table,stats]=anova1([table_data{b}{:}],[foo{:}],'off');
    c=multcompare(stats,'alpha',crit,'display','off');
    tmp{nrows*b+6,1}=color_red(p,p<crit);
    if(b==1)  tmp{6,1}='ANOVA';  end
    for g2=1:size(c,1)
      if(b==1)
        tmp(7:9,g2)=repmat({strcat(grouplistcopy{c(g2,1)},'-',grouplistcopy{c(g2,2)})},3,1);
      end
      foo=(sign(c(g2,3))==sign(c(g2,5)));
      tmp{nrows*b+7,g2}=color_red(c(g2,3),foo);
      tmp{nrows*b+8,g2}=color_red(c(g2,4),foo);
      tmp{nrows*b+9,g2}=color_red(c(g2,5),foo);
    end
    if(b==1)
      tmp(7:9,size(c,1)+1)={'(5%, Tukey post-hoc)'; '(mean)'; '(95%)'};
    end

    foo=cellfun(@(x,y) y*ones(1,length(x)),table_data{b},num2cell(1:length(table_data{b})),...
        'uniformoutput',false);
    [p,table,stats]=kruskalwallis([table_data{b}{:}],[foo{:}],'off');
    c=multcompare(stats,'alpha',crit,'display','off');
    tmp{nrows*b+10,1}=color_red(p,p<crit);
    if(b==1)  tmp{10,1}='Kruskal-Wallis';  end
    for g2=1:size(c,1)
      if(b==1)
        tmp(11:13,g2)=repmat({strcat(grouplistcopy{c(g2,1)},'-',grouplistcopy{c(g2,2)})},3,1);
      end
      foo=(sign(c(g2,3))==sign(c(g2,5)));
      tmp{nrows*b+11,g2}=color_red(c(g2,3),foo);
      tmp{nrows*b+12,g2}=color_red(c(g2,4),foo);
      tmp{nrows*b+13,g2}=color_red(c(g2,5),foo);
    end
    if(b==1)
      tmp(11:13,size(c,1)+1)={'(5%, Tukey post-hoc)'; '(mean)'; '(95%)'};
    end
  end
end
tmp{end+1,1}='';

if(~isnan(fid))
  fprintf(fid,'\n%% statistics\n');
  for j=1:(size(tmp,1)/nrows)
    for i=1:nrows
      if(j==1)  fprintf(fid,'%% ');  end
      tmp(nrows*(j-1)+i,:);
      ans(cellfun(@isempty,ans))='';
      cellfun(@(x) regexprep(x,'<[^<]*>','*'),ans,'uniformoutput',false);
      fprintf(fid,'%s, ',ans{:});
      fprintf(fid,'\n');
    end
  end
end


% ---
function handles=behavior_barchart_plot(handles)

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

handlesexperimentlist=[handles.experimentlist{:}];

cumsum_num_indi_per_exp=[0 cumsum(handles.individuals_feature)];
cumsum_num_exp_per_group=[0 cumsum(cellfun(@length,handles.experimentlist))];
mat2cell(cumsum_num_exp_per_group(1:end-1),1,ones(1,length(cumsum_num_exp_per_group)-1));
cellfun(@(x,y) x+y,handles.experimentvalue,ans,'uniformoutput',false);
selected_exp=[ans{:}];
cumsum_num_selexp_per_group=[0 cumsum(cellfun(@length,handles.experimentvalue))];

%ggee=1:length(handlesexperimentlist);
ggee=selected_exp;
individual=handles.individualidx;
if isnumeric(individual)
  ggee=find(cumsum_num_indi_per_exp<individual,1,'last');
  individual=individual-cumsum_num_indi_per_exp(ggee);
end

fid=nan;  if(handles.dump2csv)  fid=fopen('most_recent_figure.csv','w');  end

handles.type='behavior bar chart';

hf=figure('toolbar','figure');

bb=handles.behaviorvalue;
if(bb==(length(handles.behaviorlist)+1))  bb=1:(bb-1);  end

behavior_logic=handles.behaviorlogic;
score_file2=handles.scorefiles{handles.behaviorvalue2};
score_file3=[];
if(handles.behaviorvalue3>1)
  score_file3=handles.scorefiles{handles.behaviorvalue3-1};
end
%sexdata=handles.sexdata;
%perwhat=handles.behaviorbarchart_style;
behaviornot=handles.behaviornot;

h={};
table_data={};
for b=bb

  tstr='';
  ystr='';  if(behaviornot)  ystr='NOT ';  end
  ystr=[ystr char(strrep(handles.behaviorlist(b),'_','-'))];
  switch(handles.behaviorlogic)
    case 2
      ystr=[ystr ' AND '];
    case 3
      ystr=[ystr ' AND NOT '];
    case 4
      ystr=[ystr ' OR '];
    case 5
      ystr=[ystr ' OR NOT '];
  end
  if(handles.behaviorlogic>1)
    ystr=[ystr char(strrep(handles.behaviorlist(handles.behaviorvalue2),'_','-'))];
  end
  ystr=[ystr ' (%)'];
  xstr='group';

  if(handles.dump2csv)  print_csv_help(fid,handles.type,tstr,xstr,ystr);  end

  if(length(bb)>1)
    ceil(sqrt(length(bb)));
    ha=subplot(ceil(length(bb)/ans),ans,b,'parent',hf);
  else
    ha=subplot(1,1,1,'parent',hf);
  end
  hold(ha,'on');

  score_file=handles.scorefiles{b};

  num_indi=0;
  collated_data=cell(1,length(ggee));
  parfor gei=1:numel(ggee)
  %for gei=1:numel(ggee)
    ge = ggee(gei);

    %if(ischar(individual)&&(~ismember(ge,selected_exp)))  continue;  end

    behavior_data=load(fullfile(handlesexperimentlist{ge},score_file));
    behavior_data=update_t01s_from_postprocessed(behavior_data);
    behavior_data2=[];
    if(behavior_logic>1)
      behavior_data2=load(fullfile(handlesexperimentlist{ge},score_file2));
      behavior_data2=update_t01s_from_postprocessed(behavior_data2);
    end
    behavior_data3=[];
    if(handles.behaviorvalue3>1)
      behavior_data3=load(fullfile(handlesexperimentlist{ge},score_file3));
      behavior_data3=update_t01s_from_postprocessed(behavior_data3);
    end

    [behavior_data,behavior_data2,behavior_data3,~,sex_data]=...
        cull_short_trajectories(handles,behavior_data,behavior_data2,behavior_data3,[],handles.sexdata{ge});
    num_indi=num_indi+length(behavior_data.allScores.scores);

    traj_len=behavior_data.allScores.tEnd-behavior_data.allScores.tStart;

    frames_labelled=nan(1,length(behavior_data.allScores.t0s)); %#ok<PFTUS>
    frames_total=nan(1,length(behavior_data.allScores.t0s));
    sex=nan(1,length(behavior_data.allScores.t0s));

    for i=1:length(behavior_data.allScores.t0s)  % individual
      tmp1 = compute_behavior_logic(behavior_data.allScores, i);
      tmp1 = tmp1(behavior_data.allScores.tStart(i) : behavior_data.allScores.tEnd(i));

      tmp2=[];
      if(behavior_logic>1)
        tmp2 = compute_behavior_logic(behavior_data2.allScores, i);
        tmp2 = tmp2(behavior_data2.allScores.tStart(i) : behavior_data2.allScores.tEnd(i));
      end

      tmp3=[];
      if(handles.behaviorvalue3>1)
        tmp3 = compute_behavior_logic(behavior_data3.allScores, i);
        tmp3 = tmp3(behavior_data3.allScores.tStart(i) : behavior_data3.allScores.tEnd(i));
        if(handles.behaviornormalizenot)
          tmp3=~tmp3;
        end
      end

      if(behaviornot)  tmp1=~tmp1;  end

      partition_idx=[];
      switch(behavior_logic)
        case(1)
          partition_idx=tmp1;
        case(2)
          partition_idx=tmp1 & tmp2;
        case(3)
          partition_idx=tmp1 & ~tmp2;
        case(4)
          partition_idx=tmp1 | tmp2;
        case(5)
          partition_idx=tmp1 | ~tmp2;
      end

      sex(i)=sum(sex_data{i}(1:length(partition_idx))) > (length(partition_idx)/2);
      frames_labelled(i)=sum(partition_idx);
      if(handles.behaviorvalue3==1)
        frames_total(i)=length(partition_idx);
      else
        frames_total(i)=sum(tmp3);
      end
    end

    collated_data{gei}={frames_labelled frames_total sex traj_len};
  end

  %idx=cellfun(@isempty,collated_data);
  %collated_data=collated_data(~idx);

  if(num_indi==0)
    delete(hf);
    set(handles.Status,'string','Ready.','foregroundcolor','g');
    set(handles.figure1,'pointer','arrow');
    uiwait(errordlg('no valid data.  check minimum trajectory length.'));  drawnow;
    return;
  end

  exp_separators=[];  maxy=0;  k=[];  m=0;  ii=0;  table_data{end+1}=[];
  for i = individual
    if(iscell(i))  i=char(i);  end
    switch(i)
      case 'A'
        frames_labelled=cellfun(@(x) x{1},collated_data,'uniformoutput',false);
        frames_total=cellfun(@(x) x{2},collated_data,'uniformoutput',false);
        traj_len=cellfun(@(x) x{4},collated_data,'uniformoutput',false);
      case 'M'
        frames_labelled=cellfun(@(x) x{1}(x{3}==1),collated_data,'uniformoutput',false);
        frames_total=cellfun(@(x) x{2}(x{3}==1),collated_data,'uniformoutput',false);
        traj_len=cellfun(@(x) x{4}(x{3}==1),collated_data,'uniformoutput',false);
      case 'F'
        frames_labelled=cellfun(@(x) x{1}(x{3}==0),collated_data,'uniformoutput',false);
        frames_total=cellfun(@(x) x{2}(x{3}==0),collated_data,'uniformoutput',false);
        traj_len=cellfun(@(x) x{4}(x{3}==0),collated_data,'uniformoutput',false);
      otherwise
        frames_labelled=cellfun(@(x) x{1}(i),collated_data,'uniformoutput',false);
        frames_total=cellfun(@(x) x{2}(i),collated_data,'uniformoutput',false);
        traj_len=cellfun(@(x) x{4}(i),collated_data,'uniformoutput',false);
    end

    if(handles.dump2csv)
      fprintf(fid,'\n%% raw data\n');
      for g=1:length(handles.grouplist)
        if ~isnumeric(individual)
          idx=(cumsum_num_selexp_per_group(g)+1):(cumsum_num_selexp_per_group(g+1));
        else
          find(cumsum_num_exp_per_group<ggee,1,'last');
          if(ans~=g)  continue;  end
          idx=1;
        end
        fprintf(fid,['%% individual ' i '\n']);
        fprintf(fid,['%% group ' handles.grouplist{g} '\n']);
        for e=1:length(idx)
          fprintf(fid,'%% experiment %s\n',handlesexperimentlist{selected_exp(idx(e))});
          fprintf(fid,'%% frames labelled\n');
          fprintf(fid,'%d, ',frames_labelled{idx(e)});
          fprintf(fid,'\n');
          fprintf(fid,'%% frames total\n');
          fprintf(fid,'%d, ',frames_total{idx(e)});
          fprintf(fid,'\n');
          fprintf(fid,'\n');
        end
      end
    end

    for g=1:length(handles.grouplist)
      color=handles.colors(g,:);
      ii=ii+1;

      if ischar(i)
        idx=(cumsum_num_selexp_per_group(g)+1):(cumsum_num_selexp_per_group(g+1));
      else
        find(cumsum_num_exp_per_group<ggee,1,'last');
        if(ans~=g)  continue;  end
        idx=1;
      end

      xticklabels{ii}=handles.grouplist{g};

      switch(handles.behaviorbarchart_style)
        case 1  % per group
          table_data{end}(ii)=100*sum([frames_labelled{idx}])./sum([frames_total{idx}]);
          h{ii}=bar(ha,ii,table_data{end}(ii));
          set(h{ii},'facecolor',color);
        case 2  % per experiment, error bars
          table_data{end}{ii}=100*cellfun(@sum,frames_labelled(idx))./cellfun(@sum,frames_total(idx));
          [ct(g),dp(g),dn(g)]=...
              calculate_ct_d(table_data{end}{ii},handles.centraltendency,handles.dispersion);
          h{ii}=errorbarplot(ha,ii,ct(g),ct(g)-dn(g),dp(g)-ct(g),color);
        case 3  % per experiment, box
          table_data{end}{ii}=100*cellfun(@sum,frames_labelled(idx))./cellfun(@sum,frames_total(idx));
          h{ii}=boxplot(ha,table_data{end}{ii},'positions',ii,'widths',0.5,'colors',color);
          findobj(h{ii},'tag','Outliers');
          set(ans,'markeredgecolor',color);
        case 4  % per fly, grouped
          cumsum(cellfun(@length,frames_labelled(idx)))';
          exp_separators=[exp_separators; ans+sum(k)];
          table_data{end}{ii}=100.*[frames_labelled{idx}]./[frames_total{idx}];
          maxy=max([maxy table_data{end}{ii}]);
          h{ii}=bar(ha,(1:length(table_data{end}{ii}))+sum(k),table_data{end}{ii},...
              'barwidth',1,'edgecolor','none');
          set(h{ii},'facecolor',color);
          k(end+1)=length(table_data{end}{ii});
        case 5  % per fly, stern-style
          table_data{end}{ii}=cell(1,length(frames_labelled(idx)));
          h{ii}=[];
          for e=1:length(idx)
            table_data{end}{ii}{e}=100.*frames_labelled{idx(e)}./frames_total{idx(e)};
            [ct,dp,dn]=calculate_ct_d(table_data{end}{ii}{e},...
                handles.centraltendency,handles.dispersion);
            plot(ha,m,ct,'o','color',color);
            h{ii}(end+1)=plot(ha,[m m],[dp dn],'-','color',color);
            plot(ha,m+(1:length(table_data{end}{ii}{e})),table_data{end}{ii}{e},'.','color',color,'markersize',12);
            m=m+16+length(table_data{end}{ii}{e});
          end
          [ct,dp,dn]=calculate_ct_d([table_data{end}{ii}{:}],...
              handles.centraltendency,handles.dispersion);
          plot(ha,m,ct,'o','color',color,'markersize',9);
          h{ii}(end+1)=plot(ha,[m m],[dp dn],'-','color',color,'linewidth',3);
          m=m+24;
          k(end+1)=24+16*length(table_data{end}{ii})+length([table_data{end}{ii}{:}]);
        case 6  % per fly, trajectory length
          cumsum(cellfun(@length,traj_len(idx)))';
          exp_separators=[exp_separators; ans+sum(k)];
          table_data{end}{ii}=[traj_len{idx}];
          maxy=max([maxy table_data{end}{ii}]);
          h{ii}=bar(ha,(1:length(table_data{end}{ii}))+sum(k),table_data{end}{ii},...
              'barwidth',1,'edgecolor','none');
          set(h{ii},'facecolor',color);
          k(end+1)=length(table_data{end}{ii});
      end
    end
  end

  if(ismember(handles.behaviorbarchart_style,[1 2 4 6]) && (length(individual)==2))
    for g=(length(h)/2+1):length(h)
      findobj([h{g}],'type','patch');
      hatchfill(ans,'single',45,5,handles.colors(g-length(h)/2,:));
    end
  end
  if(ismember(handles.behaviorbarchart_style,[3 5]) && (length(individual)==2))
    for g=1:length(h)
      findobj([h{g}],'type','line');
      if (g<=(length(h)/2))
        set(ans,'linestyle','-');
      else
        set(ans,'linestyle','--');
      end
    end
  end

  if(handles.dump2csv)
    fprintf(fid,'\n%% summary data\n');
    fprintf(fid,['%% xdata\n']);  fprintf(fid,'%s, ',xticklabels{:});  fprintf(fid,'\n');
  end
  switch(handles.behaviorbarchart_style)
    case 1  % per group
      if(handles.dump2csv)
        fprintf(fid,['%% ydata, per group\n']);
        fprintf(fid,'%g, ',table_data{end});
        fprintf(fid,'\n');
      end
    case 2  % per experiment, error bars
      if(handles.dump2csv)
        fprintf(fid,['%% ydata, CT+D\n']);  fprintf(fid,'%g, ',dp);  fprintf(fid,'\n');
        fprintf(fid,['%% ydata, CT-D\n']);  fprintf(fid,'%g, ',dn);  fprintf(fid,'\n');
        fprintf(fid,['%% ydata, CT\n']);    fprintf(fid,'%g, ',ct);  fprintf(fid,'\n');
      end
    case 3  % per experiment, box
      if(handles.dump2csv)
        fprintf(fid,['%% ydata, percentiles 1, 5, 25, 50, 75, 95, 99\n']);
        for i=1:length(table_data{end})
          fprintf(fid,'%g, ',prctile([table_data{end}{i}],[1 5 25 50 75 95 99]));
          fprintf(fid,'\n');
        end
      end
    case {4,6}  % per fly, grouped
      l=exp_separators(1:2:(end-1));
      r=exp_separators(2:2:end);
      hh=patch(0.5+[l r r l l]',repmat([0 0 maxy*1.05 maxy*1.05 0]',1,floor(length(exp_separators)/2)),...
          [0.95 0.95 0.95],'parent',ha);
      set(hh,'edgecolor','none');
      set(ha,'children',circshift(get(ha,'children'),-1));
      k=round(cumsum(k)-k/2);
      if(handles.dump2csv)
        fprintf(fid,['%% ydata\n']);
        for i=1:length(table_data{end})
          fprintf(fid,'%g, ',[table_data{end}{i}]);
          fprintf(fid,'\n');
        end
      end
    case 5  % per fly, stern-style
      k=round(cumsum(k)-k/2);
      if(handles.dump2csv)
        fprintf(fid,['%% ydata\n']);
        for i=1:length(table_data{end})
          for j=1:length(table_data{end}{i})
            fprintf(fid,'%g, ',[table_data{end}{i}{j}]);
            fprintf(fid,'\n');
          end
          fprintf(fid,'\n');
        end
      end
  end

  if(isempty(k))  k=1:length(frames_labelled);  end
  %if(length(individual)==2)  k=k*2-0.5;  end
  ylabel(ha,ystr,'interpreter','none');
  set(ha,'xtick',k,'xticklabel',xticklabels);
  axis(ha,'tight');  vt=axis;
  axisalmosttight([],ha);  vat=axis;
  if(handles.behaviorbarchart_style==4)
    axis(ha,[vat(1) vat(2) 0 vt(4)]);
  else
    axis(ha,[vat(1) vat(2) 0 vat(4)]);
  end

  if(handles.dump2mat)
    find(b==bb);
    raw_data{ans}=collated_data;
  end
end

if(handles.dump2mat)
  save('most_recent_figure.mat','handles','raw_data');
end

if(iscell(individual))
  h2(1)=plot(0,0,'k-');   hh2{1}='males';
  h2(2)=plot(0,0,'k--');  hh2{2}='females';
  set(h2,'visible','off');
  legend(ha,h2,hh2,'interpreter','none');
end
%idx=find(~cellfun(@isempty,h));
%if ischar(individual)
%  legend(cellfun(@(x) x(1),h(idx)),handles.grouplist);
%else
%  legend(cellfun(@(x) x(1),h(idx)),handles.individuallist(handles.individualvalue));
%end

%if((ismember(handles.behaviorbarchart_style,[2 3])) && (individual<4))
uicontrol(hf,'style','pushbutton','string','Params','position',[5 5 60 20],...
    'callback',@figure_params_callback);
if((ismember(handles.behaviorbarchart_style,[2 3 4])) && ischar(individual) && (length(handles.grouplist)>1))
  uicontrol(hf,'style','pushbutton','string','Stats','position',[70 5 50 20],...
      'callback',@figure_stats_callback);
  handles.statistics=calculate_statistics(table_data,handles.behaviorlist(bb),handles.grouplist,...
      fid,handles.pvalue);
%  set(handles.Table,'Data',tmp);
%  set(handles.Table,'ColumnWidth','auto');
%  set(handles.Table,'ColumnName',{''});
%else
%  set(handles.Table,'Data',[]);
%  set(handles.Table,'ColumnName',{});
end

if(handles.dump2csv)  fclose(fid);  end

guidata(hf,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% ---
% returns a logical time series of when the behavior occurred
function ret_val = compute_behavior_logic(allScores, i)

ret_val=zeros(1,allScores.tEnd(i));
ret_val(allScores.t0s{i})=1;
ret_val(allScores.t1s{i})=-1;
ret_val=logical(cumsum(ret_val));
ret_val=ret_val(1:allScores.tEnd(i));

        
% ---
function handles=behavior_timeseries_plot(handles)

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

handlesexperimentlist=[handles.experimentlist{:}];

cumsum_num_indi_per_exp=[0 cumsum(handles.individuals_feature)];
cumsum_num_exp_per_group=[0 cumsum(cellfun(@length,handles.experimentlist))];
mat2cell(cumsum_num_exp_per_group(1:end-1),1,ones(1,length(cumsum_num_exp_per_group)-1));
cellfun(@(x,y) x+y,handles.experimentvalue,ans,'uniformoutput',false);
selected_exp=[ans{:}];
cumsum_num_selexp_per_group=[0 cumsum(cellfun(@length,handles.experimentvalue))];

%ggee=1:length(handlesexperimentlist);
ggee=selected_exp;
individual=handles.individualidx;
if isnumeric(individual)
  ggee=find(cumsum_num_indi_per_exp<individual,1,'last');
  individual=individual-cumsum_num_indi_per_exp(ggee);
end

fid=nan;  if(handles.dump2csv)  fid=fopen('most_recent_figure.csv','w');  end

handles.type='behavior time series';

hf=figure('toolbar','figure');

bb=handles.behaviorvalue;
if(bb==(length(handles.behaviorlist)+1))  bb=1:(bb-1);  end

behavior_logic=handles.behaviorlogic;
score_file2=handles.scorefiles{handles.behaviorvalue2};
score_file3=[];
if(handles.behaviorvalue3>1)
  score_file3=handles.scorefiles{handles.behaviorvalue3-1};
end
%sexdata=handles.sexdata;
convolutionwidth=round(handles.convolutionwidth*handles.fps);
style=handles.behaviortimeseries_style;
centraltendency=handles.centraltendency;
dispersion=handles.dispersion;
xoffset=handles.xoffset;
behaviornot=handles.behaviornot;

h=[];  hh={};
for b=bb

  if(length(bb)>1)
    ceil(sqrt(length(bb)));
    ha=subplot(ceil(length(bb)/ans),ans,b,'parent',hf);
  else
    ha=subplot(1,1,1,'parent',hf);
  end
  hold(ha,'on');

  score_file=handles.scorefiles{b};

  num_indi=0;
  raw_data=cell(length(ggee),length(individual));
  behavior_cumulative=cell(length(ggee),length(individual));
  parfor gei=1:numel(ggee),
%   for gei=1:numel(ggee),
    ge = ggee(gei);
    
    %if(ischar(individual)&&(~ismember(ge,selected_exp)))  continue;  end

    behavior_data=load(fullfile(handlesexperimentlist{ge},score_file));
    behavior_data=update_t01s_from_postprocessed(behavior_data);
    behavior_data2=[];
    if(behavior_logic>1)
      behavior_data2=load(fullfile(handlesexperimentlist{ge},score_file2));
      behavior_data2=update_t01s_from_postprocessed(behavior_data2);
    end
    behavior_data3=[];
    if(handles.behaviorvalue3>1)
      behavior_data3=load(fullfile(handlesexperimentlist{ge},score_file3));
      behavior_data3=update_t01s_from_postprocessed(behavior_data3);
    end

    [behavior_data,behavior_data2,behavior_data3,~,sex_data]=...
        cull_short_trajectories(handles,behavior_data,behavior_data2,behavior_data3,[],handles.sexdata{ge});
    num_indi=num_indi+length(behavior_data.allScores.scores);

    if(xoffset==1)
      parfor_tmp=zeros(2,max(behavior_data.allScores.tEnd));
    else
      parfor_tmp=zeros(2,max(behavior_data.allScores.tEnd)-min(behavior_data.allScores.tStart)+1);
    end

    ii2=0;
    parfor_tmp2=cell(1,length(individual));
    for i2 = individual
      ii2=ii2+1;
      if(iscell(i2))  i2=char(i2);  end
      for i=1:length(behavior_data.allScores.t0s)   % individual
        if(isnumeric(i2)&&(i2~=i))  continue;  end

        tmp1 = compute_behavior_logic(behavior_data.allScores, i);

        tmp2=[];
        if(behavior_logic>1)
          tmp2 = compute_behavior_logic(behavior_data2.allScores, i);
        end

        tmp3=[];
        if(handles.behaviorvalue3>1)
          tmp3 = compute_behavior_logic(behavior_data3.allScores, i);
          if(handles.behaviornormalizenot)
            tmp3=~tmp3;
          end
        end

        if(behaviornot)  tmp1=~tmp1;  end

        partition_idx=zeros(2,length(tmp1));
        switch(behavior_logic)
          case(1)
            partition_idx(1,:)=tmp1;
          case(2)
            partition_idx(1,:)=tmp1 & tmp2;
          case(3)
            partition_idx(1,:)=tmp1 & ~tmp2;
          case(4)
            partition_idx(1,:)=tmp1 | tmp2;
          case(5)
            partition_idx(1,:)=tmp1 | ~tmp2;
        end
        if(handles.behaviorvalue3==1)
          partition_idx(2,behavior_data.allScores.tStart(i):behavior_data.allScores.tEnd(i))=1;
        else
          partition_idx(2,:)=tmp3;
        end

        [ones(1,behavior_data.allScores.tEnd(i))];
        switch(i2)
          case('M')
            [ones(1,behavior_data.allScores.tStart(i)-1) sex_data{i}];
          case('F')
            [ones(1,behavior_data.allScores.tStart(i)-1) ~sex_data{i}];
        end
        partition_idx(1,:) = partition_idx(1,:) & ans;
        partition_idx(2,:) = partition_idx(2,:) & ans;
        
        idx1 = find(partition_idx(1,:));
        idx2 = find(partition_idx(2,:));
        switch(xoffset)
          case(1)
          case(2)
            idx1 = idx1 - behavior_data.allScores.tStart(i) + 1;
            idx2 = idx2 - behavior_data.allScores.tStart(i) + 1;
          case(3)
            idx1 = idx1 - min(behavior_data.allScores.tStart) + 1;
            idx2 = idx2 - min(behavior_data.allScores.tStart) + 1;
        end
        parfor_tmp(1,idx1)=parfor_tmp(1,idx1)+1;
        parfor_tmp(2,idx2)=parfor_tmp(2,idx2)+1;
      end
      parfor_tmp2{ii2}(1,:)=parfor_tmp(1,:)./parfor_tmp(2,:);
      find(parfor_tmp(2,:)==0);  parfor_tmp(2,ans)=nan;
      parfor_tmp2{ii2}(2,:)=[nan(1,floor((convolutionwidth-1)/2)) ...
         conv(parfor_tmp(1,:),ones(1,convolutionwidth),'valid') ./ ...
         conv(parfor_tmp(2,:),ones(1,convolutionwidth),'valid') ...
         nan(1,ceil((convolutionwidth-1)/2))];
    end
    raw_data(gei,:)=cellfun(@(x) x(1,:),parfor_tmp2,'uniformoutput',false);
    behavior_cumulative(gei,:)=cellfun(@(x) x(2,:),parfor_tmp2,'uniformoutput',false);
  end
  raw_data=reshape(raw_data,1,prod(size(raw_data)));
  behavior_cumulative=reshape(behavior_cumulative,1,prod(size(behavior_cumulative)));

  if(num_indi==0)
    delete(hf);
    set(handles.Status,'string','Ready.','foregroundcolor','g');
    set(handles.figure1,'pointer','arrow');
    uiwait(errordlg('no valid data.  check minimum trajectory length.'));  drawnow;
    return;
  end

  max(cellfun(@(x) size(x,2),behavior_cumulative));
  cellfun(@(x) [x nan(size(x,1),ans-size(x,2))],behavior_cumulative,'uniformoutput',false);
  behavior_cumulative=cat(1,ans{:});

  tstr='';
  time_base=(1:size(behavior_cumulative,2))./handles.fps;
  xstr='time (sec)';
  if(time_base(end)>300)
    time_base=time_base./60;
    xstr='time (min)';
  end
  if(time_base(end)>300)
    time_base=time_base./60;
    xstr='time (hr)';
  end
  if(time_base(end)>120)
    time_base=time_base./24;
    xstr='time (d)';
  end
  ystr='';  if(behaviornot)  ystr='NOT ';  end
  ystr=[ystr char(strrep(handles.behaviorlist(b),'_','-'))];
  switch(handles.behaviorlogic)
    case 2
      ystr=[ystr ' AND '];
    case 3
      ystr=[ystr ' AND NOT '];
    case 4
      ystr=[ystr ' OR '];
    case 5
      ystr=[ystr ' OR NOT '];
  end
  if(handles.behaviorlogic>1)
    ystr=[ystr char(strrep(handles.behaviorlist(handles.behaviorvalue2),'_','-'))];
  end
  ystr=[ystr ' (%)'];

  if(handles.dump2csv)  print_csv_help(fid,handles.type,tstr,xstr,ystr);  end

  ii=0;
  for i = individual
    ii=ii+1;
    if(iscell(i))  i=char(i);  end
    for g=1:length(handles.grouplist)
      color=handles.colors(g,:);

      if(handles.dump2csv)  fprintf(fid,['%% group ' handles.grouplist{g} '\n']);  end

      if ischar(i)
        idx=(cumsum_num_selexp_per_group(g)+1):(cumsum_num_selexp_per_group(g+1));
      else
        find(cumsum_num_exp_per_group<ggee,1,'last');
        if(ans~=g)  continue;  end
        idx=1;
      end
      idx2=idx+(ii-1)*numel(ggee);
      linestyle='-';  if(ii>1)  linestyle='--';  end

      plot_it(ha,time_base,100.*behavior_cumulative(idx2,:),...
          style,centraltendency,dispersion,color,1,linestyle,fid,handlesexperimentlist(idx));
      if (ii==1)
        h(g)=ans;
        hh{g}=handles.grouplist{g};
      end
    end

    if(handles.dump2csv)  fprintf(fid,'\n%% raw data\n');  end
    for g=1:length(handles.grouplist)
      if ischar(i)
        idx=(cumsum_num_selexp_per_group(g)+1):(cumsum_num_selexp_per_group(g+1));
      else
        find(cumsum_num_exp_per_group<ggee,1,'last');
        if(ans~=g)  continue;  end
        idx=1;
      end
      idx2=idx+(ii-1)*numel(ggee);

      if(handles.dump2csv)
        fprintf(fid,['%% group ' handles.grouplist{g} '\n']);
        for e=1:length(idx)
          fprintf(fid,'%% experiment %s\n',handlesexperimentlist{selected_exp(idx(e))});
          print_csv_data(fid,raw_data{idx2(e)});
          fprintf(fid,'\n');
        end
      end
    end
  end

  if(handles.dump2mat)
    find(b==bb);
    raw_data2{ans}=raw_data;
  end

  xlabel(ha,xstr,'interpreter','none');
  ylabel(ha,ystr,'interpreter','none');
  axis(ha,'tight');  zoom(ha,'reset');
end

if(handles.dump2mat)
  raw_data=raw_data2;
  save('most_recent_figure.mat','handles','raw_data');
end

if(iscell(individual))
  h(end+1)=plot(0,0,'k-');   hh{end+1}='males';
  h(end+1)=plot(0,0,'k--');  hh{end+1}='females';
  set(h((end-1):end),'visible','off');
end

idx=find(h>0);
if ~isnumeric(individual)
%  legend(ha,h(idx),[cellfun(@(x) [x ' ' handles.individuallist{handles.individualvalue}],...
%      handles.grouplist,'uniformoutput',false)],'interpreter','none');
  legend(ha,h(idx),hh(idx),'interpreter','none');
%else
%  legend(ha,h(idx),handles.individuallist(handles.individualvalue),'interpreter','none');
end

uicontrol(hf,'style','pushbutton','string','Params','position',[5 5 60 20],...
    'callback',@figure_params_callback);

if(handles.dump2csv)  fclose(fid);  end

guidata(hf,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% --- 
function [bout_lengths sex inter_bout_lengths inter_sex]=...
    calculate_boutstats(behavior_data,behavior_logic,behavior_data2,sexdata,behaviornot)

bout_lengths=cell(1,length(behavior_data.allScores.t0s));
inter_bout_lengths=cell(1,length(behavior_data.allScores.t0s));
sex=cell(1,length(behavior_data.allScores.t0s));
inter_sex=cell(1,length(behavior_data.allScores.t0s));
for i=1:length(behavior_data.allScores.t0s)  % individual
  tmp1 = compute_behavior_logic(behavior_data.allScores, i);
  tmp1 = tmp1(behavior_data.allScores.tStart(i) : behavior_data.allScores.tEnd(i));

  tmp2=[];
  if(behavior_logic>1)
    tmp2 = compute_behavior_logic(behavior_data2.allScores, i);
    tmp2 = tmp2(behavior_data2.allScores.tStart(i) : behavior_data2.allScores.tEnd(i));
  end

  if(behaviornot)  tmp1=~tmp1;  end

  partition_idx=[];
  switch(behavior_logic)
    case(1)
      partition_idx=tmp1;
    case(2)
      partition_idx=tmp1 & tmp2;
    case(3)
      partition_idx=tmp1 & ~tmp2;
    case(4)
      partition_idx=tmp1 | tmp2;
    case(5)
      partition_idx=tmp1 | ~tmp2;
  end

  % inter-bout NOT <behavior> is not quite the same as bout <behavior>

  partition_idx=[0 partition_idx 0];
  start=1+find(~partition_idx(1:(end-1)) &  partition_idx(2:end))-1;
  stop =  find( partition_idx(1:(end-1)) & ~partition_idx(2:end))-1;
  if(length(start)>0)
    bout_lengths{i}=stop-start+1;
    sex{i}=zeros(1,length(bout_lengths{i}));
    for j=1:length(bout_lengths{i})
      sex{i}(j)=sum(sexdata{i}(start(j):stop(j))) > (bout_lengths{i}(j)/2);
    end
    if(length(start)>1)
      inter_bout_lengths{i}=start(2:end)-stop(1:(end-1));
      inter_sex{i}=zeros(1,length(inter_bout_lengths{i}));
      for j=1:length(inter_bout_lengths{i})
        inter_sex{i}(j)=sum(sexdata{i}(stop(j):start(j+1))) > (inter_bout_lengths{i}(j)/2);
      end
    end
  end
end


% ---
function handles=bout_stats_plot(handles)

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

handlesexperimentlist=[handles.experimentlist{:}];

cumsum_num_indi_per_exp=[0 cumsum(handles.individuals_feature)];
cumsum_num_exp_per_group=[0 cumsum(cellfun(@length,handles.experimentlist))];
mat2cell(cumsum_num_exp_per_group(1:end-1),1,ones(1,length(cumsum_num_exp_per_group)-1));
cellfun(@(x,y) x+y,handles.experimentvalue,ans,'uniformoutput',false);
selected_exp=[ans{:}];
cumsum_num_selexp_per_group=[0 cumsum(cellfun(@length,handles.experimentvalue))];

%ggee=1:length(handlesexperimentlist);
ggee=selected_exp;
individual=handles.individualidx;
if isnumeric(individual)
  ggee=find(cumsum_num_indi_per_exp<individual,1,'last');
  individual=individual-cumsum_num_indi_per_exp(ggee);
end

fid=nan;  if(handles.dump2csv)  fid=fopen('most_recent_figure.csv','w');  end

handles.type='bout stats';

hf=figure('toolbar','figure');

bb=handles.behaviorvalue;
if(bb==(length(handles.behaviorlist)+1))  bb=1:(bb-1);  end

behavior_logic=handles.behaviorlogic;
score_file2=handles.scorefiles{handles.behaviorvalue2};
%sexdata=handles.sexdata;
behaviornot=handles.behaviornot;

h=[];
table_data={};
for b=bb

  if(length(bb)>1)
    ceil(sqrt(length(bb)));
    ha=subplot(ceil(length(bb)/ans),ans,b,'parent',hf);
  else
    ha=subplot(1,1,1,'parent',hf);
  end
  hold(ha,'on');

  score_file=handles.scorefiles{b};

  num_indi=0;
  collated_data=cell(1,length(ggee));
  parfor gei=1:length(ggee)
  %for gei=1:length(ggee)
    ge = ggee(gei);

    behavior_data=load(fullfile(handlesexperimentlist{ge},score_file));
    behavior_data=update_t01s_from_postprocessed(behavior_data);
    behavior_data2=[];
    if(behavior_logic>1)
      behavior_data2=load(fullfile(handlesexperimentlist{ge},score_file2));
      behavior_data2=update_t01s_from_postprocessed(behavior_data2);
    end

    [behavior_data,behavior_data2,~,~,sex_data]=...
        cull_short_trajectories(handles,behavior_data,behavior_data2,[],[],handles.sexdata{ge});
    num_indi=num_indi+length(behavior_data.allScores.scores);

    [bout_lengths sex inter_bout_lengths inter_sex]=...
        calculate_boutstats(behavior_data,behavior_logic,behavior_data2,sex_data,behaviornot);
    collated_data{gei}={bout_lengths sex inter_bout_lengths inter_sex};
  end

  if(num_indi==0)
    delete(hf);
    set(handles.Status,'string','Ready.','foregroundcolor','g');
    set(handles.figure1,'pointer','arrow');
    uiwait(errordlg('no valid data.  check minimum trajectory length.'));  drawnow;
    return;
  end

  tstr='';  if(behaviornot)  tstr='NOT ';  end
  tstr=[tstr char(strrep(handles.behaviorlist(b),'_','-'))];
  switch(handles.behaviorlogic)
    case 2
      tstr=[tstr ' AND '];
    case 3
      tstr=[tstr ' AND NOT '];
    case 4
      tstr=[tstr ' OR '];
    case 5
      tstr=[tstr ' OR NOT'];
  end
  if(handles.behaviorlogic>1)
    tstr=[tstr char(strrep(handles.behaviorlist(handles.behaviorvalue2),'_','-'))];
  end
  ystr='bout length (sec)';
  if(handles.boutstats_style2==2)  ystr=['inter-' ystr];  end
  xstr='group';

  if(handles.dump2csv)  print_csv_help(fid,handles.type,tstr,xstr,ystr);  end

  idx=cellfun(@isempty,collated_data);
  collated_data=collated_data(~idx);

  exp_separators=[];  maxy=0;  k=[];  m=0;  ii=0;  table_data{end+1}=[];
  for i = individual
    if(iscell(i))  i=char(i);  end
    idx=handles.boutstats_style2*2-1;
    switch(i)
      case 'A'
        length_data=cellfun(@(x) x{idx},collated_data,'uniformoutput',false);
      case {'M'}
        for ge=1:length(collated_data)
          length_data{ge}=cellfun(@(x,y) x(y==1),...
              collated_data{ge}{idx},collated_data{ge}{idx+1},'uniformoutput',false);
        end
      case {'F'}
        for ge=1:length(collated_data)
          length_data{ge}=cellfun(@(x,y) x(y==0),...
              collated_data{ge}{idx},collated_data{ge}{idx+1},'uniformoutput',false);
        end
      otherwise
        length_data=cellfun(@(x) x{idx}(individual),collated_data,'uniformoutput',false);
    end

    if(handles.dump2csv)
      fprintf(fid,'\n%% raw data\n');
      for g=1:length(handles.grouplist)
        if ~isnumeric(individual)
          idx=(cumsum_num_selexp_per_group(g)+1):(cumsum_num_selexp_per_group(g+1));
        else
          find(cumsum_num_exp_per_group<ggee,1,'last');
          if(ans~=g)  continue;  end
          idx=1;
        end
        fprintf(fid,['%% group ' handles.grouplist{g} '\n']);
        for e=1:length(idx)
          fprintf(fid,'%% experiment %s\n',handlesexperimentlist{selected_exp(idx(e))});
          for i2=1:length(length_data{idx(e)})
            fprintf(fid,'%% individual %d\n',i2);
            print_csv_data(fid,length_data{idx(e)}{i2}./handles.fps);
            fprintf(fid,'\n');
          end
        end
      end
    end

    for g=1:length(handles.grouplist)
      color=handles.colors(g,:);
      ii=ii+1;

      if ischar(i)
        idx=(cumsum_num_selexp_per_group(g)+1):(cumsum_num_selexp_per_group(g+1));
      else
        find(cumsum_num_exp_per_group<ggee,1,'last');
        if(ans~=g)  continue;  end
        idx=1;
      end

      xticklabels{ii}=handles.grouplist{g};

      switch(handles.boutstats_style)
        case 1  % per experiment, error bars
          table_data{end}{ii}=cellfun(@(x) nanmean([x{:}]./handles.fps),length_data(idx));
          [ct(g),dp(g),dn(g)]=...
              calculate_ct_d(table_data{end}{ii},handles.centraltendency,handles.dispersion);
          h{ii}=errorbarplot(ha,ii,ct(g),ct(g)-dn(g),dp(g)-ct(g),color);
        case 2  % per fly, grouped
          cumsum(cellfun(@length,length_data(idx)))';
          exp_separators=[exp_separators; ans+sum(k)];
          table_data{end}{ii}=cellfun(@nanmean,[length_data{idx}])./handles.fps;
          maxy=max([maxy table_data{end}{ii}]);
          h{ii}=bar(ha,(1:length(table_data{end}{ii}))+sum(k),table_data{end}{ii},...
              'barwidth',1,'edgecolor','none');
          set(h{ii},'facecolor',color);
          k(end+1)=length(table_data{end}{ii});
      end
    end
  end

  if(ismember(handles.behaviorbarchart_style,[1 2 4 6]) && (length(individual)==2))
    for g=(length(h)/2+1):length(h)
      findobj([h{g}],'type','patch');
      hatchfill(ans,'single',45,5,handles.colors(g-length(h)/2,:));
    end
  end

  if(handles.dump2csv)
    fprintf(fid,'\n%% summary data\n');
    fprintf(fid,['%% xdata\n']);  fprintf(fid,'%s, ',xticklabels{:});  fprintf(fid,'\n');
  end
  switch(handles.boutstats_style)
    case 1  % per experiment, error bars
      if(handles.dump2csv)
        fprintf(fid,['%% ydata, CT+D\n']);  fprintf(fid,'%g, ',dp);  fprintf(fid,'\n');
        fprintf(fid,['%% ydata, CT-D\n']);  fprintf(fid,'%g, ',dn);  fprintf(fid,'\n');
        fprintf(fid,['%% ydata, CT\n']);    fprintf(fid,'%g, ',ct);  fprintf(fid,'\n');
      end
    case 2  % per fly, grouped
      l=exp_separators(1:2:(end-1));
      r=exp_separators(2:2:end);
      hh=patch(0.5+[l r r l l]',repmat([0 0 maxy*1.05 maxy*1.05 0]',1,floor(length(exp_separators)/2)),...
          [0.95 0.95 0.95],'parent',ha);
      set(hh,'edgecolor','none');
      set(ha,'children',circshift(get(ha,'children'),-1));
      k=round(cumsum(k)-k/2);
      if(handles.dump2csv)
        fprintf(fid,['%% ydata\n']);
        for i=1:length(table_data{end})
          fprintf(fid,'%g, ',[table_data{end}{i}]);
          fprintf(fid,'\n');
        end
      end
  end

  if(handles.dump2mat)
    find(b==bb);
    raw_data{ans}=collated_data;
  end

  if(isempty(k))  k=1:length(length_data);  end
  title(ha,tstr,'interpreter','none');
  ylabel(ha,ystr,'interpreter','none');
  set(ha,'xtick',k,'xticklabel',xticklabels);
  axis(ha,'tight');  vt=axis;
  axisalmosttight([],ha);  vat=axis;
  if(handles.boutstats_style==2)
    axis(ha,[vat(1) vat(2) 0 vt(4)]);
  else
    axis(ha,[vat(1) vat(2) 0 vat(4)]);
  end
  %if(handles.dump2csv)  fprintf(fid,'\n');  end
end

if(handles.dump2mat)
  save('most_recent_figure.mat','handles','raw_data');
end

if(iscell(individual))
  h2(1)=plot(0,0,'k-');   hh2{1}='males';
  h2(2)=plot(0,0,'k--');  hh2{2}='females';
  set(h2,'visible','off');
  legend(ha,h2,hh2,'interpreter','none');
end

uicontrol(hf,'style','pushbutton','string','Params','position',[5 5 60 20],...
    'callback',@figure_params_callback);
if(ischar(individual) && (length(handles.grouplist)>1))
  uicontrol(hf,'style','pushbutton','string','Stats','position',[70 5 50 20],...
      'callback',@figure_stats_callback);
  handles.statistics=calculate_statistics(table_data,handles.behaviorlist(bb),handles.grouplist,...
      fid,handles.pvalue);
end

if(handles.dump2csv)  fclose(fid);  end

guidata(hf,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% --- 
function [during not_during raw_during raw_not_during]=...
    calculate_social(behavior_data,behavior_logic,behavior_data2,feature_data)

if(iscell(feature_data.data{1}))
  vals=unique([feature_data.data{:}]);
  if(length(vals)>2)  error('uhoh');  end
  feature_data.data=cellfun(@(x) strcmp(x,vals{1}),feature_data.data,'uniformoutput',false);
end

during={};  not_during={};
raw_during={};  raw_not_during={};
for i=1:length(behavior_data.allScores.t0s)  % individual
  tmp1=false(1,length(feature_data.data{i}));
  for j=1:length(behavior_data.allScores.t0s{i})
    tmp1((behavior_data.allScores.t0s{i}(j):(behavior_data.allScores.t1s{i}(j)-1))...
        -behavior_data.allScores.tStart(i)+1)=true;
  end
  if(behavior_logic>1)
    tmp2=false(1,length(feature_data.data{i}));
    for j=1:length(behavior_data2.allScores.t0s{i})
      tmp2((behavior_data2.allScores.t0s{i}(j):(behavior_data2.allScores.t1s{i}(j)-1))...
          -behavior_data2.allScores.tStart(i)+1)=true;
    end
  end
  switch(behavior_logic)
    case(1)
      partition_idx=tmp1;
    case(2)
      partition_idx=tmp1 & tmp2;
    case(3)
      partition_idx=tmp1 & ~tmp2;
    case(4)
      partition_idx=tmp1 | tmp2;
  end
  if length(partition_idx)==(length(feature_data.data{i}+1))  % n-1 features
    partition_idx=partition_idx(1:end-1);
  end
  partition_idx=[0 partition_idx 0];
  start=1+find(~partition_idx(1:(end-1)) &  partition_idx(2:end))-1;
  stop =  find( partition_idx(1:(end-1)) & ~partition_idx(2:end))-1;
  during{i}={};  not_during{i}={};
  raw_during{i}={};  raw_not_during{i}={};
  if(length(start)>0)
    for j=1:length(start)
      raw_during{i}{j}=feature_data.data{i}(start(j):stop(j));
      [m,f,c]=mode(raw_during{i}{j});
      during{i}{j}={c{1} f./length(raw_during{i}{j})*100};
    end
    if(length(start)>1)
      for j=1:(length(start)-1)
        raw_not_during{i}{j}=feature_data.data{i}(stop(j):start(j+1));
        [m,f,c]=mode(raw_not_during{i}{j});
        not_during{i}{j}={c{1} f./length(raw_not_during{i}{j})*100};
      end
    end
  end
end



% ---
function ret_val=print_stats2(arg)

sprintf('%d,',arg{1});
ret_val=[ans(1:(end-1)) ' (' num2str(arg{2},3) '%)']; %#ok<COLND>



% ---
function [table_data raw_table_data]=plot_social(experiment_value,experiment_list,...
    behavior_value,behavior_list,behavior_logic,behavior_value2,behavior_list2,...
    feature_value,feature_list,perframe_dir,color)

during={};  not_during={};
parfor e=1:length(experiment_value)
%for e=1:length(experiment_value)
  behavior_data=load(fullfile(experiment_list{experiment_value(e)},...
        [behavior_list{behavior_value} '.mat']));
  behavior_data=update_t01s_from_postprocessed(behavior_data);
  if(behavior_logic>1)
    behavior_data2=load(fullfile(experiment_list{experiment_value(e)},...
        [behavior_list2{behavior_value2} '.mat']));
    behavior_data2=update_t01s_from_postprocessed(behavior_data2);
  else
    behavior_data2=[];
  end
  feature_data=load(fullfile(experiment_list{experiment_value(e)},perframe_dir,...
        [feature_list{feature_value} '.mat']));

  [during{e} not_during{e} raw_during{e} raw_not_during{e}]=...
      calculate_social(behavior_data,behavior_logic,behavior_data2,feature_data);
end

k=1;
table_data={};
raw_table_data={};
for e=1:length(experiment_value)
  for i=1:length(during{e})
    table_data{k,1}=['Exp #' color num2str(experiment_value(e)) ', Indi #' num2str(i)];
    if ~isempty(during{e}{i})
      cellfun(@(x) x{1},during{e}{i},'uniformoutput',false);
      raw_table_data{k,2}=cat(1,ans{:});
      [m,f,c]=mode(raw_table_data{k,2});
      sprintf('%d,',c{1});
      table_data{k,2}=[ans(1:(end-1)) ' (' num2str(f/length(raw_table_data{k,2})*100,3) '%)'];
    end
    if ~isempty(not_during{e}{i})
      cellfun(@(x) x{1},not_during{e}{i},'uniformoutput',false);
      raw_table_data{k,3}=cat(1,ans{:});
      [m,f,c]=mode(raw_table_data{k,3});
      sprintf('%d,',c{1});
      table_data{k,3}=[ans(1:(end-1)) ' (' num2str(f/length(raw_table_data{k,3})*100,3) '%)'];
    end
    if ~isempty(during{e}{i})
      tmp=cellfun(@print_stats2,during{e}{i},'uniformoutput',false);
      [table_data{k,4:2:(2+2*length(during{e}{i}))}]=deal(tmp{:});
      [raw_table_data{k,4:2:(2+2*length(during{e}{i}))}]=deal(raw_during{e}{i}{:});
    end
    if ~isempty(not_during{e}{i})
      tmp=cellfun(@print_stats2,not_during{e}{i},'uniformoutput',false);
      [table_data{k,5:2:(3+2*length(not_during{e}{i}))}]=deal(tmp{:});
      [raw_table_data{k,5:2:(3+2*length(not_during{e}{i}))}]=deal(raw_not_during{e}{i}{:});
    end
    k=k+1;
  end
end


% --- Executes on button press in SocialStats.
function SocialStats_Callback(hObject, eventdata, handles)
% hObject    handle to SocialStats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

experiment_value=get(handles.ExperimentList,'Value');
experiment_list=get(handles.ExperimentList,'String');
experiment_value2=get(handles.ExperimentList2,'Value');
experiment_list2=get(handles.ExperimentList2,'String');
behavior_value=get(handles.BehaviorList,'Value');
behavior_list=get(handles.BehaviorList,'String');
behavior_logic=get(handles.BehaviorLogic,'Value');
behavior_value2=get(handles.BehaviorList2,'Value');
behavior_list2=get(handles.BehaviorList2,'String');
feature_value=get(handles.FeatureList,'Value');
feature_list=get(handles.FeatureList,'String');

if(~strncmp(feature_list(feature_value),'closestfly_',11))
  feature_value=find(cellfun(@(x) strncmp(x,'closestfly_',11),feature_list),1);
  if(isempty(feature_value))
    error('no closestfly feature found');
  end
  set(handles.FeatureList,'Value',feature_value);
end

%axes(handles.Axes);  cla;  hold on;
figure;  hold on;

if(length(experiment_value)>0)
  [table_data raw_table_data]=plot_social(experiment_value,experiment_list,...
      behavior_value,behavior_list,behavior_logic,behavior_value2,behavior_list2,...
      feature_value,feature_list,handles.perframe_dir,'R');
end
if(length(experiment_value2)>0)
  [table_data2 raw_table_data2]=plot_social(experiment_value2,experiment_list2,...
      behavior_value,behavior_list,behavior_logic,behavior_value2,behavior_list2,...
      feature_value,feature_list,handles.perframe_dir,'B');
end

if((length(experiment_value)>0) && (length(experiment_value2)>0))
  tmp(1:2:2*size(table_data,1),1:size(table_data, 2))=table_data;
  tmp(2:2:2*size(table_data2,1),1:size(table_data2,2))=table_data2;
  raw_tmp(1:2:2*size(raw_table_data,1),1:size(raw_table_data, 2))=raw_table_data;
  raw_tmp(2:2:2*size(raw_table_data2,1),1:size(raw_table_data2,2))=raw_table_data2;
  ii=max(size(table_data,2),size(table_data2,2));
  set(handles.Table,'Data',tmp);
  set(handles.Table,'RowStriping','on','BackgroundColor',[1 0.8 0.8; 0.8 0.8 1]);
elseif(length(experiment_value)>0)
  tmp=table_data;
  raw_tmp=raw_table_data;
  ii=size(table_data,2);
  set(handles.Table,'Data',tmp);
  set(handles.Table,'RowStriping','off','BackgroundColor',[1 0.8 0.8]);
else
  tmp=table_data2;
  raw_tmp=raw_table_data2;
  ii=size(table_data2,2);
  set(handles.Table,'Data',tmp);
  set(handles.Table,'RowStriping','off','BackgroundColor',[0.8 0.8 1]);
end

handles.raw_table_data=raw_tmp;
tmp={'Individual' 'Bout Mode' 'Inter Bout Mode'};
cellstr(num2str((1:(((ii-3)/2)+1))','B #%d'))';
[tmp{4:2:ii}]=deal(ans{:});
cellstr(num2str((1: ((ii-3)/2)   )','IB #%d'))';
[tmp{5:2:ii}]=deal(ans{:});
%for i=1:(ii-3)/2
%  tmp{2*i+2}=['B #' num2str(i)];
%  tmp{2*i+3}=['IB #' num2str(i)];
%end
set(handles.Table,'ColumnName',tmp);
set(handles.Table,'ColumnWidth',{100 75 100});

handles.table='social_stats';
guidata(hObject,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% --- 
function CellSelectionCallback(hObject,eventdata)

if(size(eventdata.Indices,1)==0)  return;  end

handles2=guidata(gcf);
handles=guidata(handles2.figure1);

if(strcmp(handles2.table,'social_stats'))
  if(eventdata.Indices(end,2)==1)  return;  end

  figure;  hold on;
  tmp={};
    experiment_value=get(handles.ExperimentList,'Value');
    experiment_value2=get(handles.ExperimentList2,'Value');
    if((length(experiment_value)>0) && (length(experiment_value2)>0))
      if(mod(eventdata.Indices(end,1),2))
        color='r';
      else
        color='b';
      end
    elseif(length(experiment_value)>0)
      color='r';
    else
      color='b';
    end

    data=handles.raw_table_data{eventdata.Indices(end,1),eventdata.Indices(end,2)};
    tmp=unique(handles.raw_table_data{eventdata.Indices(end,1),2});
    hist(data,min(tmp):max(tmp));

  xlabel('closest fly (#)');
  axis tight;
end

if((~strcmp(handles2.table,'dprime')) && (~strcmp(handles2.table,'zscore')))  return;  end

if(((strcmp(handles2.table,'dprime')) && ismember(eventdata.Indices(end,2),[1 3])) || ...
   ((strcmp(handles2.table,'zscore')) && (eventdata.Indices(end,2)==1)))
  group=handles2.table_data(eventdata.Indices(end,1),eventdata.Indices(end,2));
  switch(group)
    case(-1)
      answer=questdlg('____ __ during / all-frames comparisons from table?','','Remove all','Retain only','cancel','cancel');
    case(-2)
      answer=questdlg('____ __ during / not-during comparisons from table?','','Remove all','Retain only','cancel','cancel');
    otherwise
      answer=questdlg(['____ __ ' handles.grouplist{group} ' groups from table?'],'','Remove all','Retain only','cancel','cancel');
  end
  if(ismember(answer,{'cancel',''}))  return;  end
  tmp=get(handles2.Table,'Data');
  if(strcmp(handles2.table,'dprime'))
    if(strcmp(answer,'Remove all'))
      idx=find((handles2.table_data(:,1)~=group)&(handles2.table_data(:,3)~=group));
    else
      idx=find((handles2.table_data(:,1)==group)|(handles2.table_data(:,3)==group));
    end
  else
    if(strcmp(answer,'Remove all'))
      idx=find(handles2.table_data(:,1)~=group);
    else
      idx=find(handles2.table_data(:,1)==group);
    end
  end
  set(handles2.Table,'Data',tmp(idx,:));
  handles2.table_data=handles2.table_data(idx,:);
end

if(((strcmp(handles2.table,'dprime')) && ismember(eventdata.Indices(end,2),[2 4])) || ...
   ((strcmp(handles2.table,'zscore')) && (eventdata.Indices(end,2)==4)))
  thresh=handles2.table_data(eventdata.Indices(end,1),eventdata.Indices(end,2));
  answer=questdlg(['Remove all rows for which the sample size is ____ than or equal to ' num2str(thresh) ...
      ' from table?'],'','less','greater','cancel','cancel');
  if(ismember(answer,{'cancel',''}))  return;  end
  tmp=get(handles2.Table,'Data');
  if(strcmp(handles2.table,'dprime'))
    if(strcmp(answer,'less'))
      idx=find((handles2.table_data(:,2)>thresh)&(handles2.table_data(:,4)>thresh));
    else
      idx=find((handles2.table_data(:,2)<thresh)&(handles2.table_data(:,4)<thresh));
    end
  else
    if(strcmp(answer,'less'))
      idx=find(handles2.table_data(:,4)>thresh);
    else
      idx=find(handles2.table_data(:,4)<thresh);
    end
  end
  set(handles2.Table,'Data',tmp(idx,:));
  handles2.table_data=handles2.table_data(idx,:);
end

if(((strcmp(handles2.table,'dprime')) && (eventdata.Indices(end,2)==5)) || ...
   ((strcmp(handles2.table,'zscore')) && (eventdata.Indices(end,2)==2)))
  not_txt='';
  if(handles2.table_data(eventdata.Indices(end,1),eventdata.Indices(end,2))<0)
    not_txt='not ';
  end
  beh='all frames';
  handles2.table_data(eventdata.Indices(end,1),eventdata.Indices(end,2));
  if(ans~=0)
    beh=handles.behaviorlist{abs(ans)};
  end
  answer=questdlg(['____ __ ' not_txt beh ...
      ' behaviors from table?'],'','Remove all','Retain only','cancel','cancel');
  if(ismember(answer,{'cancel',''}))  return;  end
  tmp=get(handles2.Table,'Data');
  if(strcmp(answer,'Remove all'))
    idx=find(handles2.table_data(:,eventdata.Indices(end,2)) ~= ...
        handles2.table_data(eventdata.Indices(end,1),eventdata.Indices(end,2)));
  else
    idx=find(handles2.table_data(:,eventdata.Indices(end,2)) == ...
        handles2.table_data(eventdata.Indices(end,1),eventdata.Indices(end,2)));
  end
  set(handles2.Table,'Data',tmp(idx,:));
  handles2.table_data=handles2.table_data(idx,:);
end

if(((strcmp(handles2.table,'dprime')) && (eventdata.Indices(end,2)==6)) || ...
   ((strcmp(handles2.table,'zscore')) && (eventdata.Indices(end,2)==3)))
  answer=questdlg(['____ __ ' ...
      handles.featurelist{handles2.table_data(eventdata.Indices(end,1),eventdata.Indices(end,2))} ...
      ' features from table?'],'','Remove all','Retain only','cancel','cancel');
  if(ismember(answer,{'cancel',''}))  return;  end
  tmp=get(handles2.Table,'Data');
  if(strcmp(answer,'Remove all'))
    idx=find(handles2.table_data(:,eventdata.Indices(end,2)) ~= ...
        handles2.table_data(eventdata.Indices(end,1),eventdata.Indices(end,2)));
  else
    idx=find(handles2.table_data(:,eventdata.Indices(end,2)) == ...
        handles2.table_data(eventdata.Indices(end,1),eventdata.Indices(end,2)));
  end
  set(handles2.Table,'Data',tmp(idx,:));
  handles2.table_data=handles2.table_data(idx,:);
end

if(((strcmp(handles2.table,'dprime')) && (eventdata.Indices(end,2)==7)) || ...
   ((strcmp(handles2.table,'zscore')) && (eventdata.Indices(end,2)==5)))
  handles.analysis='feature_histogram';
  handles.behaviorlogic=1;
  handles.individual=1;
  handles.comparison=1;
  if(strcmp(handles2.table,'dprime'))
    handles.behaviornot=round(0.5*(1+-sign(handles2.table_data(eventdata.Indices(end,1),5))));
    handles.behaviorvalue=max(1,abs(handles2.table_data(eventdata.Indices(end,1),5)));
    handles.featurevalue=handles2.table_data(eventdata.Indices(end,1),6);
    if(handles2.table_data(eventdata.Indices(end,1),3)<0)
      handles.comparison=-handles2.table_data(eventdata.Indices(end,1),3);
    end
  end
  if(strcmp(handles2.table,'zscore'))
    handles.behaviornot=round(0.5*(1+-sign(handles2.table_data(eventdata.Indices(end,1),2))));
    handles.behaviorvalue=max(1,abs(handles2.table_data(eventdata.Indices(end,1),2)));
    handles.featurevalue=handles2.table_data(eventdata.Indices(end,1),3);
    %handles.comparison=0;
  end
  button_comparison_set(handles);
  %FeatureHistogram_Callback(hObject, eventdata, handles);
  update_figure(handles);
  handles=feature_histogram_plot(handles);
end

if((strcmp(handles2.table,'dprime')) && (eventdata.Indices(end,2)==8))
  if(isnan(handles2.table_data(eventdata.Indices(end,1),eventdata.Indices(end,2))))
    answer=questdlg(['Remove all rows for which d'' all frames is NaN?'],'','Yes','No','No');
    if(ismember(answer,{'No',''}))  return;  end
    tmp=get(handles2.Table,'Data');
    idx=find(~isnan(handles2.table_data(:,8)));
    set(handles2.Table,'Data',tmp(idx,:));
    handles2.table_data=handles2.table_data(idx,:);
  else
    crit=handles2.table_data(eventdata.Indices(end,1),7) ./ handles2.table_data(eventdata.Indices(end,1),8);
    answer=questdlg(['Remove all rows for which the ratio of d'' to d'' all frames is less than ' num2str(crit) '?'],...
        '','Yes','No','No');
    if(ismember(answer,{'No',''}))  return;  end
    tmp=get(handles2.Table,'Data');
    idx=find((handles2.table_data(:,7)./handles2.table_data(:,8))>=crit);
    set(handles2.Table,'Data',tmp(idx,:));
    handles2.table_data=handles2.table_data(idx,:);
  end
end

update_figure(handles);

guidata(handles.figure1,handles);
guidata(handles2.figure2,handles2);


% ---
function menu_classify_forcecompute_set(handles)

if(handles.classify_forcecompute)
  set(handles.MenuClassifyForceCompute,'Checked','on');
else
  set(handles.MenuClassifyForceCompute,'Checked','off');
end


% --------------------------------------------------------------------
function MenuClassifyForceCompute_Callback(hObject, eventdata, handles)
% hObject    handle to MenuClassifyForceCompute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.classify_forcecompute=~handles.classify_forcecompute;
menu_classify_forcecompute_set(handles);
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuClassify_Callback(hObject, eventdata, handles)
% hObject    handle to MenuClassify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in CloseAll.
function CloseAll_Callback(hObject, eventdata, handles)
% hObject    handle to CloseAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

questdlg('Close all figures?','','Yes','No','No');
if(strcmp(ans,'No'))  return;  end
hh=findobj('type','figure');
for h=1:length(hh)
  if(hh(h)~=handles.figure1)  close(hh(h));  end
end


%[file,path]=uiputfile('*.txt','Select file to save table to');
%if (file==0)  return;  end
%fid=fopen(fullfile(path,file),'w');
%tmp=get(handles.Table,'columnname');
%tmp=transpose(tmp);
%tmp2=get(handles.Table,'data');
%tmp2(2:end+1,:)=tmp2;
%tmp2(1,:)=tmp;
%tmp2=cellfun(@char,tmp2,'uniformoutput',false);
%for i=1:size(tmp2,1)
%  for j=1:size(tmp2,2)
%    fprintf(fid,'%s\t',char(tmp2(i,j)));
%  end
%  fprintf(fid,'\r\n');
%end
%fclose(fid);


% --------------------------------------------------------------------
function MenuHelp_Callback(hObject, eventdata, handles)
% hObject    handle to MenuHelp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MenuHelpDocumentation_Callback(hObject, eventdata, handles)
% hObject    handle to MenuHelpDocumentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

html_file = 'http://jaaba.sourceforge.net/PlottingResults.html';
if isdeployed,
  [stat,msg] = myweb_nocheck(html_file);
  if stat ~= 0,
    errordlg({'Please see documentation at http://jaaba.sourceforge.net'
      'Error opening webpage within MATLAB:'
      msg});
  end
else
  web('-browser',html_file);
end


% --------------------------------------------------------------------
function MenuHelpBugReport_Callback(hObject, eventdata, handles)
% hObject    handle to MenuHelpBugReport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

{...
'Having problems?  Before contacting us first try:'...
''...
'1. Going to the file menu and choosing Update.  Does that fix it?'...
'2. Going to the file menu and choosing Reset.  Add all of your experiments and try again.  Fixed?'...
'3. If you quit out of JAABAPlot and restart, does it still happen?  Does quiting out of Matlab help?  Rebooting your computer?'...
''...
'If it still doesn''t work, please e-mail our Google newsgroup with the following information:'...
''...
'1. Are you using JAABA on Windows, Mac, or Linux?'...
'2. Are you using the compiled version or running it within Matlab?'...
'3. What version of JAABA are you using?'...
'4. What exactly is the problem?  What did you expect JAABAPlot to do?'...
'5. What is the sequence of steps leading up to the problem?  That is, the order of buttons you push.'...
'6. Does it happen every time?'...
'7. JAABAPlot creates a file called most_recent_config.mat.  Please attach it.'...
'8. Please also cut and paste any red error messages on the Matlab command line.  Or a screenshot if you''re using the compiled version.'...
''...
'JAABAPlot is still a work in progress and we appreciate your feedback.'...
};
msgbox(ans,'Bug Report');


% --------------------------------------------------------------------
function MenuFileReset_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFileReset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

questdlg('Reset configuration to default?','','Yes','No','No');
if(strcmp(ans,'No'))  return;  end
handles=initialize(handles);
update_figure(handles);
set(handles.figure1,'pointer','arrow');
guidata(hObject, handles);

% --------------------------------------------------------------------
function MenuFileMultithreading_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFileReset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


numCores = feature('numCores');
c = parcluster;
NumWorkers = c.NumWorkers;
maxJobs = min(numCores,NumWorkers);
prompts = {sprintf('N. threads for computation (btn 1 and %d)',maxJobs)};

while true,
    defaults = {num2str(handles.computation_threads)};
    res = inputdlg(prompts,'Multi-threading Preferences',1,defaults);
    if isempty(res), return, end;
    errs = {};
    
    ischange = false;
    
    computation_threads = str2double(res{1});
    if isnan(computation_threads) || computation_threads < 1 || computation_threads > maxJobs || rem(computation_threads,1) ~= 0,
        errs{end+1} = 'Number of threads devoted to computation must be a positive integer less than or equal to the number of CPU cores';  %#ok<AGROW>
    else
        if(handles.computation_threads ~= computation_threads)
            ischange = true;
        end
    end
    
    if ischange && isempty(errs),
        
        % remove extra computation threads
        if matlabpool('size') > computation_threads,
            set(handles.Status,'string',sprintf('Shrinking matlab pool to %d workers',computation_threads),'foregroundcolor','b');
            set(handles.figure1,'pointer','watch');
            pause(2);
            matlabpool close;
            matlabpool('open',computation_threads);
        end
        
        % add extra computation threads
        if matlabpool('size') < computation_threads,
            set(handles.Status,'string',sprintf('Growing matlab pool to %d workers',computation_threads),'foregroundcolor','b');
            set(handles.figure1,'pointer','watch');
            drawnow;
            if matlabpool('size') > 0,
                matlabpool close;
            end
            pause(1);
            matlabpool('open',computation_threads);
            pause(1);
        end
        
        
        handles.computation_threads = computation_threads;
        %     ClearStatus(handles);
        set(handles.Status,'string','Ready','foregroundcolor','g');
        set(handles.figure1,'pointer','arrow');
        drawnow;
        
    end
    
    if isempty(errs),
        break;
    else
        uiwait(warndlg(errs,'Bad multithreading options'));
    end
    
end
guidata(hObject,handles);
return


% --------------------------------------------------------------------
function MenuFileUpdate_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFileUpdate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

handles=update_experiment_data(handles,true,true,true);
update_figure(handles);
guidata(hObject, handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% --------------------------------------------------------------------
function MenuFileLoad_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFileLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file,path]=uigetfile('*.mat','Select configuration file to open');
if(isnumeric(file) && isnumeric(path) && (file==0) && (path==0))  return;  end
handles=load_configuration_file(fullfile(path,file),hObject,eventdata,handles);
update_figure(handles);
guidata(hObject, handles);


% --------------------------------------------------------------------
function MenuFileSave_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFileSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file,path]=uiputfile('*.mat','Select file to save configuration to');
if(isnumeric(file) && isnumeric(path) && (file==0) && (path==0))  return;  end

filename = fullfile(path,file);
[success,msg] = SaveConfiguration(handles,filename);
if ~success,
  warning('Error saving to file %s: %s',filename,msg);
end


% --------------------------------------------------------------------
function MenuFile_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% ---
function figure_stats_callback(src,evt)

handles=guidata(src);

len=max(cellfun(@(x) length(regexprep(num2str(x),'<[^<]*>','*')),handles.statistics));
hf=figure('menubar','none','toolbar','none','numbertitle','off',...
    'name',['statistics for Figure ' num2str(gcbf)]);
ht=uitable('data',handles.statistics,'columnwidth',num2cell(8*len),...
    'rowname',[],'columnname',[],'rowstriping','off');
extT=get(ht,'extent');
posF=get(hf,'position');
%set(hf,'position',[posF(1) posF(2) 16*(posF(4)<extT(4))+extT(3) posF(4)]);
set(hf,'position',[posF(1) posF(2) min(posF(3),16*(posF(4)<extT(4))+extT(3)) min(posF(4),extT(4))]);
%set(ht,'position',[0 0 16*(posF(4)<extT(4))+extT(3) posF(4)]);
set(ht,'units','normalized','position',[0 0 1 1]);



% ---
function figure_params_callback(src,evt)

handles=guidata(src);

CT={'Mean' 'Median' 'Mode'};
D={'Std. Dev.' 'Std. Err.' '5%-95%' '25%-75%'};
xoffset={'none', 'start', 'min(start)'};

tmp={};

switch handles.type
  case 'feature histogram'
    tmp{end+1}=['style = ' handles.featurehistogram_stylelist{handles.featurehistogram_style}];
    tmp{end+1}=['style2 = ' handles.featurehistogram_stylelist2{handles.featurehistogram_style2}];
    tmp{end+1}=['allframes=' num2str(handles.comparison==1)];
    tmp{end+1}=['notduring=' num2str(handles.comparison==2)];
    tmp{end+1}=['logbinsize=' num2str(handles.logbinsize)];
    tmp{end+1}=['nbins=' num2str(handles.nbins)];
    tmp{end+1}=['central tendency = ' CT{handles.centraltendency}];
    tmp{end+1}=['dispersion = '  D{handles.dispersion}];
  case 'feature time series'
    tmp{end+1}=['style = ' handles.featuretimeseries_stylelist{handles.featuretimeseries_style}];
    tmp{end+1}=['style2 = ' handles.featuretimeseries_stylelist2{handles.featuretimeseries_style2}];
    tmp{end+1}=['subtractmean=' num2str(handles.subtractmean)];
    tmp{end+1}=['conv. width = '  num2str(handles.convolutionwidth) ' sec'];
    tmp{end+1}=['radius=' num2str(handles.windowradius)];
    tmp{end+1}=['x-offset = '  xoffset{handles.xoffset}];
    tmp{end+1}=['central tendency = ' CT{handles.centraltendency}];
    tmp{end+1}=['dispersion = '  D{handles.dispersion}];
  case 'behavior bar chart'
    tmp{end+1}=['style = ' handles.behaviorbarchart_stylelist{handles.behaviorbarchart_style}];
    tmp{end+1}=['central tendency = ' CT{handles.centraltendency}];
    tmp{end+1}=['dispersion = '  D{handles.dispersion}];
    tmp2='';  if(handles.behaviornormalizenot)  tmp2='not ';  end
    tmp3='all frames';  if(handles.behaviorvalue3>1)  tmp3=handles.behaviorlist(handles.behaviorvalue3-1);  end
    tmp{end+1}=['normalize = ' tmp2 char(tmp3)];
  case 'behavior time series'
    tmp{end+1}=['style = ' handles.behaviortimeseries_stylelist{handles.behaviortimeseries_style}];
    tmp{end+1}=['conv. width = '  num2str(handles.convolutionwidth) ' sec'];
    tmp{end+1}=['x-offset = '  xoffset{handles.xoffset}];
    tmp{end+1}=['central tendency = ' CT{handles.centraltendency}];
    tmp{end+1}=['dispersion = '  D{handles.dispersion}];
    tmp2='';  if(handles.behaviornormalizenot)  tmp2='not ';  end
    tmp3='all frames';  if(handles.behaviorvalue3>1)  tmp3=handles.behaviorlist(handles.behaviorvalue3-1);  end
    tmp{end+1}=['normalize = ' tmp2 char(tmp3)];
  case 'bout stats'
    tmp{end+1}=['style = ' handles.boutstats_stylelist{handles.boutstats_style}];
    tmp{end+1}=['style2 = ' handles.boutstats_stylelist2{handles.boutstats_style2}];
    tmp{end+1}=['central tendency = ' CT{handles.centraltendency}];
    tmp{end+1}=['dispersion = '  D{handles.dispersion}];
end
tmp{end+1}=['minimum trajectory length = '  num2str(handles.minimumtrajectorylength)];
tmp{end+1}=['individual = '  num2str(handles.individuallist{handles.individualvalue})];
tmp{end+1}='';

for g=1:length(handles.grouplist)
  if(length(handles.experimentvalue{g})==0)  continue;  end
  tmp{end+1}=['group ' handles.grouplist{g}];
  for e=1:length(handles.experimentvalue{g})
    [~,tmp{end+1},~]=fileparts(handles.experimentlist{g}{handles.experimentvalue{g}(e)});
  end
  tmp{end+1}='';
end

hf=figure('menubar','none','toolbar','none','numbertitle','off',...
    'name',['parameters for Figure ' num2str(gcbf)]);
ht=uitable('data',tmp','columnwidth',num2cell(8*max(cellfun(@length,tmp))),...
    'rowname',[],'columnname',[],'rowstriping','off');
extT=get(ht,'extent');
posF=get(hf,'position');
set(hf,'position',[posF(1) posF(2) min(posF(3),16*(posF(4)<extT(4))+extT(3)) min(posF(4),extT(4))]);
set(ht,'units','normalized','position',[0 0 1 1]);



%BOUT STATISTICS -- not converted to multivariate yet
%
%Create a table of bout lengths (BL) and inter-bout lengths (IBL) using the
%"Bout Stats" button.  As for behavior stats, each row is a behavior;  the
%columns break it down by sex and individual;  the experiment groups are
%striped;  logical operators can be applied; and the Prefs contextual menu
%controls statistics.
%
%Selecting a cell plots a histogram, normalized to unit area, of bout and
%inter-bout lengths.  The "LogY" and "Stats" buttons can be used to scale
%the axis and overlay summary statistics, respectively.  Multiple selected
%cells are overlayed.
%
%
%SOCIAL STATISTICS -- not converted to multivariate yet
%
%Create a table of the closest individual using the "Social Stats" button.
%Choose which behavior to analysis using the pull-down menu in the Behavior
%panel, and which metric to use for closest using the pull-down menu in the
%Feature panel.  In the resulting table, each row is an individual, and the
%columns show the mode and percentage of frames of the closest individual
%for each bout, as well as the overall mode of the per-bout modes.
%
%Selecting a cell in the table plots a histogram of the closest indidividual.



function MinimumTrajectoryLength_Callback(hObject, eventdata, handles)
% hObject    handle to MinimumTrajectoryLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MinimumTrajectoryLength as text
%        str2double(get(hObject,'String')) returns contents of MinimumTrajectoryLength as a double

handles.minimumtrajectorylength=str2num(get(hObject,'String'));
handles.interestingfeaturehistograms_cache=[];
handles.interestingfeaturetimeseries_cache=[];
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function MinimumTrajectoryLength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MinimumTrajectoryLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ConvolutionWidth_Callback(hObject, eventdata, handles)
% hObject    handle to ConvolutionWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ConvolutionWidth as text
%        str2double(get(hObject,'String')) returns contents of ConvolutionWidth as a double

%inputdlg({'Convolution width:'},'',1,{num2str(handles.convolutionwidth)});
%if(isempty(ans))  return;  end
handles.convolutionwidth=str2num(get(hObject,'String'));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function ConvolutionWidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ConvolutionWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PValue_Callback(hObject, eventdata, handles)
% hObject    handle to PValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PValue as text
%        str2double(get(hObject,'String')) returns contents of PValue as a double

%inputdlg({'P value:'},'',1,{num2str(handles.pvalue)});
%if(isempty(ans))  return;  end
handles.pvalue=str2num(get(hObject,'String'));
guidata(hObject,handles);




% --- Executes during object creation, after setting all properties.
function PValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NBins_Callback(hObject, eventdata, handles)
% hObject    handle to NBins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NBins as text
%        str2double(get(hObject,'String')) returns contents of NBins as a double

%handles.nbins=str2num(char(inputdlg({'Number of bins:'},'',1,...
%    {num2str(handles.nbins)})));
handles.nbins=str2num(get(hObject,'String'));
guidata(hObject,handles);




% --- Executes during object creation, after setting all properties.
function NBins_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NBins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function WindowRadius_Callback(hObject, eventdata, handles)
% hObject    handle to WindowRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WindowRadius as text
%        str2double(get(hObject,'String')) returns contents of WindowRadius as a double

handles.windowradius=str2num(get(hObject,'String'));
guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function WindowRadius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WindowRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in StatisticsList.
function StatisticsList_Callback(hObject, eventdata, handles)
% hObject    handle to StatisticsList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns StatisticsList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from StatisticsList

handles.statistic=get(handles.Statistic,'value');


% --- Executes during object creation, after setting all properties.
function StatisticsList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StatisticsList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in StyleList.
function StyleList_Callback(hObject, eventdata, handles)
% hObject    handle to StyleList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns StyleList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from StyleList

get(hObject,'Value');
switch(handles.analysis)
  case 'feature_histogram'
    handles.featurehistogram_style=ans;
  case 'feature_timeseries'
    handles.featuretimeseries_style=ans;
  case 'behavior_barchart'
    handles.behaviorbarchart_style=ans;
  case 'behavior_timeseries'
    handles.behaviortimeseries_style=ans;
  case 'bout_stats'
    handles.boutstats_style=ans;
  case 'interesting_feature_histograms'
end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function StyleList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StyleList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in StyleList2.
function StyleList2_Callback(hObject, eventdata, handles)
% hObject    handle to StyleList2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns StyleList2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from StyleList2

get(hObject,'Value');
switch(handles.analysis)
  case {'feature_histogram','interesting_feature_histograms'}
    handles.featurehistogram_style2=ans;
  case 'feature_timeseries'
    handles.featuretimeseries_style2=ans;
    update_figure(handles);
  case 'behavior_barchart'
  case 'behavior_timeseries'
  case 'bout_stats'
    handles.boutstats_style2=ans;
end
handles.interestingfeaturehistograms_cache=[];
handles.interestingfeaturetimeseries_cache=[];
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function StyleList2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StyleList2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% ---
function button_comparison_set(handles)

set(handles.AllFrames,'Value',0);
set(handles.NotDuring,'Value',0);
switch(handles.comparison)
  case(1), set(handles.AllFrames,'Value',1);
  case(2), set(handles.NotDuring,'Value',1);
end


% --- Executes on button press in AllFrames.
function AllFrames_Callback(hObject, eventdata, handles)
% hObject    handle to AllFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(handles.comparison==1)
  handles.comparison=0;
else
  handles.comparison=1;
end
button_comparison_set(handles);
guidata(hObject,handles);


% --- Executes on button press in NotDuring.
function NotDuring_Callback(hObject, eventdata, handles)
% hObject    handle to NotDuring (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(handles.comparison==2)
  handles.comparison=0;
else
  handles.comparison=2;
end
button_comparison_set(handles);
guidata(hObject,handles);


% --- Executes on button press in OmitInf.
function OmitInf_Callback(hObject, eventdata, handles)
% hObject    handle to OmitInf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.omitinf=~handles.omitinf;
%button_omitinf_set(handles.omitinf);
%handles.interestingfeaturehistograms_cache=[];
guidata(hObject,handles);


% --- Executes on button press in OmitNaN.
function OmitNaN_Callback(hObject, eventdata, handles)
% hObject    handle to OmitNaN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.omitnan=~handles.omitnan;
%button_omitnan_set(handles.omitnan);
%handles.interestingfeaturehistograms_cache=[];
guidata(hObject,handles);


% --- Executes on button press in AbsDPrimeZScore.
function AbsDPrimeZScore_Callback(hObject, eventdata, handles)
% hObject    handle to AbsDPrimeZScore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.absdprimezscore=~handles.absdprimezscore;
%button_absdprime_set(handles.absdprime);
%handles.interestingfeaturehistograms_cache=[];
guidata(hObject,handles);


% ---
function button_comparison2_set(handles)

set(handles.DPrime,'Value',0);
set(handles.ZScore,'Value',0);
switch(handles.comparison2)
  case(0), set(handles.DPrime,'Value',1);
  case(1), set(handles.ZScore,'Value',1);
end


% --- Executes on button press in DPrime.
function DPrime_Callback(hObject, eventdata, handles)
% hObject    handle to DPrime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.comparison2=0;
button_comparison2_set(handles);
guidata(hObject,handles);


% --- Executes on button press in ZScore.
function ZScore_Callback(hObject, eventdata, handles)
% hObject    handle to ZScore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.comparison2=1;
button_comparison2_set(handles);
guidata(hObject,handles);


% --- Executes on button press in DumpToCSV.
function DumpToCSV_Callback(hObject, eventdata, handles)
% hObject    handle to DumpToCSV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.dump2csv=~handles.dump2csv;
guidata(hObject,handles);


% --- Executes on button press in DumpToMAT.
function DumpToMAT_Callback(hObject, eventdata, handles)
% hObject    handle to DumpToMAT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.dump2mat=~handles.dump2mat;
guidata(hObject,handles);


% --- Executes on button press in ForceCompute.
function ForceCompute_Callback(hObject, eventdata, handles)
% hObject    handle to ForceCompute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in SubtractMean.
function SubtractMean_Callback(hObject, eventdata, handles)
% hObject    handle to SubtractMean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.subtractmean=~handles.subtractmean;
%button_subtractmean_set(handles.subtractmean);
guidata(hObject,handles);


% --- Executes on button press in LogBinSize.
function LogBinSize_Callback(hObject, eventdata, handles)
% hObject    handle to LogBinSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.logbinsize=~handles.logbinsize;
%button_logbinsize_set(handles.logbinsize);
guidata(hObject,handles);


% --- Executes on selection change in XOffset.
function XOffset_Callback(hObject, eventdata, handles)
% hObject    handle to XOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns XOffset contents as cell array
%        contents{get(hObject,'Value')} returns selected item from XOffset

handles.xoffset=get(handles.XOffset,'value');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function XOffset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to XOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Statistic.
function CentralTendency_Callback(hObject, eventdata, handles)
% hObject    handle to Statistic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Statistic contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Statistic

handles.centraltendency=get(handles.CentralTendency,'value');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function CentralTendency_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Statistic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Statistic.
function Dispersion_Callback(hObject, eventdata, handles)
% hObject    handle to Statistic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Statistic contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Statistic

handles.dispersion=get(handles.Dispersion,'value');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function Dispersion_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Statistic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
