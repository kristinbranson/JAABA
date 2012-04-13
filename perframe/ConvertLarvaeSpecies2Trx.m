function [outexpdir] = ConvertLarvaeSpecies2Trx(inexpdir,rootoutputdir,varargin)

%% set parameters

% default parameters
trxfilestr = 'trx.mat';
moviefilestr = 'movie.wmv';
makesoftlink = true;
intrxfilestr = 'trx.mat';
inmovieformat = 'wmv';
arena_radius_mm = 60;
swapxy = false;

[trxfilestr,moviefilestr,makesoftlink,intrxfilestr,inmovieformat,arena_radius_mm,swapxy] = ...
  myparse(varargin,...
  'trxfilestr',trxfilestr,...
  'moviefilestr',moviefilestr,...
  'makesoftlink',makesoftlink,...
  'intrxfilestr',intrxfilestr,...
  'inmovieformat',inmovieformat,...
  'arena_radius_mm',arena_radius_mm,...
  'swapxy',swapxy); 

%% get paths to input data

[rootdatadir,inexperiment_name] = fileparts(inexpdir); %#ok<ASGLU>

files = dir(fullfile(inexpdir,['*.',inmovieformat]));
if isempty(files),
  error('No %s movie found in %s',inmovieformat,inexpdir);
end
if numel(files) > 1,
  warning('More than one %s movie found in %s, choosing %s',inmovieformat,inexpdir,files(1).name);
  files = files(1);
end
inmoviefile = fullfile(inexpdir,files(1).name);

intrxfile = fullfile(inexpdir,intrxfilestr);
if ~exist(intrxfile,'file'),
  error('Could not find trx file %s',intrxfile);
end

%% create output directory

% make sure the root output directory exists
if ~exist(rootoutputdir,'dir'),
  mkdir(rootoutputdir);
end

% create the experiment directory
outexperiment_name = inexperiment_name;
outexpdir = fullfile(rootoutputdir,outexperiment_name);
if ~exist(outexpdir,'dir'),
  mkdir(outexpdir);
end

%% load in trx data
% 
% load(intrxfile,'trx');
% nlarvae = numel(trx); %#ok<NODEF>

%% get arena location
% 
% arena = full(trx(1).arena);
% stats = regionprops(arena,'centroid','majoraxislength');
% arena_center_x = stats.Centroid(1);
% arena_center_y = stats.Centroid(2);
% arena_radius = stats.MajorAxisLength/2;

%% swap x and y if nec
% 
% if swapxy,
%   for i = 1:nlarvae,
%     tmp = trx(i).x;
%     trx(i).x = trx(i).y; %#ok<*AGROW>
%     trx(i).y = tmp; 
%     isdata = ~cellfun(@isempty,trx(i).spine);
%     trx(i).spine(isdata) = cellfun(@(x) x(:,[2,1]),trx(i).spine(isdata),'UniformOutput',false); 
%     trx(i).contour = cellfun(@(x) x(:,[2,1]),trx(i).contour,'UniformOutput',false); 
%   end
% end

%% remove extra stuff at the end of fields

% fns = {'x','y','a','b','theta','area'};
% for i = 1:numel(trx),
%   isdiscrepency = false;
%   for j = 1:numel(fns),
%     fn = fns{j};
%     if numel(trx(i).(fn)) ~= trx(i).nframes,
%       if ~isdiscrepency,
%         fprintf('Warning: for larva %d, size of fields and nframes = %d, endframe=%d - firstframe=%d -> %d do not match: ',...
%           i,trx(i).nframes,trx(i).lastframe,trx(i).firstframe,trx(i).lastframe-trx(i).firstframe+1);
%         isdiscrepency = true;
%       end
%       fprintf('|%s|=%d ',fn,numel(trx(i).(fn)));
%       trx(i).(fn) = trx(i).(fn)(1:trx(i).nframes); %#ok<AGROW>
%     end
%   end
%   if isdiscrepency,
%     fprintf('\n');
%   end
% end

%% convert to mm

% pxpermm = arena_radius / arena_radius_mm;
% 
% for i = 1:nlarvae,
%   trx(i).x_mm = (trx(i).x - arena_center_x) / pxpermm;
%   trx(i).y_mm = (trx(i).y - arena_center_y) / pxpermm;
%   trx(i).a_mm = trx(i).a / pxpermm;
%   trx(i).b_mm = trx(i).b / pxpermm;
%   trx(i).theta_mm = trx(i).theta;
%   trx(i).pxpermm = pxpermm;
%   % convert spines to matrix
%   if i == 1,
%     nspinepts = max(cellfun(@(x) size(x,1),trx(i).spine));
%   end
%   trx(i).xspine = nan(nspinepts,trx(i).nframes);
%   trx(i).yspine = nan(nspinepts,trx(i).nframes);
%   isdata = ~cellfun(@isempty,trx(i).spine);
%   trx(i).xspine(:,isdata) = cell2mat(cellfun(@(x) x(:,1),trx(i).spine(isdata),'UniformOutput',false));
%   trx(i).yspine(:,isdata) = cell2mat(cellfun(@(x) x(:,2),trx(i).spine(isdata),'UniformOutput',false));
%   trx(i).xspine_mm = (trx(i).xspine - arena_center_x) / pxpermm;
%   trx(i).yspine_mm = (trx(i).yspine - arena_center_y) / pxpermm;
%   
%   trx(i).xcontour = cellfun(@(x) x(:,1),trx(i).contour,'UniformOutput',false);
%   trx(i).ycontour = cellfun(@(x) x(:,2),trx(i).contour,'UniformOutput',false);
%   trx(i).xcontour_mm = cellfun(@(x) (x - arena_center_x)/pxpermm,trx(i).xcontour,'UniformOutput',false);
%   trx(i).ycontour_mm = cellfun(@(x) (x - arena_center_y)/pxpermm,trx(i).ycontour,'UniformOutput',false);
%   
%   trx(i).area_mm = trx(i).area / pxpermm^2;
%   
% end
% 
% %% create dt, compress arena, rename lastframe as endframe
% 
% for i = 1:nlarvae,
%   trx(i).dt = repmat(1/trx(i).fps,[1,trx(i).nframes-1]);
%   trx(i).arena = struct('arena_center_mm_x',0,'arena_center_mm_y',0,'arena_radius_mm',arena_radius_mm,...
%     'arena_center_px_x',arena_center_x,'arena_center_px_y',arena_center_y,'arena_radius_px',arena_radius);
%   trx(i).endframe = trx(i).lastframe;
% end
% 
% %% remove redundant fields
% 
% trx = rmfield(trx,{'spine','contour','lastframe'}); %#ok<NASGU>
% 
% %% create the trx file
% 
% outtrxfile = fullfile(outexpdir,trxfilestr);
% save(outtrxfile,'trx');

%% copy over the movie file

outmoviefile = fullfile(outexpdir,moviefilestr);
if isunix && makesoftlink,
  if exist(outmoviefile,'file'),
    delete(outmoviefile);
  end
  cmd = sprintf('ln -s %s %s',inmoviefile,outmoviefile);
  unix(cmd);
else
  [success,msg] = copyfile(inmoviefile,outmoviefile);
  if ~success,
    error('Error copying file %s to %s: %s',inmoviefile,outmoviefile,msg);
  end
end

