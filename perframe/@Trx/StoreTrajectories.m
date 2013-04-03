function StoreTrajectories(obj,n,traj,dooverwrite)

if ~exist('dooverwrite','var'),
  dooverwrite = true;
end

traj_fns = {'x','y','theta','a','b','timestamps','area',...
  'x_mm','y_mm','a_mm','b_mm','theta_mm','dt','sex','area_mm',...
  'wing_anglel','wing_angler',...
  'xspine','yspine','xspine_mm','yspine_mm'};

for i = 1:numel(traj_fns),
  fn = traj_fns{i};
  if ~isfield(traj,fn), continue; end
  obj.SetPerFrameData(fn,{traj.(fn)},n);
  
  filename = obj.GetPerFrameFile(fn,n);
  data = {traj.(fn)}; %#ok<*NASGU>
  switch fn,
    case {'x','y','a','b'},
      units = parseunits('px');
    case {'x_mm','y_mm','a_mm','b_mm'},
      units = parseunits('mm');
    case {'theta','theta_mm','wing_anglel','wing_angler'},
      units = parseunits('rad');
    case 'dt',
      units = parseunits('s');
    case 'timestamps',
      units = parseunits('days');
    case 'sex',
      units = parseunits('unit');
  end
  try
    
    fileexists = exist(filename,'file');
    if ~fileexists && isunix,
      [res,link] = unix(sprintf('readlink %s',filename));
      if ~res && ~isempty(link),
        warning('Deleting broken soft link from %s to %s.\n',filename,link);
        unix(sprintf('rm %s',filename));
      end
    end
    
    if dooverwrite && fileexists,
      try
        delete(filename);
      catch ME,
        warning('Could not delete file %s: %s',filename,getReport(ME));
      end
    end
    if dooverwrite || ~fileexists,
      save(filename,'data','units');
    end
  catch %#ok<CTCH>
    if ~exist(filename,'file'),
      error('Could not save to file %s',filename);
    else
      fprintf('Could not save to file %s, but file exists, so skipping\n',filename);
    end
  end
end
% conversion from pixels to mm
if isfield(traj,'pxpermm'),
  obj.pxpermm(n) = traj(1).pxpermm;
else
  obj.pxpermm(n) = nan;
end
flies = obj.exp2flies{n};
obj.firstframes(flies) = [traj.firstframe];
obj.endframes(flies) = [traj.endframe];
obj.nframes(flies) = [traj.nframes];

% Adding the arena parameters if they exist for rectangular arena.
setlandmarkparams = false;
if isfield(traj(1),'arena') 
  
  if isfield(traj(1).arena,'tl')
    count = 1;
    for ndx = flies(:)'
      obj.landmark_params{n}.tl_x(ndx) = traj(count).arena.tl(1);
      obj.landmark_params{n}.tl_y(ndx) = traj(count).arena.tl(2);
      obj.landmark_params{n}.tr_x(ndx) = traj(count).arena.tr(1);
      obj.landmark_params{n}.tr_y(ndx) = traj(count).arena.tr(2);
      obj.landmark_params{n}.bl_x(ndx) = traj(count).arena.bl(1);
      obj.landmark_params{n}.bl_y(ndx) = traj(count).arena.bl(2);
      obj.landmark_params{n}.br_x(ndx) = traj(count).arena.br(1);
      obj.landmark_params{n}.br_y(ndx) = traj(count).arena.br(2);
      count = count+1;
    end
    setlandmarkparams = true;
  elseif isfield(traj(1).arena,'arena_radius_mm')
    count = 1;
    for ndx = flies(:)'
      obj.landmark_params{n}.arena_radius_mm(ndx) = traj(count).arena.arena_radius_mm;
      obj.landmark_params{n}.arena_center_mm_x(ndx) = traj(count).arena.arena_center_mm_x;
      obj.landmark_params{n}.arena_center_mm_y(ndx) = traj(count).arena.arena_center_mm_y;
      count = count+1;
    end
    setlandmarkparams = true;
  end
end
if ~setlandmarkparams,
  fns = fieldnames(obj.default_landmark_params);
  obj.landmark_params{n} = obj.default_landmark_params;
  for i = 1:numel(fns),
    obj.landmark_params{n}.(fns{i}) = repmat(obj.landmark_params{n}.(fns{i}),[1,numel(flies)]);
  end
end

% fps
if isfield(traj,'fps'),
  obj.fps(n) = traj(1).fps;
else
  obj.fps(n) = nan;
end

% if isfield(traj,'landmark_params')
%   obj.landmark_params{n} = traj.landmark_params;
% else
%   obj.landmark_params{n} = obj.default_landmark_params;
% end
% roi
if isfield(traj,'roi'),
  for i = 1:numel(flies),
    obj.roi(flies(i)) = traj(i).roi;
  end
else
  obj.roi(flies) = 1;
end
