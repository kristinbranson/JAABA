function StoreTrajectories(obj,n,traj,dooverwrite)

if ~exist('dooverwrite','var'),
  dooverwrite = true;
end

traj_fns = {'x','y','theta','a','b','timestamps',...
  'x_mm','y_mm','a_mm','b_mm','theta_mm','dt','sex',...
  'wing_anglel','wing_angler'};

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
    if dooverwrite && exist(filename,'file'),
      try
        delete(filename);
      catch ME,
        warning('Could not delete file %s: %s',filename,getReport(ME));
      end
    end
    if dooverwrite || ~exist(filename,'file'),
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
if isfield(traj(1).arena,'tl')
  count = 1;
  for ndx = flies(:)'
    obj.tl_x(ndx) = traj(count).arena.tl(1);
    obj.tl_y(ndx) = traj(count).arena.tl(2);
    obj.tr_x(ndx) = traj(count).arena.tr(1);
    obj.tr_y(ndx) = traj(count).arena.tr(2);
    obj.bl_x(ndx) = traj(count).arena.bl(1);
    obj.bl_y(ndx) = traj(count).arena.bl(2);
    obj.br_x(ndx) = traj(count).arena.br(1);
    obj.br_y(ndx) = traj(count).arena.br(2);
    count = count+1;
  end
end

% fps
if isfield(traj,'fps'),
  obj.fps(n) = traj(1).fps;
else
  obj.fps(n) = nan;
end

if isfield(traj,'landmark_params')
  obj.landmark_params{n} = traj.landmark_params;
else
  obj.landmark_params{n} = obj.default_landmark_params;
end