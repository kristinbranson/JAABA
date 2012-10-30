function [success,msg,...
  arenatype,arenacenterx,arenacentery,...
  arenaradius,arenawidth,arenaheight,...
  pxpermm] = ...
  ReadArenaParameters_Ctrax(varargin)

success = false;
msg = ''; %#ok<NASGU>

[inmoviefile,...
  arenatype,arenacenterx,arenacentery,...
  arenaradius,arenawidth,arenaheight,...
  pxpermm,intrxfile,annfile] = myparse(varargin,...
  'inmoviefile','',...
  'arenatype','None','arenacenterx',0,'arenacentery',0,...
  'arenaradius',123,'arenawidth',123,'arenaheight',123,...
  'pxpermm',1,'intrxfile','','annfile','');

% try to read from the trx file first
if isempty(intrxfile),
  msg = 'Input trxfile not yet set';
  return;
end

try
  [trx,~,success1] = load_tracks(intrxfile,inmoviefile,...
    'dosave',false,'annname',annfile,'verbose',false);
  if ~success1,
    msg = sprintf('Could not load tracks from trxfile %s',intrxfile);
    return;
  end
catch ME,
  msg = sprintf('Could not load from trxfile %s: %s',intrxfile,getReport(ME));
  return;
end

readarenafrom = '';
readpxpermmfrom = '';
readoff = false;
if isfield(trx,'pxpermm'),
  newpxpermm = nanmean([trx.pxpermm]);
  if ~isnan(newpxpermm) && newpxpermm > 0,
    pxpermm = newpxpermm;
    readpxpermmfrom = 'intrxfile';
  end
  if all(isfield(trx,{'x','x_mm','y','y_mm'})),
    offx = nanmean(([trx.x] - [trx.x_mm])*pxpermm);
    offy = nanmean(([trx.y] - [trx.y_mm])*pxpermm);
    readoff = true;
  end
end

if isfield(trx,'arena'),
  
  % in px
  if all(isfield(trx(1).arena,{'x','y','r'})),
    arenacenterx = trx(1).arena.x;
    arenacentery = trx(1).arena.y;
    arenaradius = trx(1).arena.r;
    arenatype = 'Circle';
    readarenafrom = 'intrxfile';
  elseif ~isempty(readpxpermmfrom) && readoff && ...
      all(isfield(trx(1).arena,{'arena_radius_mm','arena_center_mm_x','arena_center_mm_y'})),
    arenacenterx = trx(1).arena_center_mm_x * pxpermm + offx;
    arenacentery = trx(1).arena_center_mm_y * pxpermm + offy;
    arenaradius = trx(1).arena_radius_mm * pxpermm;
    readarenafrom = 'intrxfile';
  elseif all(isfield(arena,{'tl_px','tr_px','bl_px','br_px'})),
    tl_px = trx(1).arena.tl;
    tr_px = trx(1).arena.tr;
    bl_px = trx(1).arena.bl;
    br_px = trx(1).arena.br;
    arenawidth = mean([abs(tr_px(1)-tl_px(1)),abs(br_px(1)-bl_px(1))]);
    arenaheight = mean([abs(tr_px(2)-br_px(2)),abs(tl_px(2)-bl_px(2))]);
    arenacenterx = mean([tl_px(1),tr_px(1),bl_px(1),br_px(1)]);
    arenacentery = mean([tl_px(2),tr_px(2),bl_px(2),br_px(2)]);
    arenatype = 'Rectangle';
    readarenafrom = 'intrxfile';
  elseif ~isempty(readpxpermmfrom) && readoff && all(isfield(arena,{'tl','tr','bl','br'})),
    tl_px = trx(1).arena.tl(:)'*pxpermm + [offx,offy];
    tr_px = trx(1).arena.tr(:)'*pxpermm + [offx,offy];
    bl_px = trx(1).arena.bl(:)'*pxpermm + [offx,offy];
    br_px = trx(1).arena.br(:)'*pxpermm + [offx,offy];
    arenawidth = mean([abs(tr_px(1)-tl_px(1)),abs(br_px(1)-bl_px(1))]);
    arenaheight = mean([abs(tr_px(2)-br_px(2)),abs(tl_px(2)-bl_px(2))]);
    arenacenterx = mean([tl_px(1),tr_px(1),bl_px(1),br_px(1)]);
    arenacentery = mean([tl_px(2),tr_px(2),bl_px(2),br_px(2)]);
    arenatype = 'Rectangle';
    readarenafrom = 'intrxfile';
  end
  
end

if isempty(readarenafrom) && ~isempty(annfile),
  % try to read from annfile
  [newarenacenterx,newarenacentery,newarenaradius,dosetcirculararena] = ...
    read_ann('arena_center_x','arena_center_y','arena_radius','do_set_circular_arena');
  if dosetcirculararena && newarenaradius > 0 && ~any(isnan([newarenaradius,newarenacenterx,newarenacentery])),
    arenacenterx = newarenacenterx;
    arenacentery = newarenacentery;
    arenaradius = newarenaradius;
    arenatype = 'Circle';
    readarenafrom = 'annfile';
  end
end

success = ~isempty(readarenafrom);
if ~isempty(readarenafrom);
  msg = sprintf('Read arena from %s, ',readarenafrom);
else
  msg = sprintf('Did not read arena, ');
end
if ~isempty(readpxpermmfrom);
  msg = [msg, sprintf('read pxpermm from %s.',readpxpermmfrom)];
else
  msg = [msg, sprintf('Did not read pxpermm.')];
end