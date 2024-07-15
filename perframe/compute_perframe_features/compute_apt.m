function [all_data, units,apt_info] = compute_apt(trx,n,fn)

flies = trx.exp2flies{n};
nflies = numel(flies);
all_data = cell(1,nflies);
all_units = cell(1,nflies);
parts = strsplit(fn,'_');
view_str = regexp(parts{1},'view(\d*)','tokens');
view = str2double(view_str{1}{1});
fn_type = parts{2};
comp_type = parts{3};

apt_info = trx(flies(1)).trkInfo;
if isfield(apt_info,'crop_loc') && ~isempty(apt_info.crop_loc) && numel(apt_info.crop_loc)==4
  crop_loc = apt_info.crop_loc;
  x_center = mean(crop_loc(1:2));
  y_center = mean(crop_loc(3:4));
  global_center = [y_center x_center];
else
  global_center = double(cell2mat(apt_info.params.imsz))/2;
end

if startsWith(fn_type,'social')
  [all_data,units] = compute_apt_social(trx,n,fn,apt_info);
  return;
end

for fndx = 1:nflies
  if view > 1 
    x = trx(flies(fndx)).(sprintf('x_view%d',view));
    y = trx(flies(fndx)).(sprintf('y_view%d',view));
    theta = trx(flies(fndx)).(sprintf('theta_view%d',view));    
  else
    x = trx(flies(fndx)).x;
    y = trx(flies(fndx)).y;
    theta = trx(flies(fndx)).theta;
  end
  dt = trx(flies(fndx)).dt;
  nframes = trx(flies(fndx)).nframes;
  pxpermm = trx.pxpermm;
  n_parts = apt_info.params.n_classes;
  mod_apt_data = trx(flies(fndx)).kpts;
  mod_apt_data = reshape(mod_apt_data,n_parts,[],nframes);
  
  switch fn_type
    case 'global'
      part_num = str2double(parts{end});
      [data,cur_units] = compute_global(mod_apt_data,dt,theta,comp_type,part_num,pxpermm,global_center);
    case 'body'
      part_num = str2double(parts{end});
      [data,cur_units] = compute_body(mod_apt_data,x,y,dt,theta,comp_type,part_num,pxpermm);
    case 'pair'
      part2 = str2double(parts{end});
      part1 = str2double(parts{end-1});
      x1 = permute(mod_apt_data(part1,1,:),[1,3,2]);
      y1 = permute(mod_apt_data(part1,2,:),[1,3,2]);
      x2 = permute(mod_apt_data(part2,1,:),[1,3,2]);
      y2 = permute(mod_apt_data(part2,2,:),[1,3,2]);
      [x1,y1] = convert_to_body(x1,y1,x,y,theta);
      [x2,y2] = convert_to_body(x2,y2,x,y,theta);
      [data,cur_units] = pair_fn(x1,y1,x2,y2,dt,comp_type,pxpermm);
    case 'triad'
      part3 = str2double(parts{end});
      part2 = str2double(parts{end-1});
      part1 = str2double(parts{end-2});
      x1 = permute(mod_apt_data(part1,1,:),[1,3,2]);
      y1 = permute(mod_apt_data(part1,2,:),[1,3,2]);
      x2 = permute(mod_apt_data(part2,1,:),[1,3,2]);
      y2 = permute(mod_apt_data(part2,2,:),[1,3,2]);
      x3 = permute(mod_apt_data(part3,1,:),[1,3,2]);
      y3 = permute(mod_apt_data(part3,2,:),[1,3,2]);
      [data,cur_units] = triad_fn(x1,x2,x3,y1,y2,y3,dt,comp_type,pxpermm);

    otherwise
      error('Undefined computation type type %s',comp_type);
  end
  all_data{fndx} = data;
  all_units{fndx} = cur_units;
end
units = all_units{1};


function [data,units] = compute_global(apt_data, dt, theta, comp_type, part,pxpermm,global_center)
x = permute(apt_data(part,1,:),[1,3,2]);
y = permute(apt_data(part,2,:),[1,3,2]);
relative_fns = {'x','y','dx','dy','sin','cos','dtheta'};
if any(strcmp(relative_fns,comp_type))
  i_x = global_center(2);
  i_y = global_center(1);
  [data,units] = relative(x,y,i_x,i_y,dt,fn,pxpermm);
else
  switch comp_type
    case 'velmag'
     [data,units] = velmag(x, y, dt,pxpermm);
    case 'dphi'
     [data,units] = dphi(x, y, dt, theta);    
    case 'dvelmag'
     [data,units] = dvelmag(x, y, dt,pxpermm);    
    otherwise
      error('Undefined computation type type %s',comp_type);
  end
end

function [b_x,b_y] = convert_to_body(x,y,trx_x,trx_y,theta)
% convert global coordinates into body coordinates
b_x = nan(size(x)); b_y = nan(size(y));
for t = 1:numel(x)
  T = [1,0,0
    0,1,0
    -trx_x(t),-trx_y(t),1];
  R = [cos(theta(t)-pi/2),-sin(theta(t)-pi/2),0
    sin(theta(t)-pi/2),cos(theta(t)-pi/2),0
    0,0,1];
  A = T*R;
  M = [x(t) y(t) 1]*A;
  b_x(t) = M(1); b_y(t) = M(2);
end


function [data,units] = compute_body(apt_data, trx_x, trx_y, dt, theta, comp_type, part,pxpermm)
relative_fns = {'x','y','dx','dy','sin','cos','dtheta'};
x = permute(apt_data(part,1,:),[1,3,2]);
y = permute(apt_data(part,2,:),[1,3,2]);
[x,y] = convert_to_body(x,y,trx_x,trx_y,theta);

switch comp_type
  case 'velmag'
   [data,units] = velmag(x, y, dt, pxpermm);
  case 'dphi'
   [data,units] = dphi(x, y, dt, theta);    
  case 'dvelmag'
   [data,units] = dvelmag(x, y, dt, pxpermm);    
  case relative_fns
    [data,units] = relative(x,y,zeros(size(x)),zeros(size(y)),dt, comp_type, pxpermm);
  case 'distcenter'
    [data,units] = dist_center(x,y,zeros(size(x)),zeros(size(y)), pxpermm);
  case 'ddistcenter'
    [data,units] = ddist_center(x,y,zeros(size(x)),zeros(size(y)),dt, pxpermm);    
  otherwise
    error('Undefined computation type type %s',comp_type);
end

function [data, units] = velmag(x,y,dt,pxpermm)
if numel(x) == 1
  data = 0;
else
  dx = diff(x,1,2);
  dy = diff(y,1,2);
  data = sqrt(dx.^2 + dy.^2)./dt;
end
data = data/pxpermm;
units = parseunits('mm/s');


function [data, units] = dvelmag(x,y,dt,pxpermm)
if numel(x) == 1
  data = [];
elseif numel(x) == 2
  data = 0;
else
  dx = diff(x,1,2);
  dy = diff(y,1,2);
  tmp = sqrt(diff(dx./dt,1,2).^2 + diff(dy./dt,1,2).^2)./dt(2:end);
  data = [0,tmp];
end
data = data/pxpermm;
units = parseunits('mm/s/s');


function [data, units] = dphi(x,y,dt,theta)
if numel(x) > 1
  dy1 = [y(2)-y(1),(y(3:end)-y(1:end-2))/2,y(end)-y(end-1)];
  dx1 = [x(2)-x(1),(x(3:end)-x(1:end-2))/2,x(end)-x(end-1)];
  badidx = dy1 == 0 & dx1 == 0;
  phi = atan2(dy1,dx1);
  phi(badidx) = theta(badidx);
  data = modrange(diff(phi),-pi,pi)./dt;
else
  data = [];
end
units = parseunits('rad/s');


function [data, units] = dist_center(x,y,x_body,y_body,pxpermm)
data = sqrt((x-x_body).^2 + (y-y_body).^2);
data = data/pxpermm;
units = parseunits('mm');


function [data, units] = ddist_center(x,y, x_body, y_body,dt,pxpermm)
if numel(x)>1
  dist = sqrt((x-x_body).^2 + (y-y_body).^2);
  data = diff(dist)./dt;
else
  data = 0;
end
data = data/pxpermm;
units = parseunits('mm/s');


function [data, units] = relative(x1, y1, x2, y2, dt, fn,pxpermm)

switch fn
  case 'x'
    data = x1-x2;
    data = data/pxpermm;
    units = parseunits('mm');
  case 'y'
    data = y1-y2;
    data = data/pxpermm;
    units = parseunits('mm');
  case 'dx'
    data = diff(x1-x2)./dt;
    data = data/pxpermm;
    units = parseunits('mm/s');
  case 'dy'
    data = diff(y1-y2)./dt;
    data = data/pxpermm;
    units = parseunits('mm');
  case 'cos'
    len = sqrt((x1-x2).^2 + (y1-y2).^2);    
    data = (x1-x2)./len;
    data(len==0) = 0;
    units = parseunits('');
  case 'sin'
    len = sqrt((x1-x2).^2 + (y1-y2).^2);
    data = (y1-y2)./len;
    data(len==0) = 0;
    units = parseunits('');
  case 'dtheta'
    if numel(x1) > 2
      theta = atan2(y1-y2, x1-x2);
      data = modrange(diff(theta,1,2),-pi,pi)./dt;
    else
      data = [];
    end
    units = parseunits('rad/s');
end



function [data,units] = pair_fn(x1,y1,x2,y2,dt, fn,pxpermm)

relative_fns = {'x','y','dx','dy','sin','cos','dtheta'};
if any(strcmp(relative_fns,fn))
  [data,units] = relative(x1,y1,x2,y2,dt,fn,pxpermm);
elseif strcmp(fn,'velmag')
  [data,units] = velmag(x1-x2,y1-y2,dt,pxpermm);
elseif strcmp(fn,'dvelmag')
  [data,units] = dvelmag(x1-x2,y1-y2,dt,pxpermm);
elseif strcmp(fn,'dist')
  [data,units] = dist_center(x1,y1,x2,y2,pxpermm);
elseif strcmp(fn,'ddist')
  [data,units] = ddist_center(x1,y1,x2,y2,dt,pxpermm);
elseif strcmp(fn,'dphi')
  theta = atan2(y1-y2,x1-x2);
  [data,units] = dphi(x1-x2, y1-y2, dt, theta);
elseif strcmp(fn,'u1')
  theta = atan2(y1-y2,x1-x2);
  dx = diff(x1,1,2)./dt;
  dy = diff(y1,1,2)./dt;
  phi = [atan2(dy,dx)];
  angle_vel = pi - theta(1:end-1) + phi; 
  % This is the angle between velocity vector and the vector from x2 to x1.
  data = sqrt(dx.^2 + dy.^2).*cos(angle_vel);
  data = data/pxpermm;
  units = parseunits('mm/s');
elseif strcmp(fn,'v1')
  theta = atan2(y1-y2,x1-x2);
  dx = diff(x1,1,2)./dt;
  dy = diff(y1,1,2)./dt;
  phi = [atan2(dy,dx)];
  angle_vel = pi - theta(1:end-1) + phi; 
  % This is the angle between velocity vector and the vector from x2 to x1.
  data = sqrt(dx.^2 + dy.^2).*sin(angle_vel);
  data = data/pxpermm;
  units = parseunits('mm/s');
elseif strcmp(fn,'u2')
  theta = atan2(y2-y1,x2-x1);
  dx = diff(x2,1,2)./dt;
  dy = diff(y2,1,2)./dt;
  phi = [atan2(dy,dx)];
  angle_vel = pi - theta(1:end-1) + phi; 
  % This is the angle between velocity vector and the vector from x2 to x1.
  data = sqrt(dx.^2 + dy.^2).*cos(angle_vel);
  data = data/pxpermm;
  units = parseunits('mm/s');
elseif strcmp(fn,'v2')
  theta = atan2(y2-y1,x2-x1);
  dx = diff(x2,1,2)./dt;
  dy = diff(y2,1,2)./dt;
  phi = atan2(dy,dx);
  angle_vel = pi - theta(1:end-1) + phi; 
  % This is the angle between velocity vector and the vector from x2 to x1.
  data = sqrt(dx.^2 + dy.^2).*sin(angle_vel);
  data = data/pxpermm;
  units = parseunits('mm/s');
elseif strcmp(fn,'areaswept')
  area = zeros(1,numel(x1));
  for ndx = 1:numel(x1)-1
    curx = [x1(ndx),x1(ndx+1),x2(ndx+1),x2(ndx)];
    cury = [y1(ndx), y1(ndx+1), y2(ndx+1), y2(ndx)];
    area(ndx) = polyarea(curx,cury);
  end
  data = area/pxpermm/pxpermm;
  units = parseunits('mm*mm');
else
  error('Unknown function for APT pair %s',fn);
end

function [data,units] = triad_fn(x1,x2,x3, y1,y2, y3, dt,fn,pxpermm)
theta1 = atan2(y1-y2,x1-x2);
theta2 = atan2(y3-y2,x3-x2);
theta = modrange(theta1-theta2,-pi,pi);
switch fn
  case 'cos'
    data = cos(theta);
    units = parseunits('');
  case 'sin'
    data = sin(theta);
    units = parseunits('');
  case 'dangle'
    data = diff(theta)./dt;
    units = parseunits('rad/s');
  case 'area'
    area = zeros(1,numel(x1));
    for ndx = 1:numel(x1)
      area(ndx) = polyarea([x1(ndx),x2(ndx),x3(ndx)],[y1(ndx),y2(ndx),y3(ndx)]);
    end
    data = area;
    data = data/pxpermm/pxpermm;
    units = parseunits('mm*mm');
  case 'darea'
    area = zeros(1,numel(x1));
    for ndx = 1:numel(x1)
      area(ndx) = polyarea([x1(ndx),x2(ndx),x3(ndx)],[y1(ndx),y2(ndx),y3(ndx)]);
    end
    data = diff(area);
    data = data/pxpermm/pxpermm;
    units = parseunits('mm*mm/s');
  case 'dlen'
    len1 = sqrt((x1-x2).^2 + (y1-y2).^2);
    len2 = sqrt((x3-x2).^2 + (y3-y2).^2);
    data = len1-len2;
    data = data/pxpermm;
    units = parseunits('mm');
end

