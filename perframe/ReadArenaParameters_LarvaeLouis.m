function [success,msg,pxpermm] = ReadArenaParameters_LarvaeLouis(varargin) 

success = false;
msg = '';

[inconfigfile,pxpermm] = myparse_nocheck(varargin,'inconfigfile','','pxpermm',nan);

if ~exist(inconfigfile,'file'),
  msg = sprintf('Config file %s does not exist',inconfigfile);
  return;
end

% read in relevant parameters from config file
fid = fopen(inconfigfile,'r');
if fid < 0,
  msg = sprintf('Could not open file %s for reading',inconfigfile);
  return;
end

readpxpermm = false;
while true,
  s = fgetl(fid);
  if ~ischar(s),
    break;
  end
  s = strtrim(s);
  if isempty(s),
    continue;
  end
  m = regexp(s,'^Camera Calibration \(um per pixel\): (.*)$','tokens','once');
  if ~isempty(m),
    pxpermm = 1000/str2double(m{1});
    readpxpermm = true;
    break;
  end
end
fclose(fid);
if ~readpxpermm,
  msg = sprintf('Could not read px per mm from %s',inconfigfile);
  return;
end
msg = sprintf('Read pxpermm = %f from %s',pxpermm,inconfigfile);
success = true;