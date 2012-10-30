function [success,msg,pxpermm] = ...
  ReadArenaParameters_Qtrax(varargin)

success = false;
msg = ''; %#ok<NASGU>

[pxpermm,featfile,roifile] = myparse(varargin,...
  'pxpermm',1,'featfile','','roifile',''); %#ok<ASGLU>

% try to read from the trx file first
if isempty(roifile),
  msg = 'Input roifile not yet set';
  return;
end

try
  roi = load(roifile);
catch ME,
  msg = sprintf('Could not load from roifile %s: %s',roifile,getReport(ME));
  return;
end

if isfield(roi,'scale') && isfield(roi.scale,'x') && isfield(roi.scale,'y'),
  pxpermm = 1/mean([roi.scale.x,roi.scale.y]);
  success = true;
  msg = 'Read pxpermm from roifile';
else
  msg = 'Could not read pxpermm from roifile';
end
