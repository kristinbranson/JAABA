function [success,msg,fps] = ReadFPS_Qtrax(varargin)

success = false;
msg = ''; %#ok<NASGU>

[fps,featfile,roifile] = myparse(varargin,...
  'fps',30,'featfile','','roifile','');  %#ok<NASGU>

% try to read from the trx file first
if isempty(featfile),
  msg = 'Input featfile not yet set';
  return;
end

try
  feat = load(featfile);
catch ME,
  msg = sprintf('Could not load from featfile %s: %s',featfile,getReport(ME));
  return;
end

if isfield(feat,'fly_feat') && isfield(feat.fly_feat,'frame') && isfield(feat.fly_feat,'time'),
  fps = nanmedian(diff(tmp2.fly_feat.frame)./diff([tmp2.fly_feat.time]));
  success = true;
  msg = 'Read fps from featfile';
else
  msg = 'Could not read fps from featfile';
end
