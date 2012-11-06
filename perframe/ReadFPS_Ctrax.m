function [success,msg,fps] = ...
  ReadFPS_Ctrax(varargin)

success = false;
msg = ''; %#ok<NASGU>

[fps,intrxfile,annfile] = myparse(varargin,...
  'fps',1,'intrxfile','','annfile','');

% try to read from the trx file first
if isempty(intrxfile),
  msg = 'Input trxfile not yet set';
  return;
end

try
  [trx,~,success1,timestamps] = load_tracks(intrxfile,'',...
    'dosave',false,'annname',annfile,'verbose',false);
  if ~success1,
    msg = sprintf('Could not load tracks from trxfile %s',intrxfile);
    return;
  end
  readfps = false;
  if isfield(trx,'fps'),
    fps = [trx.fps];
    if nnz(~isnan(fps)) > 0,
      fps = nanmedian(fps);
      readfps = true;
    end
  end
  if ~readfps,
    if nnz(~isnan(timestamps)) == 0,
      msg = sprintf('No timestamps read in from trxfile %s',intrxfile);
      return;
    end
    fps = 1/nanmedian(diff(timestamps));
  end
catch ME,
  msg = sprintf('Could not load from trxfile %s: %s',intrxfile,getReport(ME));
  return;
end

success = true;
msg = 'Read fps from trxfile';
