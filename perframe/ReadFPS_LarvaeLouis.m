function [success,msg,fps] = ...
  ReadFPS_LarvaeLouis(varargin)

success = false;
msg = ''; %#ok<NASGU>

[fps,indatafile] = myparse(varargin,...
  'fps',1,'indatafile');

% try to read from the trx file first
if isempty(indatafile),
  msg = 'Input data file not yet set';
  return;
end

if ~exist(indatafile,'file'),
  msg = sprintf('Data file %s does not exist',indatafile);
  return;
end

try

  indata.data = [];
  fid = fopen(indatafile,'r');
  if fid < 0,
    error('Could not open file %s for reading',indatafile);
  end
  while true,
    s = fgetl(fid);
    if ~ischar(s),
      break;
    end
    ss = strsplit(s,',');
    ss = str2num(char(ss(1:17))); %#ok<ST2NM>
    indata.data(:,end+1) = ss;
  end
  %indata.data = [(0:size(indata.data,1)-1)',indata.data];
  fclose(fid);
  
  timestamps = indata.data(2,:);
  dt = diff(timestamps);
  fps = 1/median(dt);
  msg = sprintf('Read fps = %f from %s',fps,indatafile);
  success = true;

catch ME,
  msg = sprintf('Error reading timestamps from %s: %s',indatafile,getReport(ME));
  return;
end
