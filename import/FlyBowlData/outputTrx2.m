function outputTrx2(trx,idx,moviename,matname,savename,t0,t1,fns,flies)

% if flies not input, then we assume the trajectories each correspond to a
% single fly
if ~exist('flies','var'),
  flies = idx;
end

version = 0.1;
fid = fopen(savename,'wb');
if fid < 0,
  error('Could not open file %s for writing',savename);
end

if ~exist('t0','var')
  t0 = trx(idx).firstframe;
end
if ~exist('t1','var'),
  t1 = trx(idx).endframe;
end
% if ~exist('fnsignore','var'),
%   fnsignore = {'timestamps'};
% end
if ~exist('fns','var'),
  fns = trx.PerFrameFieldNames;
end

% format:
% version (double)
% moviename_length (int)
% moviename (char*)
% matname_length (int)
% matname (char*)
% nflies (int)
% fly ids (int*)
% firstframe (int)
% lastframe (int)
% sex (M/F) (char)
% fps (double)
% nfields (int)
% fieldname1_length (int)
% fieldname1 (char*)
% units_numerator1_length (int)
% units_numerator1 (char*)
% units_denominator1_length (int)
% units_denominator1 (char*)
% trx.fieldname1(firstframe) (double)
% ...
% trx.fieldname1(endframe) (double)
% ...
% fieldnamen_length (int)
% fieldnamen (char*)
% units_numeratorn_length (int)
% units_numeratorn (char*)
% units_denominatorn_length (int)
% units_denominatorn (char*)
% trx.fieldnamen(firstframe) (double)
% ...
% trx.fieldnamen(endframe) (double)

% header
fwrite(fid,version,'double');
fwrite(fid,length(moviename),'int');
fwrite(fid,moviename,'char');
fwrite(fid,length(matname),'int');
fwrite(fid,matname,'char');
fwrite(fid,length(flies),'int');
fwrite(fid,flies,'int');
fwrite(fid,t0-1,'int');
fwrite(fid,t1-1,'int'); 

% sex = unique(trx(idx).sex);
% if length(sex) > 1 % CSC 20110317: require that sex is constant per fly in registered trx 
%     error('Error: Violated assumption in %s: fly %d DOES change sex within regeistered trx\n', moviename);
% end
% 
% fwrite(fid,sex{:},'char'); % CSC 20110317: assumes that sex is constant per fly in registered trx 
% % fwrite(fid,trx(idx).sex,'char');
allSex = trx.sex;
fwrite(fid,allSex{1},'char'); % CSC 20110318: write first value only

fwrite(fid,trx.fps,'double');

% number of fields
nframes = trx(idx).nframes;
%fns = setdiff(fieldnames(trx),fnsignore);
nfields = 1; % count timestamps
for i = 1:length(fns),
  fn = fns{i};
  l = length(trx(idx).(fn));
  if l < nframes - 2 || l > nframes,
    continue;
  end
  nfields = nfields + 1;
end
fwrite(fid,nfields,'int');

% timestamps
%if ~isfield(trx,'timestamps'),
%  [stamps] = fmf_read_timestamps( moviename, t0, t1 );
%else

%  ext0 = floor((nframes-l)/2);
%  off = t0 - trx(idx).firstframe;
%  stamps = padgrab(double(trx(idx).timestamps),nan,off+t0-ext0,off+t1);
    stamps = trx(idx).timestamps(t0:t1);
  %stamps = trx(idx).timestamps;
%end

%fn = 'timestamps';
fn = 'timestamp';
fwrite(fid,length(fn),'int');
fwrite(fid,fn,'char');
units.num = 's';
units.den = '1';
fwrite(fid,length(units.num),'int');
fwrite(fid,units.num,'char');
fwrite(fid,length(units.den),'int');
fwrite(fid,units.den,'char');
fwrite(fid,stamps,'double');

% data
%fns = setdiff(fieldnames(trx),fnsignore);
for i = 1:length(fns),
  fn = fns{i};
  l = length(trx(idx).(fn));
  if l < nframes - 2 || l > nframes,
    continue;
  end
  fwrite(fid,length(fn),'int');
  fwrite(fid,fn,'char');
  units = struct;
  if ~isfield(trx.units,fn),
%  if ~isfield(trx(idx), 'units') || ~isfield(trx(idx).units,fn), % CSC 20110317: check whether field 'units' exists as well
    units.num = '1';
    units.den = '1';
  else
    if isempty(trx.units.(fn).num),
      units.num = '1';
    else
      units.num = sprintf('%s ',trx.units.(fn).num{:});
    end
    if isempty(trx.units.(fn).den),
      units.den = '1';
    else
      units.den = sprintf('%s ',trx.units.(fn).den{:});
    end
  end
  fwrite(fid,length(units.num),'int');
  fwrite(fid,units.num,'char');
  fwrite(fid,length(units.den),'int');
  fwrite(fid,units.den,'char');
    
%   ext0 = floor((nframes-l)/2);
%   %ext1 = nframes - l - ext0;
%   off = t0 - trx(idx).firstframe;
%   data = padgrab(double(trx(idx).(fn)),nan,off+t0-ext0,off+t1);  
    values = trx(idx).(fn);
    data = values(t0:t1);
  fwrite(fid,data,'double');  
end
  
fclose(fid);
