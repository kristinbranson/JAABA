function [trx,fieldnames,data,units] = ReadStructSVMTrxFile(trxfile)

fid = fopen(trxfile,'rb');
version = fread(fid,1,'double'); %#ok<NASGU>
moviename_length = fread(fid,1,'int');
moviename = fread(fid,moviename_length,'char=>char')';
matname_length = fread(fid,1,'int');
matname = fread(fid,matname_length,'char=>char')'; %#ok<NASGU>
nflies = fread(fid,1,'int'); %#ok<NASGU>
fly_ids = fread(fid,1,'int');
firstframe = fread(fid,1,'int');
firstframe = firstframe + 1;
lastframe = fread(fid,1,'int');
lastframe = lastframe + 1;
sex = fread(fid,1,'char=>char');
fps = fread(fid,1,'double');
nfields = fread(fid,1,'int');

fieldnames = cell(1,nfields);
data = cell(1,nfields);
units = cell(1,nfields);
nframes = lastframe-firstframe+1;

for i = 1:nfields,
  fieldname_length = fread(fid,1,'int');
  fieldname = fread(fid,fieldname_length,'char=>char')';
  fieldnames{i} = fieldname;
  units_numerator_length = fread(fid,1,'int');
  units_numerator = fread(fid,units_numerator_length,'char=>char')';
  units_numerator = regexp(strtrim(units_numerator),'\s','split');
  idx = strcmp(units_numerator,'1');
  units_numerator(idx) = [];
  
  units_denominator_length = fread(fid,1,'int');
  units_denominator = fread(fid,units_denominator_length,'char=>char')';
  units_denominator = regexp(strtrim(units_denominator),'\s','split');
  idx = strcmp(units_denominator,'1');
  units_denominator(idx) = [];

  
  units{i} = struct('numerator',{units_numerator},'denominator',{units_denominator});
  datacurr = fread(fid,nframes,'double');
  data{i} = datacurr(:)';  
end

fclose(fid);

trxfns = {'x','y','a','b','theta','x_mm','y_mm','a_mm','b_mm','theta_mm','timestamps'};

trx = struct;
for i = 1:numel(trxfns),
  fn = trxfns{i};
  j = find(strcmp(fn,fieldnames),1);
  if isempty(j), continue; end
  trx.(fn) = data{j};
end

if isfield(trx,'theta') && ~isfield(trx,'theta_mm') && isfield(trx,'x_mm'),
  trx.theta_mm = trx.theta;
end

j = find(strcmp('timestamp',fieldnames),1);
if ~isfield(trx,'timestamps') && ~isempty(j),
  trx.timestamps = data{j};
end
if ~isfield(trx,'timestamps'),
  trx.dt = diff(trx.timestamps);
end
if ~isfield(trx,'dt'),
  trx.dt = repmat(1/fps,[1,nframes-1]);
end
if ~isfield(trx,'timestamps'),
  trx.timestamps = cumsum([(firstframe-1)/fps,trx.dt]);
end

trx.firstframe = firstframe;
trx.off = 1-trx.firstframe;
trx.endframe = lastframe;
trx.nframes = nframes;
trx.id = fly_ids;
trx.moviename = moviename;
trx.sex = sex;
trx.fps = fps;