function rawdata = ReadRawDataFile(datafile)

metadatafns = {'mouse','exp','id','upper'};
metadataisnum = [true,false,false,false];

datafnsspecial = {'auto_Grab_success','auto_Grab_successtype'};

fid = fopen(datafile,'r');

% read headers
l = fgetl(fid);
headers = regexp(l,',','split');
idxdata = find(~cellfun(@isempty,regexp(headers,'^auto','once')) & ...
  ~ismember(headers,datafnsspecial));% & ...
  %cellfun(@isempty,regexp(headers,'successfultrial')));
metadatai = nan(1,numel(metadatafns));
for i = 1:numel(metadatafns),
  metadatai(i) = find(strcmpi(headers,metadatafns{i}));
end

successi = find(strcmpi(headers,'auto_Grab_success'),1);
successtypei = find(strcmpi(headers,'auto_Grab_successtype'),1);

datei = find(strcmpi(headers,'date'),1);
parentdiri = find(strcmpi(headers,'upper'),1);
expi = find(strcmpi(headers,'exp'),1);

datafns = headers(idxdata);

rawdata = [];


while true,
  l = fgetl(fid);
  if ~ischar(l),
    break;
  end
  if isempty(l),
    continue;
  end
  s = regexp(l,',','split');
  
  rawdatacurr = struct;
  
  m = regexp(s{datei},'^(\d+)([^\d]*)$','once','tokens');
  rawdatacurr.date = m{1};
  rawdatacurr.iscno = strcmpi(m{2},'CNO');
  m = regexp(s{parentdiri},'[\\/]([^\\/]+)$','tokens','once'); 
  rawdatacurr.session = m{1};
  m = regexp(s{expi},'_[vV](\d+)$','tokens','once');
  rawdatacurr.trial = str2double(m{1});

  
  for i = 1:numel(metadatafns),
    fn = metadatafns{i};
    rawdatacurr.(fn) = s{metadatai(i)};
    if metadataisnum(i),
      rawdatacurr.(fn) = str2double(rawdatacurr.(fn));
    end
  end
  
  rawdatacurr.success = find(strcmpi(s{successi},{'false','true'}))-1;
  rawdatacurr.successtype = s{successtypei};

  for i = 1:numel(datafns),
    fn = datafns{i};
    rawdatacurr.(fn) = str2double(s{idxdata(i)});
  end
  
  if isfield(rawdatacurr,'exp') && ~isfield(rawdatacurr,'expfull'),    
    rawdatacurr.expfull = rawdatacurr.exp;
  end
    
  rawdata = structappend(rawdata,rawdatacurr);
  
  
end

fclose(fid);