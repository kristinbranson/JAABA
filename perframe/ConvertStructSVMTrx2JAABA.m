function [success,msg] = ConvertStructSVMTrx2JAABA(trxfiles,labelfiles,inmoviefile,expdir,varargin)

global CSSVM_INTRKFILE CSSVM_TD;

success = false;
msg = {};

wingunits = struct(...
  'nwingsdetected',parseunits('unit'),...
  'wing_areal',parseunits('px^2'),...
  'wing_arear',parseunits('px^2'),...
  'wing_trough_angle',parseunits('rad'));

[trxfilestr,moviefilestr,perframedirstr,intrkfile,...
  quartermaj,quartermin,...
  infofile,pxpermm,dosoftlink,...
  duplicatelabels,fracsubsamplenegs,...
  cachetrk] = ...
  myparse(varargin,...
  'trxfilestr','trx.mat','moviefilestr','movie.seq',...
  'perframedirstr','perframe',...
  'intrkfile','',...
  'a',7,'b',2,...
  'infofile','',...
  'pxpermm',11.84,...
  'dosoftlink',false,...
  'duplicatelabels',false,...
  'fracsubsamplenegs',1,...
  'cachetrk',false);

if numel(trxfiles) ~= numel(labelfiles),
  msg = {'Number of trx and label files do not match'};
  return;
end

trx = [];
alldata = {};
allfieldnames = {};
allunits = {};
alllabels = [];
allgtlabels = [];
allbehaviors = {};

% loop through trxfiles and labelfiles, one for each fly
for i = 1:numel(trxfiles),
  
  trxfile = trxfiles{i};
  labelfile = labelfiles{i};
  
  % read the trx file
  [trxcurr,fieldnamescurr,datacurr,unitscurr] = ReadStructSVMTrxFile(trxfile);
  trx = structappend(trx,trxcurr);
  id = trxcurr.id;
  if numel(trxfiles) == 1,
    id = 1;
  end
  
  % check that it is from the same movie
  if i == 1,
    moviefile = trxcurr.moviename;
  else
    if ~strcmp(moviefile,trxcurr.moviename),
      msg = {sprintf('Movies mismatch %d %s %s', i, trxfile, labelfile) };
      return;
    end
  end

  % add on the perframe data
  [ism,idx] = ismember(fieldnamescurr,allfieldnames);
  if any(ism),
    for j = find(ism),
      alldata{idx(j)}{id} = datacurr{j};
    end
  end
  if any(~ism),
    for j = find(~ism),
      alldata{end+1} = {}; %#ok<*AGROW>
      alldata{end}{id} = datacurr{j};
      allunits{end+1} = unitscurr{j};
      allfieldnames{end+1} = fieldnamescurr{j};
    end
  end
  
  % read the labels    
  if ~isempty(labelfile) && exist(labelfile, 'file')
      [labelscurr,behaviorscurr] = ReadStructSVMLabelFile(labelfile);

      % randomly subsample negatives
      if fracsubsamplenegs < 1,

        labels0 = labelscurr;

        for behaviori = 1:numel(labelscurr),

          maxT = max(labelscurr{behaviori}.t1s);
          label = zeros(1,maxT);
          for j = 1:numel(labelscurr{behaviori}.t0s),
            if strcmpi(labelscurr{behaviori}.names{j},'none'),
              label(labelscurr{behaviori}.t0s(j):labelscurr{behaviori}.t1s(j)-1) = 2;
            else
              label(labelscurr{behaviori}.t0s(j):labelscurr{behaviori}.t1s(j)-1) = 1;
            end
          end

          dokeep = label ~= 2 | rand(1,maxT)<=fracsubsamplenegs;
          label(~dokeep) = 0;

          newlabels = struct('t0s',[],'t1s',[],'names',{{}},'flies',[],'timestamp',[],'imp_t0s',[],'imp_t1s',[],'off',[]);
          labelnames_curr = {behaviorscurr{behaviori},'None'};
          for j = 1:2,
            [i0s,i1s] = get_interval_ends(label==j);
            if ~isempty(i0s),
              n = numel(i0s);
              newlabels.t0s(end+1:end+n) = i0s;
              newlabels.t1s(end+1:end+n) = i1s;
              newlabels.imp_t0s(end+1:end+n) = i0s;
              newlabels.imp_t1s(end+1:end+n) = i1s;
              newlabels.names(end+1:end+n) = repmat(labelnames_curr(j),[1,n]);
            end
          end
          newlabels.timestamp = repmat(labelscurr{behaviori}.timestamp(end),size(newlabels.t0s));
          newlabels.flies = labelscurr{behaviori}.flies;
          newlabels.off = labelscurr{behaviori}.off;
          labelscurr{behaviori} = newlabels;

        end

      end


      % add them to alllabels
      [ism,idx] = ismember(behaviorscurr,allbehaviors);
      if any(ism),
        for j = find(ism),
          alllabels{idx(j)}.t0s{end+1} = labelscurr{j}.t0s;
          alllabels{idx(j)}.t1s{end+1} = labelscurr{j}.t1s;
          alllabels{idx(j)}.names{end+1} = labelscurr{j}.names;
          alllabels{idx(j)}.flies(:,end+1) = labelscurr{j}.flies;
          alllabels{idx(j)}.off(end+1) = labelscurr{j}.off;
          alllabels{idx(j)}.timestamp{end+1} = labelscurr{j}.timestamp;
          alllabels{idx(j)}.imp_t0s{end+1} = labelscurr{j}.imp_t0s;
          alllabels{idx(j)}.imp_t1s{end+1} = labelscurr{j}.imp_t1s;
        end
        for j = find(ism),
          allgtlabels{idx(j)}.t0s{end+1} = labels0{j}.t0s;
          allgtlabels{idx(j)}.t1s{end+1} = labels0{j}.t1s;
          allgtlabels{idx(j)}.names{end+1} = labels0{j}.names;
          allgtlabels{idx(j)}.flies(:,end+1) = labels0{j}.flies;
          allgtlabels{idx(j)}.off(end+1) = labels0{j}.off;
          allgtlabels{idx(j)}.timestamp{end+1} = labels0{j}.timestamp;
          allgtlabels{idx(j)}.imp_t0s{end+1} = labels0{j}.imp_t0s;
          allgtlabels{idx(j)}.imp_t1s{end+1} = labels0{j}.imp_t1s;
        end

      end
      if any(~ism),
        for j = find(~ism),
          alllabels{end+1} = labelscurr{j};
          alllabels{end}.t0s = {alllabels{end}.t0s};
          alllabels{end}.t1s = {alllabels{end}.t1s};
          alllabels{end}.names = {alllabels{end}.names};
          alllabels{end}.timestamp = {alllabels{end}.timestamp};
          alllabels{end}.imp_t0s = {alllabels{end}.imp_t0s};
          alllabels{end}.imp_t1s = {alllabels{end}.imp_t1s};
          allbehaviors{end+1} = behaviorscurr{j};
        end
        for j = find(~ism),
          allgtlabels{end+1} = labels0{j};
          allgtlabels{end}.t0s = {allgtlabels{end}.t0s};
          allgtlabels{end}.t1s = {allgtlabels{end}.t1s};
          allgtlabels{end}.names = {allgtlabels{end}.names};
          allgtlabels{end}.timestamp = {allgtlabels{end}.timestamp};
          allgtlabels{end}.imp_t0s = {allgtlabels{end}.imp_t0s};
          allgtlabels{end}.imp_t1s = {allgtlabels{end}.imp_t1s};
          %allbehaviors{end+1} = behaviorscurr{j};
        end
      end
  end
end

msg{end+1} = sprintf('Read %d trx from trx files',numel(trx));
msg{end+1} = sprintf('Read %d behavior labels',numel(alllabels));

% sort trx by id
[~,order] = sort([trx.id]);
trx = trx(order);

% create timestamps
timestamps = nan(1,max([trx.endframe]));
for i = 1:numel(trx),
  timestamps(trx(i).firstframe:trx(i).endframe) = trx(i).timestamps;
end

if ~isempty(infofile),
  id = load(infofile);
  pxpermm = id.hinfo.PPM;
  b = bwboundaries(id.hinfo.mask);
  if numel(b) > 1,
    msg = {'Do not know how to handle multiple ROIs'};
    return;
  end
  idx = convhull(b{1}(:,2),b{1}(:,1),'simplify',true);
  x = b{1}(idx,2);
  y = b{1}(idx,1);
  [nr,nc,~] = size(id.hinfo.mask);
  [~,tl] = min( x.^2 + y.^2 );
  [~,tr] = min( (nc-x).^2 + y.^2 );
  [~,bl] = min( x.^2 + (nr-y).^2 );
  [~,br] = min( (nc-x).^2 + (nr-y).^2 );
  for fly = 1:numel(trx),
    trx(fly).arena = struct('tl',[x(tl),y(tl)]/pxpermm,'tr',[x(tr),y(tr)]/pxpermm,...
      'bl',[x(bl),y(bl)]/pxpermm,'br',[x(br),y(br)]/pxpermm);
  end
  msg{end+1} = sprintf('Read pxperm = %f from %s',pxpermm,infofile);
  msg{end+1} = sprintf('Read rectangular arena position from %s',infofile);
  
end
for fly = 1:numel(trx),
  trx(fly).pxpermm = pxpermm;
end

perframedata = struct;
if ~isempty(intrkfile),
  
  if cachetrk,
    if ~isempty(CSSVM_INTRKFILE) && ~isempty(CSSVM_TD) && ...
        strcmp(intrkfile,CSSVM_INTRKFILE),
      td = CSSVM_TD;
    else
      td = load(intrkfile);
      CSSVM_INTRKFILE = intrkfile;
      CSSVM_TD = td;
    end
  else
    td = load(intrkfile);
  end
  nflies = numel(td.trk.sequences);
  nflies_trx = numel(trx);
  if nflies_trx >= nflies,
    msg = {sprintf('N. flies read from %s = %d < n. flies read from trxfiles = %d',intrkfile,nflies,nflies_trx)};
    return;
  end
  % may need to expand trx
  if nflies > nflies_trx,
    if nflies_trx ~= 1,
      msg{end+1} = sprintf('N. flies read from %s = %d > n. flies read from trxfiles = %d',intrkfile,nflies,nflies_trx);
    end
    trx = [trx,repmat(trx(end),[1,nflies-nflies_trx])];
    for i = 1:numel(alldata),
      alldata{i} = [alldata{i},repmat(alldata{i}{end},[1,nflies-nflies_trx])];
    end
  end
  
  if duplicatelabels,
    for i = 1:numel(alllabels),
      for fly = nflies_trx+1:nflies,
        alllabels{i}.t0s{fly} = alllabels{i}.t0s{nflies_trx};
        alllabels{i}.t1s{fly} = alllabels{i}.t1s{nflies_trx};
        alllabels{i}.names{fly} = alllabels{i}.names{nflies_trx};
        alllabels{i}.timestamp{fly} = alllabels{i}.timestamp{nflies_trx};
        alllabels{i}.imp_t0s{fly} = alllabels{i}.imp_t0s{nflies_trx};
        alllabels{i}.imp_t1s{fly} = alllabels{i}.imp_t1s{nflies_trx};
        alllabels{i}.off(fly) = alllabels{i}.off(nflies_trx);
        alllabels{i}.flies(fly,:) = fly;
      end
    end
    for i = 1:numel(allgtlabels),
      for fly = nflies_trx+1:nflies,
        allgtlabels{i}.t0s{fly} = allgtlabels{i}.t0s{nflies_trx};
        allgtlabels{i}.t1s{fly} = allgtlabels{i}.t1s{nflies_trx};
        allgtlabels{i}.names{fly} = allgtlabels{i}.names{nflies_trx};
        allgtlabels{i}.timestamp{fly} = allgtlabels{i}.timestamp{nflies_trx};
        allgtlabels{i}.imp_t0s{fly} = allgtlabels{i}.imp_t0s{nflies_trx};
        allgtlabels{i}.imp_t1s{fly} = allgtlabels{i}.imp_t1s{nflies_trx};
        allgtlabels{i}.off(fly) = allgtlabels{i}.off(nflies_trx);
        allgtlabels{i}.flies(fly,:) = fly;
      end
    end
  end
  
  perframedata.nwingsdetected = cell(1,nflies);
  perframedata.wing_areal = cell(1,nflies);
  perframedata.wing_arear = cell(1,nflies);
  perframedata.wing_trough_angle = cell(1,nflies);
  
  for fly = 1:nflies,
    firstframe = max(trx(fly).firstframe,td.trk.sequences{fly}.time_start);
    endframe = min(trx(fly).endframe,td.trk.sequences{fly}.time_end);
    nframes = endframe-firstframe+1;
    
    % shorten trx
    off = 1-trx(fly).firstframe;
    trx(fly).timestamps = trx(fly).timestamps(firstframe+off:endframe+off);
    trx(fly).dt = trx(fly).dt(firstframe+off:endframe-1+off);
    
    % from trk to trx
    off = -td.trk.sequences{fly}.time_start+1;
    trx(fly).x = td.trk.sequences{fly}.pos.x(firstframe+off:endframe+off)';
    trx(fly).y = td.trk.sequences{fly}.pos.y(firstframe+off:endframe+off)';
    trx(fly).theta = -td.trk.sequences{fly}.pos.ori(firstframe+off:endframe+off)';
    trx(fly).a = repmat(quartermaj,[1,nframes]);
    trx(fly).b = repmat(quartermin,[1,nframes]);
    
    % wings
    for t = firstframe:endframe,
      wings = td.trk.sequences{fly}.pos.wings{t+off};
      if isempty(wings),
        wing_angles = 0;
        wing_lengths = trx(fly).a(t-firstframe+1)*2.5;
      else
        
        % angle between the wings and the tail of the fly
        wing_angles = nan(1,numel(wings));
        wing_lengths = nan(1,numel(wings));
        for wingi = 1:numel(wings),
          dx = wings{wingi}(1)-trx(fly).x(t-firstframe+1);
          dy = wings{wingi}(2)-trx(fly).y(t-firstframe+1);
          wing_angles(wingi) = modrange(atan2(dy,dx)-trx(fly).theta(t-firstframe+1)-pi,-pi,pi);
          wing_lengths(wingi) = sqrt(dx^2+dy^2);
        end
      end
      [trx(fly).wing_anglel(t-firstframe+1),l] = max(wing_angles);
      [trx(fly).wing_angler(t-firstframe+1),r] = min(wing_angles);
      if ~isempty(wings),
        trx(fly).xwingl(t-firstframe+1) = wings{l}(1);
        trx(fly).ywingl(t-firstframe+1) = wings{l}(2);
        trx(fly).xwingr(t-firstframe+1) = wings{r}(1);
        trx(fly).ywingr(t-firstframe+1) = wings{r}(2);
      else
        trx(fly).xwingl(t-firstframe+1) = trx(fly).x(t-firstframe+1) + cos(-trx(fly).theta(t-firstframe+1)*wing_lengths);
        trx(fly).xwingr(t-firstframe+1) = trx(fly).xwingl(t-firstframe+1);
        trx(fly).ywingl(t-firstframe+1) = trx(fly).y(t-firstframe+1) + sin(-trx(fly).theta(t-firstframe+1)*wing_lengths);
        trx(fly).ywingr(t-firstframe+1) = trx(fly).ywingl(t-firstframe+1);
      end
      perframedata.nwingsdetected{fly}(t-firstframe+1) = numel(wings);
  
      % this is really the wing length
      perframedata.wing_areal{fly}(t-firstframe+1) = wing_lengths(l);
      perframedata.wing_arear{fly}(t-firstframe+1) = wing_lengths(r);
  
    end
      
    % set trough angle just to be between the two wings
    perframedata.wing_trough_angle{fly} = (trx(fly).wing_angler - trx(fly).wing_anglel)/2;
    
    % overwrite start and end
    trx(fly).firstframe = firstframe;
    trx(fly).endframe = endframe;
    trx(fly).nframes = nframes;
    trx(fly).off = 1-firstframe;
  
  end
  
  msg{end+1} = sprintf('Read trajectories for %d flies from %s',numel(trx),intrkfile);
end

outmoviefile = fullfile(expdir,moviefilestr);

for fly = 1:numel(trx),
  if isfield(trx,'x'),
    trx(fly).x_mm = trx(fly).x/pxpermm;
  end
  if isfield(trx,'y'),
    trx(fly).y_mm = trx(fly).y/pxpermm;
  end
  if isfield(trx,'a'),
    trx(fly).a_mm = trx(fly).a/pxpermm;
  end
  if isfield(trx,'b'),
    trx(fly).b_mm = trx(fly).b/pxpermm;
  end
  if isfield(trx,'theta'),
    trx(fly).theta_mm = trx(fly).theta;
  end
  trx(fly).moviefile = outmoviefile;
end

if ~exist(expdir,'dir'),
  mkdir(expdir);
end

% save trx
trxfile = fullfile(expdir,trxfilestr);
save(trxfile,'trx','timestamps');

% create perframe dir
perframedir = fullfile(expdir,perframedirstr);
if ~exist(perframedir,'dir'),
  mkdir(perframedir);
end

% save perframe data
for i = 1:numel(alldata),
  data = alldata{i}; %#ok<NASGU>
  fn = allfieldnames{i};
  units = allunits{i}; %#ok<NASGU>
  filename = fullfile(perframedir,[fn,'.mat']);
  save(filename,'data','units');
end


fns = fieldnames(perframedata);
for i = 1:numel(fns),
  fn = fns{i};
  s = struct('data',{perframedata.(fn)},'units',wingunits.(fn)); %#ok<NASGU>
  filename = fullfile(perframedir,[fn,'.mat']);
  save(filename,'-struct','s');
end


% save labels
for i = 1:numel(alllabels),
  filestr = sprintf('labeled_%s.mat',allbehaviors{i});
  labels = alllabels{i}; %#ok<NASGU>
  save(fullfile(expdir,filestr),'-struct','labels');
end
for i = 1:numel(allgtlabels),
  filestr = sprintf('labeled_gt_%s.mat',allbehaviors{i});
  labels = allgtlabels{i}; %#ok<NASGU>
  save(fullfile(expdir,filestr),'-struct','labels');
end

if ~isempty(inmoviefile),

  if ~exist(inmoviefile,'file'),
    msg = {sprintf('Input movie file %s does not exist',inmoviefile)};
    return;
  end
  
  % copy/soft-link movie
  if dosoftlink,
    if exist(outmoviefile,'file'),
      delete(outmoviefile);
    end
    if isunix,
      cmd = sprintf('ln -s %s %s',inmoviefile,outmoviefile);
      unix(cmd);
      % test to make sure it worked
      [status,result] = unix(sprintf('readlink %s',outmoviefile));
      result = strtrim(result);
      if status ~= 0 || ~strcmp(result,inmoviefile),
        res = questdlg(sprintf('Failed to make soft link. Copy %s to %s instead?',inmoviefile,outmoviefile));
        if ~strcmpi(res,'Yes'),
          msg = sprintf('Failed to make soft link from %s to %s.',inmoviefile,outmoviefile);
          return;
        end
        dosoftlink = false;
      end
    elseif ispc,
      if exist([outmoviefile,'.lnk'],'file'),
        delete([outmoviefile,'.lnk']);
      end
      cmd = sprintf('mkshortcut.vbs /target:"%s" /shortcut:"%s"',inmoviefile,outmoviefile);
      fprintf('Making a Windows shortcut file at "%s" with target "%s"\n',outmoviefile,inmoviefile);
      system(cmd);
      % test to make sure that worked
      [equalmoviefile,didfind] = GetPCShortcutFileActualPath(outmoviefile);
      if ~didfind || ~strcmp(equalmoviefile,inmoviefile),
        res = questdlg(sprintf('Failed to make shortcut. Copy %s to %s instead?',inmoviefile,outmoviefile));
        if ~strcmpi(res,'Yes'),
          msg = sprintf('Failed to make shortcut from %s to %s.',inmoviefile,outmoviefile);
          return;
        end
        dosoftlink = false;
      end
    else
      res = questdlg(sprintf('Unknown OS, cannot soft-link movie file %s. Copy instead?',inmoviefile));
      if ~strcmpi(res,'Yes'),
        msg = sprintf('Failed to make softlink from %s to %s.',inmoviefile,outmoviefile);
        return;
      end
      dosoftlink = false;
    end
    if dosoftlink,
      msg{end+1} = sprintf('Made a link to movie file %s at %s',inmoviefile,outmoviefile);
    end
  end
  
  if ~dosoftlink,
    if ispc,
      if exist([outmoviefile,'.lnk'],'file'),
        delete([outmoviefile,'.lnk']);
      end
    end
    if exist(outmoviefile,'file'),
      delete(outmoviefile);
    end
    [success1,msg1] = copyfile(inmoviefile,outmoviefile);
    if ~success1,
      msg = msg1;
      success = false;
      return;
    end
    msg{end+1} = sprintf('Copied movie file %s to %s',inmoviefile,outmoviefile);
  end
  
end

success = true;