function [fixeddirs,problemdirs,okdirs] = fixTrxFiles20140218(rootdirs,varargin)

TRXFILE = 'trx.mat';
FEATURESFILE = 'features.mat';
fns_fix = {'x','y','theta','a','b','x_mm','y_mm','a_mm','b_mm','theta_mm'};
  
fixeddirs = {};
problemdirs = {};
okdirs = {};

if nargin < 1 || isempty(rootdirs),
  rootdirs = uipickfiles('Prompt','Select root directories containing experiments to fix',...
    'DirsOnly',true);
  if ~ischar(rootdirs) || isempty(rootdirs),
    return;
  end
elseif ischar(rootdirs),
  rootdirs = {rootdirs};
end

[debug] = myparse(varargin,'debug',false);

for i = 1:numel(rootdirs),
  
  rootdir = rootdirs{i};
  trxfile = fullfile(rootdir,TRXFILE);
  if exist(trxfile,'file'),
    featuresfile = fullfile(rootdir,FEATURESFILE);    
    if exist(featuresfile,'file'),
      
      td = load(trxfile);
      nframes_trx = td.trx(1).nframes;
      fd = load(featuresfile);
      nframes_fea = size(fd.curFeatures,1);
      if nframes_trx > nframes_fea,
        fprintf('%s: n. frames in trx (%d) > n. frames in features (%d), fixing trx\n',rootdir,nframes_trx,nframes_fea);
        nframes_remove = nframes_trx-nframes_fea;
        for j = 1:numel(td.trx),
          for k = 1:numel(fns_fix),
            fn = fns_fix{k};
            td.trx(j).(fn) = td.trx(j).(fn)(nframes_remove+1:end);
          end
          td.trx(j).timestamps = td.trx(j).timestamps(1:end-nframes_remove);
          td.trx(j).dt = diff(td.trx(j).timestamps);
          td.trx(j).nframes = nframes_fea;
          td.trx(j).endframe = nframes_fea;
        end
        
        if debug == 0,          
          delete(trxfile);
          save(trxfile,'-struct','td');
        end
        
        fixeddirs{end+1} = rootdir; %#ok<AGROW>
      elseif nframes_trx < nframes_fea,

        warning('%s: n. frames in trx (%d) < n. frames in features (%d), not sure how to fix!\n',rootdir,nframes_trx,nframes_fea);
        problemdirs{end+1} = rootdir; %#ok<AGROW>
        
      else
        
        okdirs{end+1} = rootdir; %#ok<AGROW>
        
      end
        
    end
  else
    
    tmp = dir(rootdir);
    tmp = tmp([tmp.isdir]);
    tmp(ismember({tmp.name},{'.','..'})) = [];
    if ~isempty(tmp),
      newrootdirs = cellfun(@(x) fullfile(rootdir,x),{tmp.name},'UniformOutput',false);    
      [newfixeddirs,newproblemdirs,newokdirs] = fixTrxFiles20140218(newrootdirs,varargin{:});
      fixeddirs = [fixeddirs,newfixeddirs]; %#ok<AGROW>
      problemdirs = [problemdirs,newproblemdirs]; %#ok<AGROW>
      okdirs = [okdirs,newokdirs]; %#ok<AGROW>
    end
    
  end
  
end