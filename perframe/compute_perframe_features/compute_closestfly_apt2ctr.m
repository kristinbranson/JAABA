function [data,units,mind] = compute_closestfly_apt2ctr(trx,aptdata,aptLD,n,dosave_d)

if nargin < 5,
  dosave_d = true;
end

flies = trx.exp2flies{n};
nflies = numel(flies);
closestfly = cell(1,nflies);
mind = cell(1,nflies);

parfor i1 = 1:nflies,
% for i1 = 1:nflies,
  fly1 = flies(i1);
  fprintf('target 1 = %d\n',fly1);

  flies2 = flies(trx.roi(fly1)==trx.roi(flies));

  mind{i1} = inf(1,trx(fly1).nframes);
  closesti = ones(1,trx(fly1).nframes);
  for i2 = 1:numel(flies2),
    fly2 = flies2(i2);
    if fly1 == fly2,
      continue;
    end
    [dcurr] = dapt2ctr_pair(trx,aptdata,aptLD,fly1,fly2);
    idx = dcurr < mind{i1};
    mind{i1}(idx) = dcurr(idx);
    closesti(idx) = i2;
  end
  closestfly{i1} = flies2(closesti);
  closestfly{i1}(isnan(mind{i1})|isinf(mind{i1})) = nan;
end
   
if dosave_d,
    data = mind; %#ok<NASGU>
    units = parseunits('units'); %#ok<NASGU>
    %   filename = trx.GetPerFrameFile('dnose2ell',n);
    perframedir = trx.perframedir;
    featurename = ['dapt',num2str(aptLD),'toctr.mat'];
    filename =  fullfile(trx.expdirs{n},perframedir,featurename);
    try
        save(filename,'data','units');
    catch ME,
        warning('Could not save file %s:\n%s',filename,getReport(ME));
    end
     % save closest fly 
    data = closestfly;
    units = parseunits('unit');
    featurename = ['closestfly_dapt',num2str(aptLD),'toctr.mat'];
    filename =  fullfile(trx.expdirs{n},perframedir,featurename);
    if exist(filename,'file'),
        delete(filename);
    else
        if isunix,
            [res,link] = unix(sprintf('readlink %s',filename));
            if ~res && ~isempty(link),
                warning('Deleting broken soft link from %s to %s.\n',filename,link);
                unix(sprintf('rm %s',filename));
            end
        end
    end
    save(filename,'data','units');

end
data = closestfly;
units = parseunits('unit');