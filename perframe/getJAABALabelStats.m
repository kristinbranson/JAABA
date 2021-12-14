function info = getJAABALabelStats(expinfo,info)
%%

isinfo = nargin >= 2;

hfig = findall(0,'type','figure','Name','JAABA');
assert(numel(hfig)==1);


gdata = guidata(hfig);
data = gdata.data;
data.GetLabels();
labels = data.labels;
expnames = data.expnames;
behaviors = data.classifiernames;
nonebehavior = data.nobehaviornames;
info.allnames = [behaviors,nonebehavior];
nbehaviors = numel(info.allnames);

if nargin >= 1,
  expnames0 = cellfun(@fileBaseName,{expinfo.file_system_path},'Uni',0);
  [ism,idx] = ismember(expnames,expnames0);
  assert(all(ism));
  expinfo = expinfo(idx);  
else
  expinfo = regexp(expnames,'^(?<label>.*)_Rig(?<rig>\w)_(?<timestamp>\d{8}T\d{6})$','names','once');
  expinfo = [expinfo{:}];
end

if ~isinfo,
  [info.labels,~,labelidx0] = unique({expinfo.label});
  
  info.nlabelframes = zeros(numel(info.labels),nbehaviors);
  info.nlabelbouts = zeros(numel(info.labels),nbehaviors);
  info.nlabelflies = zeros(numel(info.labels),nbehaviors);
  info.nlabelexps = zeros(numel(info.labels),nbehaviors);
  info.ntotalexps = zeros(numel(info.labels),1);
  info.nlabelframesperexp = zeros(numel(expinfo),nbehaviors);
  info.nlabelboutsperexp = zeros(numel(expinfo),nbehaviors);
  info.nlabelfliesperexp = zeros(numel(expinfo),nbehaviors);
  
  for labeli = 1:numel(info.labels),
    
    idxcurr = find(labelidx0==labeli);
    info.ntotalexps(labeli) = numel(idxcurr);
    for jj = 1:numel(idxcurr),
      islabeledexp = false(1,nbehaviors);
      expi = idxcurr(jj);
      for flyi = 1:numel(labels(expi).t0s),
        if isempty(labels(expi).t0s(flyi)),
          continue;
        end
        if isempty(labels(expi).names{flyi}), continue; end
        [ism,labelidx] = ismember(labels(expi).names{flyi},info.allnames);
        assert(all(ism));
        counts = hist(labelidx,1:nbehaviors);
        info.nlabelflies(labeli,counts>0) = info.nlabelflies(labeli,counts>0) + 1;
        info.nlabelbouts(labeli,:) = info.nlabelbouts(labeli,:) + counts;
        info.nlabelfliesperexp(expi,counts>0) = info.nlabelfliesperexp(expi,counts>0) + 1;
        info.nlabelboutsperexp(expi,:) = info.nlabelboutsperexp(expi,:) + counts;
        islabeledexp(counts>0) = true;
        for behi = 1:nbehaviors,
          ncurr = sum((labels(expi).t1s{flyi}(labelidx==behi)-labels(expi).t0s{flyi}(labelidx==behi)));
          info.nlabelframes(labeli,behi) = info.nlabelframes(labeli,behi) + ncurr;
          info.nlabelframesperexp(expi,behi) = info.nlabelframesperexp(expi,behi) + ncurr;
        end
      end
      info.nlabelexps(labeli,islabeledexp) = info.nlabelexps(labeli,islabeledexp) + 1;
    end
    
  end
end

fprintf('Labels per experiment type:\n');
labellens = cellfun(@numel,info.labels);
maxlabellen = max(labellens);
namelens = cellfun(@numel,info.allnames);
maxnamelen = max(namelens);
for labeli = 1:numel(info.labels),
  for behi = 1:nbehaviors,
    fprintf('%s%s %s%s:\t',info.labels{labeli},repmat(' ',[1,maxlabellen-labellens(labeli)]),...
      info.allnames{behi},repmat(' ',[1,maxnamelen-namelens(behi)]));
    fprintf('%2d/%2d exps, %3d flies, %6d bouts, %6d frames\n',info.nlabelexps(labeli,behi),info.ntotalexps(labeli),...
      info.nlabelflies(labeli,behi),info.nlabelbouts(labeli,behi),info.nlabelframes(labeli,behi));
  end
end

fprintf('\nLabels per experiment:\n');
expnamelens = cellfun(@numel,expnames);
maxexpnamelen = max(expnamelens);
for expi = 1:numel(expinfo),
  labeli = find(strcmp(info.labels,expinfo(expi).label),1);
  for behi = 1:nbehaviors,
    fprintf('%s%s %s%s %s%s:\t',expinfo(expi).label,repmat(' ',[1,maxlabellen-labellens(labeli)]),expnames{expi},repmat(' ',[1,maxexpnamelen-expnamelens(expi)]),...
      info.allnames{behi},repmat(' ',[1,maxnamelen-namelens(behi)]));
    fprintf('%3d flies, %6d bouts, %6d frames\n',info.nlabelfliesperexp(expi,behi),info.nlabelboutsperexp(expi,behi),info.nlabelframesperexp(expi,behi));
  end
end