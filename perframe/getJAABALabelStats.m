function info = getJAABALabelStats()
%%
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

expinfo = regexp(expnames,'^(?<type>.*)(_Rig(?<rig>\w))?_(?<timestamp>\d{8}T\d{6})$','names','once');  
  
expinfo = [expinfo{:}];
[info.types,~,typeidx] = unique({expinfo.type});

info.nlabelframes = zeros(numel(info.types),nbehaviors);
info.nlabelbouts = zeros(numel(info.types),nbehaviors);
info.nlabelflies = zeros(numel(info.types),nbehaviors);
info.nlabelexps = zeros(numel(info.types),nbehaviors);

for typei = 1:numel(info.types),
  
  idxcurr = find(typeidx==typei);
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
      info.nlabelflies(typei,counts>0) = info.nlabelflies(typei,counts>0) + 1;
      info.nlabelbouts(typei,:) = info.nlabelbouts(typei,:) + counts;
      islabeledexp(counts>0) = true;
      for behi = 1:nbehaviors,
        info.nlabelframes(typei,behi) = info.nlabelframes(typei,behi) + sum((labels(expi).t1s{flyi}(labelidx==behi)-labels(expi).t0s{flyi}(labelidx==behi)));
      end
    end
    info.nlabelexps(typei,islabeledexp) = info.nlabelexps(typei,islabeledexp) + 1;
  end
  
end

typelens = cellfun(@numel,info.types);
maxtypelen = max(typelens);
namelens = cellfun(@numel,info.allnames);
maxnamelen = max(namelens);
for typei = 1:numel(info.types),
  for behi = 1:nbehaviors,
    fprintf('%s%s %s%s:\t',info.types{typei},repmat(' ',[1,maxtypelen-typelens(typei)]),...
      info.allnames{behi},repmat(' ',[1,maxnamelen-namelens(behi)]));
    fprintf('%2d exps, %3d flies, %6d bouts, %6d frames\n',info.nlabelexps(typei,behi),...
      info.nlabelflies(typei,behi),info.nlabelbouts(typei,behi),info.nlabelframes(typei,behi));
  end
end

