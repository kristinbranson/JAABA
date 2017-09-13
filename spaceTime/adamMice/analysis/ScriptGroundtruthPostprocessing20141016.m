% ground truth post-processed labels

addpath ../miceCode;
setuppaths;
labeledjabfile = '/tier2/hantman/Adam/Jabfiles/Liftm119g_Handopenm119g_Grabm119g_Supm119g_Atmouthm119g_Chewm119g.jab';
classifierjabfile = '/tier2/hantman/Adam/Jabfiles/Liftm119p_Handopenm119p_Grabm119p_Supm119p_Atmouthm119p_Chewm119p.jab';
rootdatadir = '/tier2/hantman';

behaviornames = {'Lift','Handopen','Grab','Sup','Atmouth','Chew'};
% behaviorext = 'm119g';
% behaviorpre = 'No_';

fps = 500;


%% load in data

nbehaviors = numel(behaviornames);

load(labeledjabfile,'-mat');
jd = x;
clear x;

for i = 1:numel(jd.expDirNames),
  expdir = fullfile(rootdatadir,strrep(jd.expDirNames{i}(4:end),'\','/'));
  %assert(exist(expdir,'dir') > 0);
  jd.expDirNames{i} = expdir;
end

%expdirs = {'M119_20140326_v01'};
expdirs = jd.expDirNames;
nexps = numel(expdirs);

% get manual labels
expnames = cell(size(expdirs));
for i = 1:numel(expdirs),
  [~,expnames{i}] = fileparts(expdirs{i});
end
jdexpnames = cell(size(jd.expDirNames));
for i = 1:numel(jd.expDirNames),
  [~,jdexpnames{i}] = fileparts(jd.expDirNames{i});
end
[~,expidx] = ismember(expnames,jdexpnames);

manualfirstframes = nan(nbehaviors,nexps);
for jj = 1:nexps,
  j = expidx(jj);
  for i = 1:nbehaviors,
    ks = find(~cellfun(@isempty,regexp(jd.labels(j).names{1},behaviornames{i},'once')));
    %assert(numel(ks)<2);    
    if ~isempty(ks),
      manualfirstframes(i,jj) = min(jd.labels(j).t0s{1}(ks));
    end
  end
end

% get automatic, post-processed first frames
data = ExpPPHeadless(expdirs,{classifierjabfile});
%ethogram_plot(expdirs(1:2),{'Liftm119p_Handopenm119p_Grabm119p_Supm119p_Atmouthm119p_Chewm119p.jab'},[],'automarks',true)

autofirstframes = nan(nbehaviors,nexps);
for i = 1:nbehaviors,
  fn = sprintf('auto_%s_0',behaviornames{i});
  autofirstframes(i,:) = [data.(fn)];
end

% count errors
isfalsepos = isnan(manualfirstframes)&~isnan(autofirstframes);
isfalseneg = ~isnan(manualfirstframes)&isnan(autofirstframes);
npos = sum(~isnan(manualfirstframes),2);
nneg = sum(isnan(manualfirstframes),2);
nfalsepos = sum(isfalsepos,2);
nfalseneg = sum(isfalseneg,2);
falseposrate = nfalsepos ./ nneg;
falsenegrate = nfalseneg ./ npos;
for i = 1:nbehaviors,
  fprintf('%s: %d/%d missed detections, %d/%d spurious detections\n',behaviornames{i},...
    nfalseneg(i),npos(i),nfalsepos(i),nneg(i));
end
err_nframes = abs(autofirstframes-manualfirstframes);
err_ms = (err_nframes/fps)*1000;

for i = 1:nbehaviors,
  idxcurr = ~isnan(err_ms(i,:));
  fprintf('%s: mean error = %f ms, stddev of error = %f ms\n',behaviornames{i},mean(err_ms(i,idxcurr)),std(err_ms(i,idxcurr),1));
end

tmp = 2*unique(round(logspace(0,log10(100),25)));
%binedges = [-tmp(end:-1:1),0,tmp];
binedges = [0,tmp];
bincenters = (binedges(1:end-1)+binedges(2:end))/2;
nbins = numel(bincenters);
counts = nan(nbehaviors,nbins);
frac = nan(nbehaviors,nbins);
for i = 1:nbehaviors,
  
  idxcurr = ~isnan(err_ms(i,:));
  counts(i,:) = hist(err_ms(i,idxcurr),bincenters);
  frac(i,:) = counts(i,:) / nnz(idxcurr);
  
end

hfig = 1;
figure(hfig);
clf;
hax = axes('Position',[.05,.05,.9,.9]);
bar(frac');
xticklabels = cell(1,nbins);
for i = 1:nbins,
  if binedges(i) >= binedges(i+1)-2,
    xticklabels{i} = sprintf('%d',binedges(i));
  else    
    xticklabels{i} = sprintf('[%d,%d)',binedges(i),binedges(i+1));
  end
end
xticklabels{end} = sprintf('>= %d',binedges(end));
set(gca,'XTick',1:nbins,'XTickLabel',xticklabels,'XLim',[0,nbins+1]);
box off;
legend(behaviornames);
xlabel('Error between automatic and manual detections of first frame of behavior (ms)');
ylabel('Fraction of trials');