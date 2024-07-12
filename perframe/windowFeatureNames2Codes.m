% [featureCodes,perframeFeatureNames,statNames,transNames] = ...
%   windowFeatureNames2Codes(featureNames,perframeFeatureNames,statNames,transNames)
% 
% Computes a compact code version of window feature names.
% Inputs:
% featureNames: featureNames{i} is a cell array with arguments defining 
% window feature i
% Optional: perframeFeatureNames, statNames, transNames: these are cell
% arrays with the names of per frame features, window statistics, and
% window transform types. Columns of the code are indices into these cells.
% If these arguments are not given, or are empty cells, then they are set
% to an alphabetical list of all per frame features, window statistics, and
% window transforms in featureNames (via the unique function). 
% Outputs:
% featureCodes: featureCodes(i,:) is the code describing window feature i.
% perframeFeatureNames, statNames, transNames: cells defining the lookup codes
% for the names of per frame features, window statistics, and window
% transform types. 
% Note that if new window feature parameter types are added to window
% features, they need to be added here! The only ones I think exist are
% change_window_radius and num_harmonic.
% Added by KB 20201209 to reduce storage size of jab files. 

function [featureCodes,perframeFeatureNames,statNames,transNames] = ...
  windowFeatureNames2Codes(featureNames,perframeFeatureNames,statNames,transNames)

if nargin < 2,
  perframeFeatureNames = {};
end
if nargin < 3,
  statNames = {};
end
if nargin < 4,
  transNames = {};
end

pff_codei = 1;
stat_codei = 2;
trans_codei = 3;
radius_codei = 4;
offset_codei = 5;
changeradius_codei = 6;
num_harmonic_codei = 7;
hist_edges_codei = [6,7];
pff_namei = 1;
stat_namei = 3;
trans_namei = 5;
radius_namei = 7;
offset_namei = 9;
extra_namei = 11;
  
featureCodes = nan(numel(featureNames),7);
pffs = cell(numel(featureNames),1);
stats = cell(numel(featureNames),1);
trans = cell(numel(featureNames),1);
radius = nan(numel(featureNames),1);
offset = nan(numel(featureNames),1);
extra = cell(numel(featureNames),1);
for i = 1:numel(featureNames),
  pffs{i} = featureNames{i}{pff_namei};
  stats{i} = featureNames{i}{stat_namei};
  trans{i} = featureNames{i}{trans_namei};
  radius(i) = featureNames{i}{radius_namei};
  offset(i) = featureNames{i}{offset_namei};
  if numel(featureNames{i}) >= extra_namei,
    % ALT is responsible for the questionable choices made in this block
    % He added them on or around 2021-03-04.
    temp = featureNames{i}{extra_namei} ;
    extra{i} = temp;
    % if ~isscalar(temp) ,
    %   fprintf('In windowFeatureNames2Codes(), featureNames{%d}{%d} is non-scalar:\n', i, extra_namei) ;
    %   disp(temp) ;
    %   fprintf('This will cause an error that is (usually) silently ignored by Matlab if it happens during loading of an object without a loadobj() method.\n') ;
    %   extra(i) = temp(1);
    % else
    %   extra(i) = temp ;
    % end
  end
end

if isempty(perframeFeatureNames),
  [perframeFeatureNames,~,featureCodes(:,pff_codei)] = unique(pffs);
else
  [ism,featureCodes(:,pff_codei)] = ismember(pffs,perframeFeatureNames);
  assert(all(ism));
end
if isempty(statNames),
  [statNames,~,featureCodes(:,stat_codei)] = unique(stats);
else
  [ism,featureCodes(:,stat_codei)] = ismember(pffs,statNames);
  assert(all(ism));
end

if isempty(transNames),
  [transNames,~,featureCodes(:,trans_codei)] = unique(trans);
else
  [ism,featureCodes(:,trans_codei)] = ismember(pffs,transNames);
  assert(all(ism));
end

featureCodes(:,radius_codei) = radius;
featureCodes(:,offset_codei) = offset;
idxextra = find(~cellfun(@isempty,extra));
%idxextra = find(~isnan(extra));
for i = idxextra',
  switch featureNames{i}{extra_namei-1},
    case 'change_window_radius',
      featureCodes(i,changeradius_codei) = extra{i};
    case 'num_harmonic'
      featureCodes(i,num_harmonic_codei) = extra{i};
    case 'hist_edges'
      featureCodes(i,hist_edges_codei) = extra{i};
    otherwise
      error('unknown extra parameter %s',featureNames{i}{extra_namei-1});
  end
end
