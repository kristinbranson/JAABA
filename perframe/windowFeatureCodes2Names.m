% featureNames =
% windowFeatureCodes2Names(featureCodes,perframeFeatureNames,statNames,transNames)  
%
% Computes the full window feature names from their compact code version. 
% Inputs:
% featureCodes: featureCodes(i,:) is the code describing window feature i.
% perframeFeatureNames, statNames, transNames: cells defining the lookup codes
% for the names of per frame features, window statistics, and window
% transform types. 
% Output:
% featureNames: featureNames{i} is a cell array with arguments defining 
% window feature i
% Note that if new window feature parameter types are added to window
% features, they need to be added here! The only ones I think exist are
% change_window_radius and num_harmonic.
% Added by KB 20201209 to reduce storage size of jab files. 

function featureNames = windowFeatureCodes2Names(featureCodes,perframeFeatureNames,statNames,transNames)

pff_codei = 1;
stat_codei = 2;
trans_codei = 3;
radius_codei = 4;
offset_codei = 5;
changeradius_codei = 6;
num_harmonic_codei = 7;
pff_namei = 1;
stat_namei = 3;
trans_namei = 5;
radius_namei = 7;
offset_namei = 9;
extra_namei = 11;

nfeatures = size(featureCodes,1);
featureNames = cell(1,nfeatures);

base = cell(1,extra_namei);
base{stat_namei-1} = 'stat';
base{trans_namei-1} = 'trans';
base{radius_namei-1} = 'radius';
base{offset_namei-1} = 'offset';
base{extra_namei-1} = 'extra';
for i = 1:nfeatures,
  featureCode = featureCodes(i,:);
  featureName = base;
  if ~isnan(featureCode(changeradius_codei)),
    featureName{extra_namei-1} = 'change_window_radius';
    featureName{extra_namei} = featureCode(changeradius_codei);
  elseif ~isnan(featureCode(num_harmonic_codei)),
    featureName{extra_namei-1} = 'num_harmonic';
    featureName{extra_namei} = featureCode(num_harmonic_codei);
  else
    featureName = base(1:offset_namei);
  end
  featureName{pff_namei} = perframeFeatureNames{featureCode(pff_codei)};
  featureName{stat_namei} = statNames{featureCode(stat_codei)};
  featureName{trans_namei} = transNames{featureCode(trans_codei)};
  featureName{radius_namei} = featureCode(radius_codei);
  featureName{offset_namei} = featureCode(offset_codei);
  featureNames{i} = featureName;
end
      
