function result= ...
  wfParamsOfOneType(wfParamsThisTypeFromWindowFeatureParams, ...
                    fallbackValues, ...
                    wfParamNames,extraParamName)

%wfParamsThisTypeFromWindowFeatureParams = windowFeatureParams.(pfName).(wfType);
result.enabled = true;
for i = 1:length(wfParamNames)
  wfParamName = wfParamNames{i};
  if isfield(wfParamsThisTypeFromWindowFeatureParams,wfParamName)
    result.values.(wfParamName) = wfParamsThisTypeFromWindowFeatureParams.(wfParamName);
  else  % fill in the default values
    % Is this really where we want to get it from?
    % 
    result.values.(wfParamName) = fallbackValues.(wfParamName);
  end
end
if ~isempty(extraParamName)
  %extraParamName = self.wfExtraParamNames{wfTypeNdx};
  result.values.(extraParamName) = wfParamsThisTypeFromWindowFeatureParams.(extraParamName);
end

end
