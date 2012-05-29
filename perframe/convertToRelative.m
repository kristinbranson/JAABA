function modX = convertToRelative(x,relativeParams)

bins = relativeParams.relativeBins;
while( any((bins(1:end-1)-bins(2:end))>0))
    ndx = find( (bins(1:end-1)-bins(2:end))>0);
    bins(ndx+1) = bins(ndx);
end

[~,modX] = histc(x,bins);
modX(modX>numel(bins))=numel(bins);

validX = ~isnan(modX) & (modX<length(bins)) & (modX>0);
extra = x(validX)-bins(modX(validX));
relExtra = extra./(bins(modX(validX)+1)-bins(modX(validX)));
modX(validX) = modX(validX) + relExtra;
modX = modX-1; % To make it go from 0 instead of 1