function success=testV050JabFileOpening()

jabFileName='/groups/branson/bransonlab/adam/may01ing.jab';
gtMode=false;
data=JLabelData('setstatusfn',@(str)(fprintf('%s\n',str)), ...
                'clearstatusfn',@()(nop()));
data.openJabFile(jabFileName,gtMode);
if ~isfield(data.trxGraphicParams,'nextra_markers')
  error('testV050JabFileOpening:nextra_markersMissing', ...
        'nextra_markers is not a field of data.trxGraphicParams');
end
data.closeJabFile();
data=[];  %#ok
success=true;

end
