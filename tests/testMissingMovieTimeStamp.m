function success=testMissingMovieTimeStamp()

jabFileName='/groups/branson/bransonlab/adam/larva_mwt_rolling.jab';
gtMode=false;
data=JLabelData('setstatusfn',@(str)(fprintf('%s\n',str)), ...
                'clearstatusfn',@()(nop()));
data.openJabFile(jabFileName,gtMode);
fileTypes=data.filetypes;
iFileType=whichstr('movie',fileTypes);
nExps=data.nexps;
for iExp=1:nExps
  if data.fileexists(iExp,iFileType) ,
    error('testMissingMovieTimeStamp:fileexistsSetWrong', ...
          'filexists is true where it should be false');
  end
  if ~isinf(data.filetimestamps(iExp,iFileType)) ,
    error('testMissingMovieTimeStamp:filetimestampSetWrong', ...
          'filetimestamps is finite where it should be -inf');
  end
end
data.closeJabFile();
data=[];  %#ok
success=true;

end
