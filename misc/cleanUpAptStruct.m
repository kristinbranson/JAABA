function result = cleanUpAptStruct(aptStructRaw)
  % Tidy up some aspects of aptStructRaw.
  % Currently, replaces path names with just (leaf) file names in .trkfilename
  % field, if it exists.
  % Returns the modified struct in result.
  result = aptStructRaw ;
  if ~isfield(aptStructRaw, 'trkfilename') ,
    return
  end
  rawtrkfilename = aptStructRaw.trkfilename ;
  % Replace full paths with just (leaf) file name.
  % I think rawtrkfilename is always a cellstr, but handle the case where its a
  % string anyway.
  if iscell(rawtrkfilename) 
    trkfilename = cellfun(@fileparts23, rawtrkfilename, 'UniformOutput', false) ;
  else
    trkfilename = fileparts23(rawtrkfilename) ;
  end
  result.trkfilename = trkfilename ;
end
