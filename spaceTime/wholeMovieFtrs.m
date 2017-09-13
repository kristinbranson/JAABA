function ftrs = wholeMovieFtrs(moviename,trackfilename,blockSize)

tracks = load_tracks(trackfilename);
[readfcn,nframes,fid,headerinfo] = get_readframe_fcn(moviename);

if fid>0,
  fclose(fid);
end

nblocks = ceil(nframes/blockSize);

parfor ndx = 1:nblocks
  [readfcn,nframes,fid,headerinfo] = get_readframe_fcn(moviename);
  % For parfor files has to be opened for every loop.
  blockStart = (ndx-1)*blockSize + 1;
  blockEnd   = min(ndx*blockSize,nframes);
  ftrs{ndx} = genFeatures(readfcn,headerinfo,blockStart,blockEnd,tracks);
  if fid>0,
    fclose(fid);
  end
  fprintf('.');
  if mod(ndx,50)==0, fprintf('\n'); end
end

