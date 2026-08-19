function success=testMOVReading()
% Test that reading frames from a .mov movie gives the same result regardless
% of the order in which frames are requested.  get_readframe_fcn reads .mov via
% the same fast keyframe-based seek path used for .mp4/.m4v, so a frame's
% content must not depend on what was read before it.  We read 10 evenly-spaced
% frames in ascending order, then read the same 10 frames in a fixed
% randomly-shuffled order, and require each shuffled read to be identical to the
% corresponding in-order read.

  % Fix the RNG so the shuffle (and thus the test) is deterministic.
  rng(42) ;

  moviePath = '/groups/branson/bransonlab/projects/JAABA/test_data/larvae/OregonRA/movie.mov' ;

  [readframe, nframes, fid] = get_readframe_fcn(moviePath) ;
  cleaner = onCleanup(@()(closeIfOpen(fid))) ;

  frameCount = 10 ;
  frameIndices = round(linspace(1, nframes, frameCount)) ;

  % Read the frames in ascending index order.
  inOrderFrames = cell(1, frameCount) ;
  for i = 1 : frameCount ,
    inOrderFrames{i} = readframe(frameIndices(i)) ;
  end

  % Read the same frames in a shuffled order.
  permutation = randperm(frameCount) ;
  shuffledFrames = cell(1, frameCount) ;
  for i = 1 : frameCount ,
    shuffledFrames{i} = readframe(frameIndices(permutation(i))) ;
  end

  % Each shuffled read must match the in-order read of the same frame index.
  for i = 1 : frameCount ,
    if ~isequal(shuffledFrames{i}, inOrderFrames{permutation(i)}) ,
      error('MOV frame %d read out of order differs from the same frame read in order', ...
            frameIndices(permutation(i))) ;
    end
  end

  success = true ;

end  % function

function closeIfOpen(fid)
% Close the file identifier fid if it refers to an open file.
  if fid > 0 ,
    fclose(fid) ;
  end
end  % function
