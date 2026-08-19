function success=testMP4ReadNframes()
% Test that mp4_read_nframes() reads the exact frame count from MP4 container
% metadata, for files with the moov box at the end (the common case, and how
% the cameras in the lab write files) and at the front (faststart).  The
% metadata count is compared against a known-good frame count for each test
% movie, so we avoid VideoReader's authoritative-but-slow NumberOfFrames, which
% is the count that mp4_read_nframes() exists to avoid computing.

  cases = struct('moviePath', ...
                   { '/groups/branson/bransonlab/projects/JAABA/test_data/mp4/mp4_nframes_test_moov_at_end.mp4', ...
                     '/groups/branson/bransonlab/projects/JAABA/test_data/mp4/mp4_nframes_test_moov_at_front.mp4' }, ...
                 'expectedFrameCount', ...
                   { 137, ...
                     89 }) ;

  for caseIndex = 1 : numel(cases) ,
    moviePath = cases(caseIndex).moviePath ;
    expectedFrameCount = cases(caseIndex).expectedFrameCount ;

    metadataFrameCount = mp4_read_nframes(moviePath) ;
    if ~isequal(metadataFrameCount, expectedFrameCount) ,
      error('mp4_read_nframes returned %d frames for %s, expected %d', ...
            metadataFrameCount, moviePath, expectedFrameCount) ;
    end

    % get_readframe_fcn() should return the same frame count (via the metadata
    % fast path) and a headerinfo struct populated with the movie's metadata,
    % without having called VideoReader's slow NumberOfFrames.
    [readframe, readframeFrameCount, fid, headerinfo] = get_readframe_fcn(moviePath) ;
    cleaner = onCleanup(@()(closeIfOpen(fid))) ;
    if ~isequal(readframeFrameCount, expectedFrameCount) ,
      error('get_readframe_fcn returned %d frames for %s, expected %d', ...
            readframeFrameCount, moviePath, expectedFrameCount) ;
    end
    requiredFields = { 'FrameRate', 'VideoFormat', 'nr', 'nc', 'nframes' } ;
    for fieldIndex = 1 : numel(requiredFields) ,
      if ~isfield(headerinfo, requiredFields{fieldIndex}) ,
        error('get_readframe_fcn headerinfo for %s is missing field %s', ...
              moviePath, requiredFields{fieldIndex}) ;
      end
    end
    firstImage = readframe(1) ;
    if size(firstImage, 1) ~= headerinfo.nr || size(firstImage, 2) ~= headerinfo.nc ,
      error('get_readframe_fcn frame size does not match headerinfo for %s', moviePath) ;
    end
    clear cleaner ;
  end

  success = true ;

end  % function

function closeIfOpen(fid)
% Close the file identifier fid if it refers to an open file.
  if fid > 0 ,
    fclose(fid) ;
  end
end  % function
