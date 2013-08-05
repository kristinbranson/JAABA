% -------------------------------------------------------------------------
function cache_thread(N,HWD,cache_filename,movie_filename)

% lastused = nan means please add framenum to cache
%          = 0 means image is invalid and removed from cache
%          = -1 means it is locked and being added to cache
%          otherwise it is the timestamp the valid image was last used

if isempty(movie_filename),
  return;
end

Mframenum = memmapfile(cache_filename, 'Writable', true, 'Format', 'double', 'Repeat', N);
Mlastused = memmapfile(cache_filename, 'Writable', true, 'Format', 'double', 'Repeat', N, 'Offset', N*8);
Mimage    = memmapfile(cache_filename, 'Writable', true, 'Format', {'uint8' HWD 'x'},  'Repeat', N, ...
    'Offset', 2*N*8);

readframe=get_readframe_fcn(movie_filename);

while true
  idx=find(isnan(Mlastused.Data));
  if(~isempty(idx))
    idx2=argmax(Mframenum.Data(idx));
    Mlastused.Data(idx(idx2))=-1;
    fnum = Mframenum.Data(idx(idx2));
    dd = uint8(readframe(fnum));
    pause(0.0003);     % BJA: why this pause??
    % MK: Cache the read frame to reduce the number of clashes with
    % UpdatePlots
    if Mframenum.Data(idx(idx2)) == fnum
      Mimage.Data(idx(idx2)).x = dd;
      Mlastused.Data(idx(idx2)) = now;
      %Mframenum.Data(idx(idx2)) = fnum;
    end
  else
    pause(0.1);
  end
end

return  %#ok