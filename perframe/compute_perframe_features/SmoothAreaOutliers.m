function areasmooth = SmoothAreaOutliers(area,filterorder,maxfreq,maxerr)

areasmooth = LowPassFilterArea(area,filterorder,maxfreq);
isoutlier = isinf(area) | isnan(area) | abs(areasmooth - area) > maxerr;
[starts,ends] = get_interval_ends(isoutlier);
ends = ends - 1;
areasmooth = area;
nframes = numel(area);
for i = 1:numel(starts),
  if starts(i) == 1 && ends(i) == nframes,
    break;
  elseif starts(i) == 1,
    areasmooth(starts(i):ends(i)) = area(ends(i)+1);
  elseif ends(i) == nframes,
    areasmooth(starts(i):ends(i)) = area(starts(i)-1);
  else
    areasmooth(starts(i):ends(i)) = (area(starts(i)-1)+area(ends(i)+1))/2;
  end
end
