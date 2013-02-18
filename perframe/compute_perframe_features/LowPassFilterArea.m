function areasmooth = LowPassFilterArea(area,filterorder,maxfreq)

f = fdesign.lowpass('N,F3db',filterorder,maxfreq);
h = design(f,'butter');
h.PersistentMemory = true;

h.filter(fliplr(area));
areasmooth = h.filter(area);
