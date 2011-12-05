function cache = InitializeCache()

cache = struct;
cache.mean = struct('radii',{[]},'data',{{}});
cache.meanRel = struct('radii',{[]},'data',{{}});
cache.std = struct('radii',{[]},'data',{{}});
cache.stdRel = struct('radii',{[]},'data',{{}});
cache.min = struct('radii',{[]},'data',{{}},...
  'idx',{struct('orig',{},'abs',{},'flip',{},'rel',{})});
cache.max = struct('radii',{[]},'data',{{}},...
  'idx',{struct('orig',{},'abs',{},'flip',{},'rel',{})});
cache.relX = [];
