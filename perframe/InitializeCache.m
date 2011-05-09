function cache = InitializeCache()

cache = struct;
cache.mean = struct('radii',{[]},'data',{{}});
cache.std = struct('radii',{[]},'data',{{}});
cache.min = struct('radii',{[]},'data',{{}});
cache.max = struct('radii',{[]},'data',{{}});
