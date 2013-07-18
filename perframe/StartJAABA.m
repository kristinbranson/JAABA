% JAABA start up script.
%
% This program is part of JAABA.
%
% JAABA: The Janelia Automatic Animal Behavior Annotator
% Copyright 2012, Kristin Branson, HHMI Janelia Farm Resarch Campus
% http://jaaba.sourceforge.net/
% bransonk@janelia.hhmi.org
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License (version 3 pasted in LICENSE.txt) for 
% more details.

SetUpJAABAPath;

try
  c=parcluster;
  if (c.NumWorkers>2) && (matlabpool('size')<1)
    % BJA: as of matlab 2013a 12 workers can be used even on a 1-core machine
    if c.NumWorkers < min(12, 2 * feature('numCores'))
      disp(['WARNING: for best performance set NumWorkers in parallel prefs / cluster profile to be at least min(12, 2 * #cores) = ' num2str(min(12, 2 * feature('numCores'))) ' on this machine']);
    end
    rc=load('.JLabelrc.mat');  % would've been nice to be able to just call LoadRC
    if(~isfield(rc,'framecache_threads'))
      rc.framecache_threads=1;
    end
    parfor_threads = min(c.NumWorkers - rc.framecache_threads, feature('numCores'));
    matlabpool('open',parfor_threads);
  end
end
% Start JAABA.
JLabel();
