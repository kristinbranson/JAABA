% JAABA start up script.
<<<<<<< Updated upstream
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

jlabelpath = fileparts(mfilename('fullpath'));
% Initialize all the paths.
baseDir = fileparts(jlabelpath);
addpath(fullfile(baseDir,'misc'));
addpath(fullfile(baseDir,'filehandling'));
addpath(fullfile(jlabelpath,'compute_perframe_features'));
addpath(fullfile(baseDir,'plot'));
addpath(fullfile(baseDir,'structsvm'));

try
  c=parcluster;
  if (c.NumWorkers>2) && (matlabpool('size')<1)
    matlabpool('open',c.NumWorkers-1);  % BJA: must save one for frame cache thread
  end
end
% Start JAABA.
JLabel();
