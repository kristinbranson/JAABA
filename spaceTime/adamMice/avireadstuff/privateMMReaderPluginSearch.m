function [ pluginList, errFlag ] = privateMMReaderPluginSearch()
%MMREADERPLUGINSEARCH Locate mmreader plugins
%   [PLUGINLIST] = PRIVATEMMREADERPLUGINSEARCH locates all available
%   mmreader plugin files and returns the fully qualified path to each in a
%   cell array PLUGINLIST.  
%
%   MMReader plugins are searched for in the mmreader plugins directory:
%
%       (matlabroot)/toolbox/matlab/audiovideo/bin/mmreader/plugins/(arch)
%
%   where (arch) is a given platform (maci, win32, etc).
%
%   PRIVATEMMREADERPLUGINSEARCH is used internally by mmreader and is not
%   intended to be used directly by an end user.

%   NH 
%   Copyright 2007-2009 The MathWorks, Inc.
%
%   $Revision: 1.1.6.6 $ $Date: 2011/07/20 00:00:22 $


% Initailize variables
pluginList = {};
errFlag = false;

osDir = computer('arch');
osExt = feature('GetSharedLibExt');

% Define the mmreader plugin path
pluginDir = fullfile(matlabroot, ...
    'toolbox', 'shared', 'multimedia', 'bin', osDir,'reader');
wildFile = ['*' osExt];

% Perform a wildecard search (i.e. *.dll) 
% for any mmreader adaptors

searchPath = fullfile(pluginDir, wildFile );
dirList = dir(searchPath);

for ii=1:length(dirList)
    [~, pluginFileName] = fileparts( dirList(ii).name );
    pluginList = {pluginList{:} fullfile(pluginDir, pluginFileName)}; %#ok<CCAT>
end


end