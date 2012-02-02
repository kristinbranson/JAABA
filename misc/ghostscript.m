function varargout = ghostscript(cmd)
%GHOSTSCRIPT  Calls a local GhostScript executable with the input command
%
% Example:
%   [status result] = ghostscript(cmd)
%
% Attempts to locate a ghostscript executable, finally asking the user to
% specify the directory ghostcript was installed into. The resulting path
% is stored for future reference.
% 
% Once found, the executable is called with the input command string.
%
% This function requires that you have Ghostscript installed on your
% system. You can download this from: http://www.ghostscript.com
%
% IN:
%   cmd - Command string to be passed into ghostscript.
%
% OUT:
%   status - 0 iff command ran without problem.
%   result - Output from ghostscript.

% Copyright: Oliver Woodford, 2009-2010

% Thanks to Jonas Dorn for the fix for the title of the uigetdir window on
% Mac OS.

% Thanks to Nathan Childress for the fix to the default location on 64-bit
% Windows systems.

% Call ghostscript
[varargout{1:nargout}] = system(sprintf('"%s" %s', gs_path, cmd));
return

function path = gs_path
% Return a valid path
% Start with the currently set path
path = current_gs_path;
% Check the path works
if check_gs_path(path)
    return
end
% Check whether the binary is on the path
if ispc
    bin = 'gswin32c.exe';
else
    bin = 'gs';
end
if check_store_gs_path(bin)
    path = bin;
    return
end
% Search the obvious places
if ispc
    default_location = 'C:\Program Files\gs\';
    dir_list = dir(default_location);
    if isempty(dir_list)
        default_location = 'C:\Program Files (x86)\gs\'; % Possible location on 64-bit systems 
        dir_list = dir(default_location);
    end
    executable = '\bin\gswin32c.exe';
    ver_num = 0;
    % If there are multiple versions, use the newest
    for a = 1:numel(dir_list)
        ver_num2 = sscanf(dir_list(a).name, 'gs%g');
        if ~isempty(ver_num2) && ver_num2 > ver_num
            path2 = [default_location dir_list(a).name executable];
            if exist(path2, 'file') == 2
                path = path2;
                ver_num = ver_num2;
            end
        end
    end
else
    path = '/usr/local/bin/gs';
end
if check_store_gs_path(path)
    return
end
% Ask the user to enter the path
while 1
    if strncmp(computer, 'MAC', 3) % Is a Mac
        % Give separate warning as the uigetdir dialogue box doesn't have a
        % title
        uiwait(warndlg('Ghostscript not found. Please locate the program.'))
    end
    base = uigetdir('/', 'Ghostcript not found. Please locate the program.');
    if isequal(base, 0)
        % User hit cancel or closed window
        break;
    end
    base = [base filesep];
    bin_dir = {'', ['bin' filesep], ['lib' filesep]};
    for a = 1:numel(bin_dir)
        path = [base bin_dir{a} bin];
        if exist(path, 'file') == 2
            break;
        end
    end
    if check_store_gs_path(path)
        return
    end
end
error('Ghostscript not found.');

function good = check_store_gs_path(path)
% Check the path is valid
good = check_gs_path(path);
if ~good
    return
end
% Update the current default path to the path found
if change_value(path, 'current_gs_path_str', [mfilename('fullpath') '.m'])
    warning('Path to ghostscript installation could not be saved. Enter it manually in ghostscript.m.');
    return
end
return

function good = check_gs_path(path)
% Check the path is valid
[good message] = system(sprintf('"%s" -h', path));
good = good == 0;
return

function current_gs_path_str = current_gs_path
current_gs_path_str = 'C:\Program Files\gs\gs8.71\bin\gswin32c.exe';
return