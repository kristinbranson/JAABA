% JAABA start up script.

% Initialize all the paths.
% try
%   if matlabpool('size')<1
%     matlabpool open;
%   end
% catch ME
%   warndlg(sprintf('Could not start matlabpool with just "matlabpool open": %s',getReport(ME)));
% end

[hsplash,hstatustext] = JAABASplashScreen();
% 
% fprintf('Current directory: %s\n',pwd);
% fprintf('JAABA location: %s\n',mfilename('fullpath'));
% fprintf('ctfroot -> %s\n',ctfroot);

try
  if ishandle(hstatustext),
    set(hstatustext,'String','Starting parallel computing workers...');
  end
  if isdeployed,
    if ispc,
      setmcruserdata('ParallelProfile','ParallelComputingConfiguration_Local_Win4.settings');
    end
  end
  if matlabpool('size') < 1,
    matlabpool('open');
  end
catch ME,
  uiwait(warndlg('Error starting parallel computing: %s',getReport(ME)));
end

% Start JAABA.

if ismac,
  cd(getenv('JAABA_RUNDIR'));
%else
%  warndlg(sprintf('Starting JAABA in %s, ctfroot is %s, mfilename is %s',pwd,ctfroot,mfilename('fullpath')));
end
try
  if ishandle(hstatustext),
    set(hstatustext,'String','Starting JAABA...');
    args = {'hsplash',hsplash,'hsplashstatus',hstatustext};
  else
    args = {};
  end
  uiwait(JLabel(args{:}));
catch ME,
  uiwait(warndlg(getReport(ME)));
end

try %#ok<TRYNC>
  if matlabpool('size')>=1
    matlabpool close;
  end
end