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

%SetUpJAABAPath;

try
  if ishandle(hstatustext),
    set(hstatustext,'String','Starting parallel computing workers...');
  end
  if isdeployed,
    if ispc || (isunix && ~ismac),
      c = parcluster;
      fprintf('Current profile: %s\n',c.Profile);
      if ~strcmp(c.Profile,'JAABAParCompProfile'),
        filename = deployedRelative2Global('JAABAParCompProfile.settings');
        if ~exist(filename,'file'),
          fprintf('Could not find file %s, not using parallel computing.\n',filename);
        else
          setmcruserdata('ParallelProfile',filename);
        end
      end
    end
  end
  if verLessThan('matlab','8.3.0.532'),
    nthreads = SetUpMatlabPool;
  else
    nthreads = SetUpParPool;
  end
catch ME,
  uiwait(warndlg(sprintf('Error starting parallel computing: %s',getReport(ME))));
  nthreads = struct;
end

% Start JAABA.

%if ismac,
%  cd(getenv('JAABA_RUNDIR'));
%else
%  warndlg(sprintf('Starting JAABA in %s, ctfroot is %s, mfilename is %s',pwd,ctfroot,mfilename('fullpath')));
%end
try
  if ishandle(hstatustext),
    set(hstatustext,'String','Starting JAABA...');
    args = {'hsplash',hsplash,'hsplashstatus',hstatustext,'nthreads',nthreads};
  else
    args = {};
  end
  JLabelHandle = JLabel(args{:});
  if ishandle(JLabelHandle)
    uiwait(JLabelHandle);
  end
catch ME,
  uiwait(warndlg(getReport(ME)));
  delete(findall(0,'type','figure'));
end

try %#ok<TRYNC>
  if matlabpool('size')>=1
    matlabpool close;
  end
end

if isdeployed,
  delete(findall(0,'type','figure'));
end