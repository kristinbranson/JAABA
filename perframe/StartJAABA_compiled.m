% JAABA start up script.

% Initialize all the paths.
if matlabpool('size')<1
  matlabpool open;
end
% Start JAABA.

if ismac,
  cd(getenv('JAABA_RUNDIR'));
%else
%  warndlg(sprintf('Starting JAABA in %s, ctfroot is %s, mfilename is %s',pwd,ctfroot,mfilename('fullpath')));
end
try
  JLabel();
catch ME,
  uiwait(warndlg(getReport(ME)));
end
