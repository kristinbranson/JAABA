% JAABA start up script.

% Initialize all the paths.
if matlabpool('size')<1
  matlabpool open;
end
% Start JAABA.
cd(getenv('JAABA_RUNDIR'));
try
  JLabel();
catch ME,
  uiwait(warndlg(getReport(ME)));
end
