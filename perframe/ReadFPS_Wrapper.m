function [success,msg,fps] = ...
  ReadFPS_Wrapper(InputDataType,varargin)

success = false;
msg = '';

[fps,leftovers] = myparse_nocheck(varargin,'fps',30);

switch InputDataType,
  case 'Ctrax',

    [success,msg,fps] = ...
      ReadFPS_Ctrax(...
      leftovers{:},...
      'fps',fps);

  case 'LarvaeRiveraAlba',
    
    success = false;
    msg = 'Not implemented';

  case 'MoTr',

    success = false;
    msg = 'FPS not stored in MoTr files';
    
  case 'Qtrax',

    [success,msg,fps] = ...
      ReadFPS_Qtrax(...
      leftovers{:},...
      'fps',fps);
    
  case 'MAGATAnalyzer',

    [success,msg,fps] = ...
      ReadFPS_MAGATAnalyzer(...
      leftovers{:},...
      'fps',fps);
    
  case 'MWT',
    
    [success,msg,fps] = ...
      ReadFPS_MWT(...
      leftovers{:},...
      'fps',fps);

  otherwise
    success = false;
    msg = sprintf('Unknown data type %s',InputDataType);
end