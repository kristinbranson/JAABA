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

  case 'MouseHouse',

    success = false;
    msg = 'FPS not stored in MouseHouse files';
    
  case 'Qtrax',

    [success,msg,fps] = ...
      ReadFPS_Qtrax(...
      leftovers{:},...
      'fps',fps);
    
  case 'MAGATAnalyzer',

    success = false;
    msg = 'Not implemented';
    
  case 'MWT',
    
    success = false;
    msg = 'Not implemented';    

  otherwise
    success = false;
    msg = sprintf('Unknown data type %s',InputDataType);
end