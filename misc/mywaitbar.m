%MYWAITBAR Display wait bar.
%   H = MYWAITBAR(X,'message', property, value, property, value, ...)
%   creates and displays a mywaitbar of fractional length X.  The
%   handle to the waitbar figure is returned in H.
%   X should be between 0 and 1.  Optional arguments property and
%   value allow to set corresponding waitbar figure properties.
%   Property can also be an action keyword 'CreateCancelBtn', in
%   which case a cancel button will be added to the figure, and
%   the passed value string will be executed upon clicking on the
%   cancel button or the close figure button.
%
%   MYWAITBAR(X) will set the length of the bar in the most recently
%   created waitbar window to the fractional length X.
%
%   MYWAITBAR(X,H) will set the length of the bar in waitbar H
%   to the fractional length X. If no handle H exists, then a new waitbar
%   figure will be created. 
%
%   MYWAITBAR(X,H,'message') will update the message text in
%   the waitbar figure, in addition to setting the fractional
%   length to X.
%
%   MYWAITBAR is typically used inside a FOR loop that performs a
%   lengthy computation.
%
%   Example:
%       h = waitbar(0,'Please wait...');
%       for i=1:1000,
%           % computation here %
%           waitbar(i/1000,h)
%       end
%
%   See also DIALOG, MSGBOX.

%   Copyright 1984-2010 The MathWorks, Inc.
%   $Revision: 1.23.4.23 $  $Date: 2011/03/09 07:06:28 $

% Generate a warning in -nodisplay and -noFigureWindows mode.

function varargout = mywaitbar(varargin)

if nargin >= 2 && ischar(varargin{end-1}) && ...
    strcmpi(varargin{end-1},'interpreter') && ...
    ischar(varargin{end}),
  interpreter = varargin{end};
  leftovers = varargin(1:end-2);
else
  interpreter = 0;
  leftovers = varargin;
end

if numel(leftovers) >= 2 && isnumeric(leftovers{2}),
  if ~ishandle(leftovers{2}),
    leftovers(2) = [];
  end
end

if ischar(interpreter),
  if numel(leftovers) >= 2 && isnumeric(leftovers{2}) && ~isempty(leftovers{2}),
    htitle = get(get(leftovers{2},'Children'),'Title');
    set(htitle,'Interpreter',interpreter);
  elseif ischar(leftovers{end}) && strcmpi(interpreter,'none'),
    ti0 = leftovers{end};
    ti = ti0;
    ti = regexprep(ti,'\\','\\\\');
    ti = regexprep(ti,'_','\_');
    leftovers{end} = ti;
    hwait = waitbar(leftovers{:});
    htitle = get(get(hwait,'Children'),'title');
    set(htitle,'Interpreter',interpreter,'String',ti0);
    drawnow;
    if nargout >= 1,
      varargout = {hwait};
      return;
    end
  end
end
varargout = cell(1,nargout);
[varargout{:}] = waitbar(leftovers{:});
