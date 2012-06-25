function fprintf_wrapper(varargin)

if nargin == 0,
  fprintf('...\n');
elseif nargin == 1,
  fprintf('%s...\n',varargin{1});
else
  s = regexprep(varargin{1},'\\([^\abfnrtvx0-9])','\\\\$1');
  fprintf([s,'...\n'],varargin{2:end});
end