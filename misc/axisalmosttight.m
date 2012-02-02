% axisalmosttight(border,varargin)

function axisalmosttight(border,varargin)

if nargin < 1 || isempty(border),
  border = 1/20;
end

axis(varargin{:},'tight');
ax = axis(varargin{:});
dx = (ax(2) - ax(1))*border;
ax(1) = ax(1) - dx;
ax(2) = ax(2) + dx;
dy = (ax(4) - ax(3))*border;
ax(3) = ax(3) - dy;
ax(4) = ax(4) + dy;

axis(varargin{:},ax);