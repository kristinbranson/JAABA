% axisalmosttight(border,varargin)

function axisalmosttight(border,varargin)

if nargin < 1 || isempty(border),
  border = 1/20;
end
if numel(varargin)>=1,
  hax = varargin{1};
else
  hax = gca;
end

ax = [inf,-inf,inf,-inf];
for i = 1:numel(hax),
  axis(hax(i),'tight');
  ax1 = axis(hax);
  ax = [min(ax(1),ax1(1)),max(ax(2),ax1(2)),min(ax(3),ax1(3)),max(ax(4),ax1(4))];
end
dx = (ax(2) - ax(1))*border;
ax(1) = ax(1) - dx;
ax(2) = ax(2) + dx;
dy = (ax(4) - ax(3))*border;
ax(3) = ax(3) - dy;
ax(4) = ax(4) + dy;

for i = 1:numel(hax),
  axis(hax(i),ax);
end