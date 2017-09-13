function hh = errorbar_x(x, y, l,u,symbol)
%ERRORBAR Error bar plot.
%   ERRORBAR_X(X,Y,L,U) plots the graph of vector X vs. vector Y with 
%   error bars specified by the vectors L and U in horizontal direction.  
%   L and U contain the lower and upper error ranges for each point 
%   in X (lower = left side, upper = right side).  Each error bar is 
%   L(i) + U(i) long and is drawn a distance of U(i) from the right 
%   and L(i) from the left the points in (X,Y).  The vectors X,Y,L 
%   and U must all be the same length.  If X,Y,L and U are matrices 
%   then each column produces a separate line.
%
%   ERRORBAR_X(X,Y,E) or ERRORBAR(Y,E) plots X with error bars [X-E X+E].
%   ERRORBAR_X(...,'LineSpec') uses the color and linestyle specified by
%   the string 'LineSpec'.  See PLOT for possibilities.
%
%   H = ERRORBAR_X(...) returns a vector of line handles.
%
%   For example,
%      x = 1:10;
%      y = sin(x);
%      e = std(y)*ones(size(x));
%      errorbar_x(x,y,e)
%   draws symmetric error bars of unit standard deviation.

%   L. Shure 5-17-88, 10-1-91 B.A. Jones 4-5-93
%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 5.19 $  $Date: 2002/06/05 20:05:14 $

%   modified for plotting error bars in a logarithmic graph by:
%   Goetz Huesken
%   e-mail: goetz.huesken(at)gmx.de
%   Date: 10/23/2006


if min(size(x))==1,
  npt = length(x);
  x = x(:);
  y = y(:);
    if nargin > 2,
        if ~isstr(l),  
            l = l(:);
        end
        if nargin > 3
            if ~isstr(u)
                u = u(:);
            end
        end
    end
else
  [npt,n] = size(x);
end

if nargin == 3
    if ~isstr(l)  
        u = l;
        symbol = '-';
    else
        symbol = l;
        l = y;
        u = y;
        y = x;
        [m,n] = size(y);
        x(:) = (1:npt)'*ones(1,n);;
    end
end

if nargin == 4
    if isstr(u),    
        symbol = u;
        u = l;
    else
        symbol = '-';
    end
end


if nargin == 2
    l = y;
    u = y;
    y = x;
    [m,n] = size(y);
    x(:) = (1:npt)'*ones(1,n);;
    symbol = '-';
end

u = abs(u);
l = abs(l);
    
if isstr(x) | isstr(y) | isstr(u) | isstr(l)
    error('Arguments must be numeric.')
end

if ~isequal(size(x),size(y)) | ~isequal(size(x),size(l)) | ~isequal(size(x),size(u)),
  error('The sizes of X, Y, L and U must be the same.');
end

m = size(y,1);                      % modification for plotting error bars in x-direction
if m == 1                           %
  tee = abs(y)/40;                  % 
else                                %
  tee = (max(y(:))-min(y(:)))/40;   % 
end                                 %
                                    %
xl = x - l;                         %
xr = x + u;                         %
ytop = y + tee;                     %
ybot = y - tee;                     %
n = size(y,2); 

% Plot graph and bars
hold_state = ishold;
cax = newplot;
next = lower(get(cax,'NextPlot'));

% build up nan-separated vector for bars
xb = zeros(npt*9,n);    % modification for plotting error bars in in x-direction
xb(1:9:end,:) = xr;     %
xb(2:9:end,:) = xl;     %
xb(3:9:end,:) = NaN;    %
xb(4:9:end,:) = xr;     %
xb(5:9:end,:) = xr;     %
xb(6:9:end,:) = NaN;    %
xb(7:9:end,:) = xl;     %
xb(8:9:end,:) = xl;     %
xb(9:9:end,:) = NaN;    %

yb = zeros(npt*9,n);    % modification for plotting error bars in in x-direction
yb(1:9:end,:) = y;      %
yb(2:9:end,:) = y;      %
yb(3:9:end,:) = NaN;    %
yb(4:9:end,:) = ytop;   %
yb(5:9:end,:) = ybot;   %
yb(6:9:end,:) = NaN;    %
yb(7:9:end,:) = ytop;   %
yb(8:9:end,:) = ybot;   %
yb(9:9:end,:) = NaN;    %

[ls,col,mark,msg] = colstyle(symbol); if ~isempty(msg), error(msg); end
symbol = [ls mark col]; % Use marker only on data part
esymbol = ['-' col]; % Make sure bars are solid

h = plot(xb,yb,esymbol); hold on
h = [h;plot(x,y,symbol)]; 

if ~hold_state, hold off; end

if nargout>0, hh = h; end
