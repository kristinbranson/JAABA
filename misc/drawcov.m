% h = drawcov(mu,S)
function varargout = drawcov(varargin)
if nargin == 1,
  mu = varargin{1}.centres;
  S = varargin{1}.covars;
elseif nargin == 2,
  mu = varargin{1};
  S = varargin{2};
else,
  error('Usage: drawcov(mix) or drawcov(mu,S)');
end;
K = size(mu,1);

hv = ishold;

hold on;
colors = 'rgb';
h = zeros(1,K);
for k = 1:K,
    if ~any(isnan(S)),
      [a,b,theta] = cov2ell(S(:,:,k));
      h(k) = ellipsedraw(a,b,mu(k,1),mu(k,2),theta,colors(k));
    end;
end;

if ~ishold,
  hold off;
end;

if nargout == 1,
    varargout{1} = h;
end;