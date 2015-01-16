function h = PlotInterpColorLine(x,y,colors,alphas,varargin)

x = x(:);
y = y(:);
x = [x;flipud(x)];
y = [y;flipud(y)];
colors = [colors;flipud(colors)];

if ischar(alphas),
  varargin = [alphas,varargin];
  alphas = [];
elseif numel(alphas) > 1,
  alphas = [alphas;flipud(alphas)];
end
  
if isempty(alphas),
  h = patch(x,y,[0,0,0],'FaceVertexCData',colors,'FaceColor','none','EdgeColor','interp',varargin{:});
elseif numel(alphas) > 1,
  h = patch(x,y,[0,0,0],'FaceVertexCData',colors,'FaceColor','none','EdgeColor','interp','FaceVertexAlphaData',...
    alphas,'EdgeAlpha','interp',varargin{:});
else
  h = patch(x,y,[0,0,0],'FaceVertexCData',colors,'FaceColor','none','EdgeColor','interp','EdgeAlpha',...
    alphas,varargin{:});
end
