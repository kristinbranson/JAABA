% h = obj.drawfly(fly,t,varargin)
% draw shape corresponding to fly position
% optional arguments:
% 'shape': either 'triangle', 'ellipse', or 'point' (default = 'triangle'). 
% 'filled': whether the shape will be filled in or not (default = false)
% 'registered': whether to plot the registered position (mm) or the raw
% pixel position
% all other extra arguments will be fed to plot/patch
function varargout = drawfly(obj,fly,t,varargin)

[shape,filled,registered,leftovers] = ...
  myparse_nocheck(varargin,...
  'shape','triangle','filled',false,'registered',false);

if registered,
  x = obj.GetPerFrameData('x_mm',fly);
  y = obj.GetPerFrameData('y_mm',fly);
  theta = obj.GetPerFrameData('theta_mm',fly);
  a = obj.GetPerFrameData('a_mm',fly);
  b = obj.GetPerFrameData('b_mm',fly);
else
  x = obj.GetPerFrameData('x',fly);
  y = obj.GetPerFrameData('y',fly);
  theta = obj.GetPerFrameData('theta',fly);
  a = obj.GetPerFrameData('a',fly);
  b = obj.GetPerFrameData('b',fly);
end

x = x(t);
y = y(t);
a = a(t);
b = b(t);
theta = theta(t);

if strcmpi(shape,'ellipse') && filled,

  h = ellipsedrawpatch(a*2,b*2,x,y,theta,leftovers{:});
  
elseif strcmpi(shape,'ellipse') && ~filled,

  h = ellipsedraw(a*2,b*2,x,y,theta,'-',leftovers{:});

elseif strcmpi(shape,'triangle'),
  % isosceles triangle not yet rotated or centered
  pts = [-a*2,-b*2
    -a*2,b*2
    a*2,0];
  
  % rotate
  costheta = cos(theta);
  sintheta = sin(theta);
  R = [costheta,sintheta;-sintheta,costheta];
  pts = pts*R;
  
  % translate
  pts(:,1) = pts(:,1) + x;
  pts(:,2) = pts(:,2) + y;

  if filled,
    h = patch(pts([1:3,1],1),pts([1:3,1],2),'b',leftovers{:});
  else
    h = plot(pts([1:3,1],1),pts([1:3,1],2),leftovers{:});
  end
    
elseif strcmpi(shape,'point'),
  h = plot(x,y,'.',leftovers{:});
end

if nargout > 0,
  varargout{1} = h;
end

