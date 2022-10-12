function handles = plot_arrow( x1,y1,x2,y2,varargin )
%
% plot_arrow - plots an arrow to the current plot
%
% format:   handles = plot_arrow( x1,y1,x2,y2 [,options...] )
%
% input:    x1,y1   - starting point
%           x2,y2   - end point
%           options - come as pairs of "property","value" as defined for "line" and "patch"
%                     controls, see matlab help for listing of these properties.
%                     note that not all properties where added, one might add them at the end of this file.
%                     
%                     additional options are:
%                     'headwidth':  relative to complete arrow size, default value is 0.07
%                     'headheight': relative to complete arrow size, default value is 0.15
%                     (encoded are maximal values if pixels, for the case that the arrow is very long)
%
% output:   handles - handles of the graphical elements building the arrow
%
% Example:  plot_arrow( -1,-1,15,12,'linewidth',2,'color',[0.5 0.5 0.5],'facecolor',[0.5 0.5 0.5] );
%           plot_arrow( 0,0,5,4,'linewidth',2,'headwidth',0.25,'headheight',0.33 );
%           plot_arrow;   % will launch demo
% 
% Copyright (c) 2003, Ohad Gal
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.


% =============================================
% for debug - demo - can be erased
% =============================================
if (nargin==0)
    figure;
    axis;
    set( gca,'nextplot','add' );
    for x = 0:0.3:2*pi
        color = [rand rand rand];
        h = plot_arrow( 1,1,50*rand*cos(x),50*rand*sin(x),...
            'color',color,'facecolor',color,'edgecolor',color );
        set( h,'linewidth',2 );
    end
    hold off;
    return
end
% =============================================
% end of for debug
% =============================================


% =============================================
% constants (can be edited)
% =============================================
alpha       = 0.15;   % head length
beta        = 0.07;   % head width
max_length  = 22;
max_width   = 10;

% =============================================
% check if head properties are given
% =============================================
% if ratio is always fixed, this section can be removed!
if ~isempty( varargin )
    for c = 1:floor(length(varargin)/2)
        try
            switch lower(varargin{c*2-1})
                % head properties - do nothing, since handled above already
            case 'headheight',alpha = max( min( varargin{c*2},1 ),0.01 );
            case 'headwidth', beta = max( min( varargin{c*2},1 ),0.01 );
            end
        catch
            fprintf( 'unrecognized property or value for: %s\n',varargin{c*2-1} );
        end
    end
end

% =============================================
% calculate the arrow head coordinates
% =============================================
den         = x2 - x1 + eps;                                % make sure no devision by zero occurs
teta        = atan( (y2-y1)/den ) + pi*(x2<x1) - pi/2;      % angle of arrow
cs          = cos(teta);                                    % rotation matrix
ss          = sin(teta);
R           = [cs -ss;ss cs];
line_length = sqrt( (y2-y1)^2 + (x2-x1)^2 );                % sizes
head_length = min( line_length*alpha,max_length );
head_width  = min( line_length*beta,max_length );
x0          = x2*cs + y2*ss;                                % build head coordinats
y0          = -x2*ss + y2*cs;
coords      = R*[x0 x0+head_width/2 x0-head_width/2; y0 y0-head_length y0-head_length];

% =============================================
% plot arrow  (= line + patch of a triangle)
% =============================================
h1          = plot( [x1,x2],[y1,y2],'k' );
h2          = patch( coords(1,:),coords(2,:),[0 0 0] );
    
% =============================================
% return handles
% =============================================
handles = [h1 h2];

% =============================================
% check if styling is required 
% =============================================
% if no styling, this section can be removed!
if ~isempty( varargin )
    for c = 1:floor(length(varargin)/2)
        try
            switch lower(varargin{c*2-1})

             % only patch properties    
            case 'edgecolor',   set( h2,'EdgeColor',varargin{c*2} );
            case 'facecolor',   set( h2,'FaceColor',varargin{c*2} );
            case 'facelighting',set( h2,'FaceLighting',varargin{c*2} );
            case 'edgelighting',set( h2,'EdgeLighting',varargin{c*2} );
                
            % only line properties    
            case 'color'    , set( h1,'Color',varargin{c*2} );
               
            % shared properties    
            case 'linestyle', set( handles,'LineStyle',varargin{c*2} );
            case 'linewidth', set( handles,'LineWidth',varargin{c*2} );
            case 'parent',    set( handles,'parent',varargin{c*2} );
                
            % head properties - do nothing, since handled above already
            case 'headwidth',;
            case 'headheight',;
                
            end
        catch
            fprintf( 'unrecognized property or value for: %s\n',varargin{c*2-1} );
        end
    end
end
