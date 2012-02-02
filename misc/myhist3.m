function [nn,ctrs] = myhist3(varargin)
%HIST3 Three-dimensional histogram of bivariate data.
%   HIST3(X) bins the elements of the M-by-2 matrix X into a 10-by-10 grid
%   of equally-spaced containers, and plots a histogram.  Each column of X
%   corresponds to one dimension in the bin grid.
%
%   HIST3(X,NBINS) plots a histogram using an NBINS(1)-by-NBINS(2) grid of
%   bins.  HIST3(X,'Nbins',NBINS) is equivalent to HIST3(X,NBINS).
%
%   HIST3(X,CTRS), where CTRS is a two-element cell array of numeric
%   vectors with monotonically non-decreasing values, uses a 2D grid of
%   bins centered on CTRS{1} in the first dimension and on CTRS{2} in the
%   second.  HIST3 assigns rows of X falling outside the range of that grid
%   to the bins along the outer edges of the grid, and ignores rows of X
%   containing NaNs.  HIST3(X,'Ctrs',CTRS) is equivalent to HIST3(X,CTRS).
%
%   HIST3(X,'Edges',EDGES), where EDGES is a two-element cell array
%   of numeric vectors with monotonically non-decreasing values, uses a 2D
%   grid of bins with edges at EDGES{1} in the first dimension and at
%   EDGES{2} in the second.  The (i,j)-th bin includes the value X(k,:) if
%
%      EDGES{1}(i) <= X(k,1) < EDGES{1}(i+1) and
%      EDGES{2}(j) <= X(k,2) < EDGES{2}(j+1).
%
%   Rows of X that that fall on the upper edges of the grid, EDGES{1}(end)
%   or EDGES{2}(end), are counted in the (I,j)-th or (i,J)-th bins, where
%   I and J are the lengths of EDGES{1} and EDGES{2}.  HIST3 does not count
%   rows of X falling outside the range of the grid.  Use -Inf and Inf in
%   EDGES to include all non-NaN values.
%
%   N = HIST3(X,...) returns a matrix containing the number of elements of
%   X that fall in each bin of the grid, and does not plot the histogram.
%   
%   [N,C] = HIST3(X,...) returns the positions of the bin centers in a
%   1-by-2 cell array of numeric vectors, and does not plot the histogram.
%
%   HIST3(AX,X,...) plots into AX instead of GCA.
%
%   HIST3(..., 'PARAM1',val1, 'PARAM2',val2, ...) allows you to specify
%   graphics parameter name/value pairs to fine-tune the plot.
%
%   Examples:
%
%      % Create the car data and make a histogram on a 7x7 grid of bins.
%      load carbig
%      X = [MPG,Weight];
%      hist3(X,[7 7]);
%      xlabel('MPG'); ylabel('Weight');
%
%      % Make a histogram with semi-transparent bars
%      hist3(X,[7 7],'FaceAlpha',.65);
%      xlabel('MPG'); ylabel('Weight');
%      set(gcf,'renderer','opengl');
%
%      % Specify bin centers, different in each direction.  Get back
%      % counts, but don't make the plot.
%      cnt = hist3(X, {0:10:50 2000:500:5000});
%
%   See also ACCUMARRAY, BAR, BAR3, HIST, HISTC.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2004/06/25 18:52:43 $

i = find(strcmpi(varargin,'weights'));
if ~isempty(i),
  weights = varargin{i(end)+1};
  varargin([i,i+1]) = [];
end

[cax,args,nargs] = axescheck(varargin{:});

if nargs < 1
    error('stats:hist3:TooFewInputs', 'Requires X.')
end
x = args{1};
 
% See if nbins/ctrs was given as the second argument, or only name/value
% pairs.
if nargs > 1 && ~ischar(args{2})
    binSpec = args{2};
    args = args(3:end); % strip off x and nbins/ctrs
else
    binSpec = [];
    args = args(2:end); % strip off x
end

% Process input parameter name/value pairs, assume unrecognized ones are
% graphics properties.
pnames = {'nbins','ctrs','edges'};
dflts =  { [],     [],       []};
[errid,errmsg,nbins,ctrs,edges,plotArgs] = statgetargs(pnames, dflts, args{:});
if ~isempty(errmsg)
    error(['stats:hist3:' errid], errmsg);
end

% Make sure they haven't mixed 'nbins'/'ctrs'/'edges' name/value pairs with
% the CTRS or NBINS unnamed second arg syntax, or used more than one of
% those parameter name/value pairs.
if (isempty(nbins)+isempty(ctrs)+isempty(edges)+isempty(binSpec)) < 3
    error('stats:hist3:AmbiguousBinSpec', 'Ambiguous specification of bins.');
elseif ~isempty(binSpec)
    if iscell(binSpec)  % hist3(x,ctrs)
        ctrs = binSpec;
    else                % hist3(x,nbins)
        nbins = binSpec;
    end
end

if ~isempty(nbins)
    % Use the specified number of bars in each direction, centers and edges
    % to be determined.
    histBehavior = true;
    if ~(isnumeric(nbins) && numel(nbins)==2)
        error('stats:hist3:BadNbins', ...
              'The number of bins must be specified with a 2-element numeric vector.');
    end
    autobins = true;
    
elseif ~isempty(ctrs)
    % Use the specified bin centers.
    histBehavior = true;
    if ~(iscell(ctrs) && numel(ctrs)==2 && isnumeric(ctrs{1}) && isnumeric(ctrs{2}))
        error('stats:hist3:BadCtrs', ...
              'Bin centers must be specified with a cell array containing two numeric vectors.');
    end
    ctrs = {ctrs{1}(:)' ctrs{2}(:)'};
    autobins = false;
    nbins = [length(ctrs{1}) length(ctrs{2})];
    
elseif ~isempty(edges)
    % Use the specified bin edges.
    histBehavior = false;
    if ~(iscell(edges) && numel(edges)==2 && isnumeric(edges{1}) && isnumeric(edges{2}))
        error('stats:hist3:BadEdges', ...
              'Bin edges must be specified with a cell array containing two numeric vectors.');
    end
    edges = {edges{1}(:)' edges{2}(:)'};
    autobins = false;
    % Just as with histc, there will be #edges bins
    nbins = [length(edges{1}) length(edges{2})];
    
else
    % Assume a 10x10 grid of bars, centers and edges to be determined.
    histBehavior = true;
    autobins = true;
    nbins = [10 10];
end

[nrows,ncols] = size(x);
if ncols ~= 2
    error('stats:hist3:WrongNumCols', 'X must be a matrix with two columns.');
end

% Special case for empty data (follows what HIST does).
if isempty(x)
    if autobins
       ctrs = {1:nbins(1) 1:nbins(2)};
    end
    n = zeros(nbins); % Nothing to count, return nbins(1) by nbins(2) zeros
    
else
    % Bin each observation in the x-direction, and in the y-direction.
    bin = zeros(nrows,2);
    for i = 1:2
        minx = min(x(:,i));
        maxx = max(x(:,i));
        
        % If only the number of bins was given, compute edges and centers
        % for equal-sized bins spanning the data.
        if autobins
            if isinf(minx) || isinf(maxx)
                error('stats:hist3:InfData', ...
                      'Bin centers or edges must be specified when data contain infinite values.');
            elseif minx == maxx
                minx = minx - floor(nbins(i)/2) - 0.5;
                maxx = maxx + ceil(nbins(i)/2) - 0.5;
            end
            binwidth{i} = (maxx - minx) / nbins(i);
            edges{i} = minx + binwidth{i}*(0:nbins(i));
            ctrs{i} = edges{i}(1:nbins(i)) + binwidth{i}/2;
            % Make histc mimic hist behavior:  everything < ctrs(1) gets
            % counted in first bin, everything > ctrs(end) gets counted in
            % last bin.  ctrs, edges, and binwidth do not reflect that, but
            % histcEdges does.
            histcEdges = [-Inf edges{i}(2:end-1) Inf];
            
        % If the bin centers were given, compute their edges and widths.
        elseif histBehavior
            c = ctrs{i};
            dc = diff(c);
            edges{i} = [c(1) c] + [-dc(1) dc dc(end)]/2;
            binwidth{i} = diff(edges{i});
            % Make histc mimic hist behavior:  everything < ctrs(1) gets
            % counted in first bin, everything > ctrs(end) gets counted in
            % last bin.  ctrs, edges, and binwidth do not reflect that, but
            % histcEdges does.
            histcEdges = [-Inf edges{i}(2:end-1) Inf];
            
        % If the bin edges were given, compute their widths and centers (if
        % asked for).
        else % if ~histBehavior
            e = edges{i};
            de = diff(e);
            histcEdges = e;
            % Make the display mimic bar's histc behavior: an implied bin
            % above edges(end), the same width as the last explicit one.
            % ctrs, edges, and binwidth need that explicitly, histcEdges
            % doesn't.
            edges{i} = [e e(end)+de(end)];
            binwidth{i} = [de de(end)];
            if nargout > 1
                c = zeros(size(de));
                c(1) = e(1) + de(1)/2;
                for j = 2:length(c)
                    c(j) = 2*e(j) - c(j-1);
                end
                % When edges are specified, it may not be possible to return
                % centers for which the edges are midpoints.  Warn if that's
                % the case.
                if any(c <= e(1:end-1)) || ...
                   abs(c(end) - (e(end)-de(end)/2)) > 1000*eps(de(end));
                    warning('stats:hist3:InconsistentEdges', ...
                            'Cannot compute centers that are consistent with EDGES.');
                    c = e(1:end-1) + de/2;
                end
                ctrs{i} = [c e(end)+de(end)/2];
            end
        end
        
        % Get the 1D bin numbers for this column of x.  Make sure +Inf
        % goes into the nth bin, not the (n+1)th.
        
        % bin(j,1) holds the x bin idx for data point j
        % bin(j,2) holds the y bin idx for data point j
        [dum,bin(:,i)] = histc(x(:,i),histcEdges,1);
        bin(:,i) = min(bin(:,i),nbins(i));
    end
    
    % Combine the two vectors of 1D bin counts into a grid of 2D bin
    % counts
    if exist('weights','var'),
      nzidx = all(bin>0,2);
      n = accumarray(bin(nzidx,:),weights(nzidx),nbins);
    else
      n = accumarray(bin(all(bin>0,2),:),1,nbins);
    end

end

if 0 < nargout
    nn = n;
    return
end

del = .001; % space between bars, relative to bar size

% Build x-coords for the eight corners of each bar.
xx = edges{1};
xx = [xx(1:nbins(1))+del*binwidth{1}; xx(2:nbins(1)+1)-del*binwidth{1}];
xx = [reshape(repmat(xx(:)',2,1),4,nbins(1)); repmat(NaN,1,nbins(1))];
xx = [repmat(xx(:),1,4) repmat(NaN,5*nbins(1),1)];
xx = repmat(xx,1,nbins(2));

% Build y-coords for the eight corners of each bar.
yy = edges{2};
yy = [yy(1:nbins(2))+del*binwidth{2}; yy(2:nbins(2)+1)-del*binwidth{2}];
yy = [reshape(repmat(yy(:)',2,1),4,nbins(2)); repmat(NaN,1,nbins(2))];
yy = [repmat(yy(:),1,4) repmat(NaN,5*nbins(2),1)];
yy = repmat(yy',nbins(1),1);

% Build z-coords for the eight corners of each bar.
zz = zeros(5*nbins(1), 5*nbins(2));
zz(5*(1:nbins(1))-3, 5*(1:nbins(2))-3) = n;
zz(5*(1:nbins(1))-3, 5*(1:nbins(2))-2) = n;
zz(5*(1:nbins(1))-2, 5*(1:nbins(2))-3) = n;
zz(5*(1:nbins(1))-2, 5*(1:nbins(2))-2) = n;

cax = newplot(cax);
holdState = ishold(cax);

% Plot the bars in a light steel blue.
cc = repmat(cat(3,.75,.85,.95), [size(zz) 1]);

% Plot the surface, using any specified graphics properties to override
% defaults.
h = surf(cax, xx, yy, zz, cc, 'tag','hist3', plotArgs{:});

if ~holdState
    % Set ticks for each bar if fewer than 16 and the centers/edges are
    % integers.  Otherwise, leave the default ticks alone.
    if (nbins(1)<16)
        if histBehavior && all(floor(ctrs{1})==ctrs{1})
            set(cax,'xtick',ctrs{1});
        elseif ~histBehavior && all(floor(edges{1})==edges{1})
            set(cax,'xtick',edges{1});
        end
    end
    if (nbins(2)<16)
        if histBehavior && all(floor(ctrs{2})==ctrs{2})
            set(cax,'ytick',ctrs{2});
        elseif ~histBehavior && all(floor(edges{2})==edges{2})
            set(cax,'ytick',edges{2});
        end
    end
    
    % Set the axis limits to have some space at the edges of the bars.
    dx = range(edges{1})*.05;
    dy = range(edges{2})*.05;
    set(cax,'xlim',[edges{1}(1)-dx edges{1}(end)+dx]);
    set(cax,'ylim',[edges{2}(1)-dy edges{2}(end)+dy]);
    
    view(cax,3);
    grid(cax,'on');
    set(get(cax,'parent'),'renderer','zbuffer');
end

if nargout > 0
    nn = n;
end
