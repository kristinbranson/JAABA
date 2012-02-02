function [idx, C, sumD, D] = eitherkmeans(X, k, varargin)
%KMEANS K-means clustering.
%   [IDX,C] = KMEANS(X, K) 
% X is a cell of length NCHOICES
% Each X{i} is a matrix of size N x P(i)
% K is a vector of length NCHOICES
% K(i) is the number of clusters using data X{i}
% This code tries to find the centers C that minimize
% \sum_{n=1}^N \min_{i=1}^NCHOICES \min_{j=1}^{K(i)} || X{i}(n,:) - C{i}(j,:) ||^2
% using a kmeans-type algorithm. 
% C is a cell of length NCHOICES. C{i} is K(i) x P(i) matrix, where
% C{i}(K(i),:) is a center. IDX(n) is the label given to data point n. 
%
%   KMEANS treats NaNs as missing data, and removes any rows of X that
%   contain NaNs. 
%
%   [IDX, C, SUMD] = KMEANS(X, K) returns the within-cluster sums of
%   point-to-centroid distances in the 1-by-K vector sumD.
%
%   [IDX, C, SUMD, D] = KMEANS(X, K) returns distances from each point
%   to every centroid in the N-by-K matrix D.
%
%   [ ... ] = KMEANS(..., 'PARAM1',val1, 'PARAM2',val2, ...) allows you to
%   specify optional parameter name/value pairs to control the iterative
%   algorithm used by KMEANS.  Parameters are:
%
%   'Distance' - Distance measure, in P-dimensional space, that KMEANS
%      should minimize with respect to.  Choices are:
%            {'sqEuclidean'} - Squared Euclidean distance
%             'cityblock'    - Sum of absolute differences, a.k.a. L1
%             'cosine'       - One minus the cosine of the included angle
%                              between points (treated as vectors)
%             'correlation'  - One minus the sample correlation between
%                              points (treated as sequences of values)
%             'Hamming'      - Percentage of bits that differ (only
%                              suitable for binary data)
%
%   'Start' - Method used to choose initial cluster centroid positions,
%      sometimes known as "seeds".  Choices are:
%                 {'sample'} - Select K observations from X at random
%                  'uniform' - Select K points uniformly at random from
%                              the range of X.  Not valid for Hamming distance.
%                  'cluster' - Perform preliminary clustering phase on
%                              random 10% subsample of X.  This preliminary
%                              phase is itself initialized using 'sample'.
%                  matrix    - A K-by-P matrix of starting locations.  In
%                              this case, you can pass in [] for K, and
%                              KMEANS infers K from the first dimension of
%                              the matrix.  You can also supply a 3D array,
%                              implying a value for 'Replicates'
%                              from the array's third dimension.
%
%   'Replicates' - Number of times to repeat the clustering, each with a
%      new set of initial centroids [ positive integer | {1}]
%
%   'Maxiter' - The maximum number of iterations [ positive integer | {100}]
%
%   'EmptyAction' - Action to take if a cluster loses all of its member
%      observations.  Choices are:
%               {'error'}    - Treat an empty cluster as an error
%                'drop'      - Remove any clusters that become empty, and
%                              set corresponding values in C and D to NaN.
%                'singleton' - Create a new cluster consisting of the one
%                              observation furthest from its centroid.
%
%   'Display' - Display level [ 'off' | {'notify'} | 'final' | 'iter' ]
%
%   Example:
%
%       X = [randn(20,2)+ones(20,2); randn(20,2)-ones(20,2)];
%       [cidx, ctrs] = kmeans(X, 2, 'dist','city', 'rep',5, 'disp','final');
%       plot(X(cidx==1,1),X(cidx==1,2),'r.', ...
%            X(cidx==2,1),X(cidx==2,2),'b.', ctrs(:,1),ctrs(:,2),'kx');
%
%   See also LINKAGE, CLUSTERDATA, SILHOUETTE.

%   KMEANS uses a two-phase iterative algorithm to minimize the sum of
%   point-to-centroid distances, summed over all K clusters.  The first
%   phase uses what the literature often describes as "batch" updates,
%   where each iteration consists of reassigning points to their nearest
%   cluster centroid, all at once, followed by recalculation of cluster
%   centroids. This phase may be thought of as providing a fast but
%   potentially only approximate solution as a starting point for the
%   second phase.  The second phase uses what the literature often
%   describes as "on-line" updates, where points are individually
%   reassigned if doing so will reduce the sum of distances, and cluster
%   centroids are recomputed after each reassignment.  Each iteration
%   during this second phase consists of one pass though all the points.
%   KMEANS can converge to a local optimum, which in this case is a
%   partition of points in which moving any single point to a different
%   cluster increases the total sum of distances.  This problem can only be
%   solved by a clever (or lucky, or exhaustive) choice of starting points.
%
% References:
%
%   [1] Seber, G.A.F., Multivariate Observations, Wiley, New York, 1984.
%   [2] Spath, H. (1985) Cluster Dissection and Analysis: Theory, FORTRAN
%       Programs, Examples, translated by J. Goldschmidt, Halsted Press,
%       New York, 226 pp.

%   Copyright 1993-2005 The MathWorks, Inc.
%   $Revision: 1.4.4.7 $  $Date: 2005/11/18 14:28:09 $

if nargin < 2
    error('stats:kmeans:TooFewInputs','At least two input arguments required.');
end

nchoices = length(X);

for choice = 1:nchoices,
  if any(isnan(X{choice}(:)))
    X{choice} = X{choice}(~any(isnan(X{choice}),2),:);
  end
end

% n points in p dimensional space
nold = size(X{1},1); p = zeros(1,nchoices);
for choice = 1:nchoices,
  [n, p(choice)] = size(X{choice});
  if nold ~= n,
    error('Number of samples should be the same for each choice');
  end
  nold = n;
end
Xsort = cell(1,nchoices); Xord = cell(1,nchoices);

pnames = {   'distance'  'start' 'replicates' 'maxiter' 'emptyaction' 'display'};
dflts =  {'sqeuclidean' 'sample'          []       100        'error'  'notify'};
[eid,errmsg,distance,start,reps,maxit,emptyact,display] ...
                       = statgetargs(pnames, dflts, varargin{:});
if ~isempty(eid)
    error(sprintf('stats:kmeans:%s',eid),errmsg);
end

if ischar(distance)
    distNames = {'sqeuclidean','cityblock','cosine','correlation','hamming'};
    i = strmatch(lower(distance), distNames);
    if length(i) > 1
        error('stats:kmeans:AmbiguousDistance', ...
              'Ambiguous ''distance'' parameter value:  %s.', distance);
    elseif isempty(i)
        error('stats:kmeans:UnknownDistance', ...
              'Unknown ''distance'' parameter value:  %s.', distance);
    end
    distance = distNames{i};
    for choice = 1:nchoices,

      switch distance 
        case 'cityblock'
          [Xsort{choice},Xord{choice}] = sort(X{choice},1);
        case 'cosine'
          Xnorm = sqrt(sum(X{choice}.^2, 2));
          if any(min(Xnorm) <= eps(max(Xnorm)))
            error('stats:kmeans:ZeroDataForCos', ...
              ['Some points have small relative magnitudes, making them ', ...
              'effectively zero.\nEither remove those points, or choose a ', ...
              'distance other than ''cosine''.']);
          end
          X{choice} = X{choice} ./ Xnorm(:,ones(1,p(choice)));
        case 'correlation'
          X{choice} = X{choice} - repmat(mean(X{choice},2),1,p(choice));
          Xnorm = sqrt(sum(X{choice}.^2, 2));
          if any(min(Xnorm) <= eps(max(Xnorm)))
            error('stats:kmeans:ConstantDataForCorr', ...
              ['Some points have small relative standard deviations, making them ', ...
              'effectively constant.\nEither remove those points, or choose a ', ...
              'distance other than ''correlation''.']);
          end
          X{choice} = X{choice} ./ Xnorm(:,ones(1,p(choice)));
        case 'hamming'
          if ~all(ismember(X{choice}(:),[0 1]))
            error('stats:kmeans:NonbinaryDataForHamm', ...
              'Non-binary data cannot be clustered using Hamming distance.');
          end
      end
      
    end % end loop over choice
else
  error('stats:kmeans:InvalidDistance', ...
    'The ''distance'' parameter value must be a string.');
end

if ischar(start)
    startNames = {'uniform','sample','cluster'};
    i = strmatch(lower(start), startNames);
    if length(i) > 1
        error('stats:kmeans:AmbiguousStart', ...
              'Ambiguous ''start'' parameter value:  %s.', start);
    elseif isempty(i)
        error('stats:kmeans:UnknownStart', ...
              'Unknown ''start'' parameter value:  %s.', start);
    elseif isempty(k) || any(k==-1)
        error('stats:kmeans:MissingK', ...
              'You must specify the number of clusters, K.');
    end
    start = startNames{i};
    if strcmp(start, 'uniform')
        if strcmp(distance, 'hamming')
            error('stats:kmeans:UniformStartForHamm', ...
                  'Hamming distance cannot be initialized with uniform random values.');
        end
        Xmins = cell(nchoices,1);
        Xmaxs = cell(nchoices,1);
        for choice = 1:nchoices,
          Xmins{choice} = min(X{choice},[],1);
          Xmaxs{choice} = max(X{choice},[],1);
        end
    end
elseif iscell(start)
  CC = start;
  start = 'numeric';
  if isempty(k)
    k = zeros(1,nchoices);
  end
  for choice = 1:nchoices,
    if k(choice) == -1
      k(choice) = size(CC{choice},1);
    elseif k(choice) ~= size(CC{choice},1);
      error('stats:kmeans:MisshapedStart', ...
        'The ''start'' matrix must have K rows.');
    elseif size(CC{choice},2) ~= p(choice)
      error('stats:kmeans:MisshapedStart', ...
        'The ''start'' matrix must have the same number of columns as X.');
    end
    if isempty(reps)
      reps = size(CC{choice},3);
    elseif reps ~= size(CC{choice},3);
      error('stats:kmeans:MisshapedStart', ...
        'The third dimension of the ''start'' array must match the ''replicates'' parameter value.');
    end
    
    % Need to center explicit starting points for 'correlation'. (Re)normalization
    % for 'cosine'/'correlation' is done at each iteration.
    if isequal(distance, 'correlation')
      CC{choice} = CC - repmat(mean(CC,2),[1,p,1]);
    end
  end
else
    error('stats:kmeans:InvalidStart', ...
          'The ''start'' parameter value must be a string or a numeric matrix or array.');
end

if ischar(emptyact)
    emptyactNames = {'error','drop','singleton'};
    i = strmatch(lower(emptyact), emptyactNames);
    if length(i) > 1
        error('stats:kmeans:AmbiguousEmptyAction', ...
              'Ambiguous ''emptyaction'' parameter value:  %s.', emptyact);
    elseif isempty(i)
        error('stats:kmeans:UnknownEmptyAction', ...
              'Unknown ''emptyaction'' parameter value:  %s.', emptyact);
    end
    emptyact = emptyactNames{i};
else
    error('stats:kmeans:InvalidEmptyAction', ...
          'The ''emptyaction'' parameter value must be a string.');
end

if ischar(display)
    i = strmatch(lower(display), strvcat('off','notify','final','iter'));
    if length(i) > 1
        error('stats:kmeans:AmbiguousDisplay', ...
              'Ambiguous ''display'' parameter value:  %s.', display);
    elseif isempty(i)
        error('stats:kmeans:UnknownDisplay', ...
              'Unknown ''display'' parameter value:  %s.', display);
    end
    display = i-1;
else
    error('stats:kmeans:InvalidDisplay', ...
          'The ''display'' parameter value must be a string.');
end

% for going to/from global label from/to local label
ktotal = sum(k);
kmax = cumsum(k);
koffset = [0,kmax(1:end-1)];
kmin = koffset+1;
% map from/to global k to/from local k
kglobtoloc = zeros(1,ktotal);
kloctoglob = cell(1,nchoices);
kglobtochoice = zeros(1,ktotal);
for choice = 1:nchoices,
  kglobtoloc(koffset(choice)+(1:k(choice))) = 1:k(choice);
  kglobtochoice(koffset(choice)+(1:k(choice))) = choice;
  kloctoglob{choice} = kmin(choice):kmax(choice);
end  

if ktotal == 1
    error('stats:kmeans:OneCluster', ...
          'The number of clusters must be greater than 1.');
elseif n < ktotal
    error('stats:kmeans:TooManyClusters', ...
          'X must have more rows than the number of clusters.');
end

% Assume one replicate
if isempty(reps)
    reps = 1;
end

%
% Done with input argument processing, begin clustering
%

dispfmt = '%6d\t%6d\t%8d\t%12g';
D = nan(n,ktotal);   % point-to-cluster distances
Del = nan(n,ktotal); % reassignment criterion
m = zeros(ktotal,1);

totsumDBest = Inf;
for rep = 1:reps
  switch start
    case 'uniform'
      for choice = 1:nchoices,
        C{choice} = unifrnd(Xmins{choice}(ones(k(choice),1),:), Xmaxs{choice}(ones(k(choice),1),:));
        % For 'cosine' and 'correlation', these are uniform inside a subset
        % of the unit hypersphere.  Still need to center them for
        % 'correlation'.  (Re)normalization for 'cosine'/'correlation' is
        % done at each iteration.
        if isequal(distance, 'correlation')
            C{choice} = C{choice} - repmat(mean(C{choice},2),1,p(choice));
        end
        if isa(X{choice},'single')
            C{choice} = single(C{choice});
        end
      end
    case 'sample'
      for choice = 1:nchoices,
        if k(choice) <= 0, continue; end;
        C{choice} = X{choice}(randsample(n,k(choice)),:);
        if ~isfloat(C{choice})      % X may be logical
            C = double(C{choice});
        end
      end
    case 'cluster'
      sampleidx = randsample(n,floor(.1*n));
      Xsubset = cell(1,nchoices);
      for choice = 1:nchoices,
        Xsubset{choice} = X{choice}(sampleidx,:);
      end
      [dum, C] = eitherkmeans(Xsubset, k, varargin{:}, 'start','sample', 'replicates',1);
    case 'numeric'
      for choice = 1:nchoices,
        C{choice} = CC{choice}(:,:,rep);
      end
  end
    
  changed = 1:ktotal;
  idx = zeros(n,1);
  totsumD = Inf;
    
  if display > 2 % 'iter'
    disp(sprintf('  iter\t phase\t     num\t         sum'));
  end
    
  %
  % Begin phase one:  batch reassignments
  %
  
  converged = false;
  iter = 0;
  while true
    % Compute the distance from every point to each cluster centroid
    for iglob = changed,
      choice = kglobtochoice(iglob);
      iloc = kglobtoloc(iglob);
      D(:,iglob) = distfun(X{choice}, C{choice}(iloc,:), distance, iter);
    end
    
    % Compute the total sum of distances for the current configuration.
    % Can't do it first time through, there's no configuration yet.
    if iter > 0
      totsumD = sum(D((idx-1)*n + (1:n)'));
      % Test for a cycle: if objective is not decreased, back out
      % the last step and move on to the single update phase
      if prevtotsumD <= totsumD
        idx = previdx;
        for iglob = changed,
          choice = kglobtochoice(iglob);
          iloc = kglobtoloc(iglob);
          [C{choice}(iloc,:), m(iglob)] = gcentroids(X{choice}, idx, iglob, distance, Xsort{choice}, Xord{choice});
        end
        iter = iter - 1;
        break;
        
      end
      if display > 2 % 'iter'
        disp(sprintf(dispfmt,iter,1,length(moved),totsumD));
      end
      if iter >= maxit, break; end
    end
    
    % Determine closest cluster for each point and reassign points to clusters
    previdx = idx;
    prevtotsumD = totsumD;
    [d, nidx] = min(D, [], 2);
    
    if iter == 0
      % Every point moved, every cluster will need an update
      moved = 1:n;
      idx = nidx;
      changed = 1:ktotal;
    else
      % Determine which points moved
      moved = find(nidx ~= previdx);
      if length(moved) > 0
        % Resolve ties in favor of not moving
        moved = moved(D((previdx(moved)-1)*n + moved) > d(moved));
      end
      if length(moved) == 0
        break;
      end
      idx(moved) = nidx(moved);
      
      % Find clusters that gained or lost members
      changed = unique([idx(moved); previdx(moved)])';
    end
    
    % Calculate the new cluster centroids and counts.
    for iglob = changed,
      choice = kglobtochoice(iglob);
      iloc = kglobtoloc(iglob);
      [C{choice}(iloc,:), m(iglob)] = gcentroids(X{choice}, idx, iglob, distance, Xsort{choice}, Xord{choice});
    end
    iter = iter + 1;
    
    % Deal with clusters that have just lost all their members
    empties = changed(m(changed) == 0);
    if ~isempty(empties)
      switch emptyact
        case 'error'
          error('stats:kmeans:EmptyCluster', ...
            'Empty cluster created at iteration %d.',iter);
        case 'drop'
          % Remove the empty cluster from any further processing
          D(:,empties) = NaN;
          changed = changed(m(changed) > 0);
          if display > 0
            warning('stats:kmeans:EmptyCluster', ...
              'Empty cluster created at iteration %d.',iter);
          end
        case 'singleton'
          if display > 0
            warning('stats:kmeans:EmptyCluster', ...
              'Empty cluster created at iteration %d.',iter);
          end
          
          for i = empties
            % Find the point furthest away from its current cluster with
            % choice of i
            % Take that point out of its cluster and use it to create
            % a new singleton cluster to replace the empty one.
            
            % choice of empty cluster
            choice = kglobtochoice(i);
            iloc = kglobtoloc(i);
            % labels with this choice
            idxchoice = kmin(choice):kmax(choice);
            
            % points with any of these labels
            idxpossible = ismember(idx,idxchoice);
            
            % point with one of these labels and largest error
            [dlarge, lonely] = max(d(idxpossible));
            idxpossible = find(idxpossible);
            lonely = idxpossible(lonely);
            
            fromglob = idx(lonely); % taking from this cluster
            fromloc = kglobtoloc(fromglob);
            
            C{choice}(iloc,:) = X{choice}(lonely,:);
            m(i) = 1;
            idx(lonely) = i;
            d(lonely) = 0;
            
            % Update clusters from which points are taken
            [C{choice}(fromloc,:), m(fromglob)] = ...
              gcentroids(X{choice}, idx, fromglob, distance, Xsort{choice}, Xord{choice});
            changed = unique([changed fromglob]);
          end
      end
    end
  end % phase one

  % Initialize some cluster information prior to phase two
  switch distance
    case 'cityblock'
      for choice = 1:nchoices,
        Xmid{choice} = zeros([k(choice),p(choice),2]);
      end
      for i = 1:ktotal
        if m(i) > 0
          iloc = kglobtoloc(i);
          choice = kglobtochoice(i);
          % Separate out sorted coords for points in i'th cluster,
          % and save values above and below median, component-wise
          Xsorted{choice} = reshape(Xsort{choice}(idx(Xord{choice})==i), m(i), p(choice));
          nn = floor(.5*m(i));
          if mod(m(i),2) == 0
            Xmid{choice}(iloc,:,1:2) = Xsorted{choice}([nn, nn+1],:)';
          elseif m(i) > 1
            Xmid{choice}(iloc,:,1:2) = Xsorted{choice}([nn, nn+2],:)';
          else
            Xmid{choice}(iloc,:,1:2) = Xsorted{choice}([1, 1],:)';
          end
        end
      end
    case 'hamming'
      for choice = 1:nchoices,
        Xsum{choice} = zeros(k(choice),p(choice));
      end
      for i = 1:ktotal
        iloc = kglobtoloc(i);
        choice = kglobtochoice(i);
        if m(i) > 0
          % Sum coords for points in i'th cluster, component-wise
          Xsum{choice}(iloc,:) = sum(X{choice}(idx==i,:), 1);
        end
      end
  end
    
  %
  % Begin phase two:  single reassignments
  %
  changed = find(m' > 0);
  lastmoved = 0;
  nummoved = 0;
  iter1 = iter;
  while iter < maxit
    % Calculate distances to each cluster from each point, and the
    % potential change in total sum of errors for adding or removing
    % each point from each cluster.  Clusters that have not changed
    % membership need not be updated.
    %
    % Singleton clusters are a special case for the sum of dists
    % calculation.  Removing their only point is never best, so the
    % reassignment criterion had better guarantee that a singleton
    % point will stay in its own cluster.  Happily, we get
    % Del(i,idx(i)) == 0 automatically for them.
    switch distance
      case 'sqeuclidean'
        for i = changed
          
          iloc = kglobtoloc(i);
          choice = kglobtochoice(i);
          
          mbrs = (idx == i);
          sgn = 1 - 2*mbrs; % -1 for members, 1 for nonmembers
          if m(i) == 1
            sgn(mbrs) = 0; % prevent divide-by-zero for singleton mbrs
          end
          Del(:,i) = (m(i) ./ (m(i) + sgn)) .* sum((X{choice} - C{choice}(repmat(iloc,n,1),:)).^2, 2);
        end
      case 'cityblock'
        for i = changed
          
          iloc = kglobtoloc(i);
          choice = kglobtochoice(i);
          
          if mod(m(i),2) == 0 % this will never catch singleton clusters
            ldist = Xmid{choice}(repmat(iloc,n,1),:,1) - X;
            rdist = X{choice} - Xmid{choice}(repmat(iloc,n,1),:,2);
            mbrs = (idx == i);
            sgn = repmat(1-2*mbrs, 1, p(choice)); % -1 for members, 1 for nonmembers
            Del(:,i) = sum(max(0, max(sgn.*rdist, sgn.*ldist)), 2);
          else
            Del(:,i) = sum(abs(X{choice} - C{choice}(repmat(iloc,n,1),:)), 2);
          end
        end
      case {'cosine','correlation'}
        % The points are normalized, centroids are not, so normalize them
        for i = changed,
          iloc = kglobtoloc(i);
          choice = kglobtochoice(i);
          normC{choice}(iloc) = sqrt(   sum( C{choice}(iloc,:).^2,2 )  );
          if normC{choice}(iloc) < eps(class(normC{choice}(iloc))) % small relative to unit-length data points
            error('stats:kmeans:ZeroCentroid', ...
              'Zero cluster centroid created at iteration %d.',iter);
          end
        end
        % This can be done without a loop, but the loop saves memory allocations
        for i = changed
          iloc = kglobtoloc(i);
          choice = kglobtochoice(i);
          XCi = X{choice} * C{choice}(iloc,:)';
          mbrs = (idx == i);
          sgn = 1 - 2*mbrs; % -1 for members, 1 for nonmembers
          Del(:,i) = 1 + sgn .*...
            (m(i).*normC{choice}(iloc) - sqrt((m(i).*normC{choice}(iloc)).^2 + 2.*sgn.*m(i).*XCi + 1));
        end
      case 'hamming'
        for i = changed
          iloc = kglobtoloc(i);
          choice = kglobtochoice(i);
          if mod(m(i),2) == 0 % this will never catch singleton clusters
            % coords with an unequal number of 0s and 1s have a
            % different contribution than coords with an equal
            % number
            unequal01 = find(2*Xsum{choice}(iloc,:) ~= m(i));
            numequal01 = p(choice) - length(unequal01);
            mbrs = (idx == i);
            Di = abs(X{choice}(:,unequal01) - C{choice}(repmat(iloc,n,1),unequal01));
            Del(:,i) = (sum(Di, 2) + mbrs*numequal01) / p(choice);
          else
            Del(:,i) = sum(abs(X{choice} - C{choice}(repmat(iloc,n,1),:)), 2) / p(choice);
          end
        end
    end

    % Determine best possible move, if any, for each point.  Next we
    % will pick one from those that actually did move.
    previdx = idx;
    prevtotsumD = totsumD;
    [minDel, nidx] = min(Del, [], 2);
    moved = find(previdx ~= nidx);
    if length(moved) > 0
      % Resolve ties in favor of not moving
      moved = moved(Del((previdx(moved)-1)*n + moved) > minDel(moved));
    end
    if length(moved) == 0
      % Count an iteration if phase 2 did nothing at all, or if we're
      % in the middle of a pass through all the points
      if (iter - iter1) == 0 | nummoved > 0
        iter = iter + 1;
        if display > 2 % 'iter'
          disp(sprintf(dispfmt,iter,2,nummoved,totsumD));
        end
      end
      converged = true;
      break;
    end
    
    % Pick the next move in cyclic order
    moved = mod(min(mod(moved - lastmoved - 1, n) + lastmoved), n) + 1;
    
    % If we've gone once through all the points, that's an iteration
    if moved <= lastmoved
      iter = iter + 1;
      if display > 2 % 'iter'
        disp(sprintf(dispfmt,iter,2,nummoved,totsumD));
      end
      if iter >= maxit, break; end
      nummoved = 0;
    end
    nummoved = nummoved + 1;
    lastmoved = moved;
    
    oidx = idx(moved);
    nidx = nidx(moved);
    totsumD = totsumD + Del(moved,nidx) - Del(moved,oidx);
    oidxloc = kglobtoloc(oidx);
    nidxloc = kglobtoloc(nidx);
    ochoice = kglobtochoice(oidx);
    nchoice = kglobtochoice(nidx);
    
    % Update the cluster index vector, and rhe old and new cluster
    % counts and centroids
    idx(moved) = nidx;
    m(nidx) = m(nidx) + 1;
    m(oidx) = m(oidx) - 1;
    switch distance
      case 'sqeuclidean'
        C{nchoice}(nidxloc,:) = C{nchoice}(nidxloc,:) + (X{nchoice}(moved,:) - C{nchoice}(nidxloc,:)) / m(nidx);
        C{ochoice}(oidxloc,:) = C{ochoice}(oidxloc,:) - (X{ochoice}(moved,:) - C{ochoice}(oidxloc,:)) / m(oidx);
      case 'cityblock'
        for i = [oidx nidx]
          choice = kglobtochoice(i);
          iloc = kglobtoloc(i);
          % Separate out sorted coords for points in each cluster.
          % New centroid is the coord median, save values above and
          % below median.  All done component-wise.
          Xsorted = reshape(Xsort{choice}(idx(Xord{choice})==i), m(i), p);
          nn = floor(.5*m(i));
          if mod(m(i),2) == 0
            C{choice}(iloc,:) = .5 * (Xsorted(nn,:) + Xsorted(nn+1,:));
            Xmid{choice}(iloc,:,1:2) = Xsorted([nn, nn+1],:)';
          else
            C{choice}(iloc,:) = Xsorted(nn+1,:);
            if m(i) > 1
              Xmid{choice}(iloc,:,1:2) = Xsorted([nn, nn+2],:)';
            else
              Xmid{choice}(iloc,:,1:2) = Xsorted([1, 1],:)';
            end
          end
        end
      case {'cosine','correlation'}
        C{nchoice}(nidxloc,:) = C{nchoice}(nidxloc,:) + (X{nchoice}(moved,:) - C{nchoice}(nidxloc,:)) / m(nidx);
        C{ochoice}(oidxloc,:) = C{ochoice}(oidxloc,:) - (X{ochoice}(moved,:) - C{ochoice}(oidxloc,:)) / m(oidx);
      case 'hamming'
        % Update summed coords for points in each cluster.  New
        % centroid is the coord median.  All done component-wise.
        Xsum{nchoice}(nidxloc,:) = Xsum{nchoice}(nidxloc,:) + X{nchoice}(moved,:);
        Xsum{ochoice}(oidxloc,:) = Xsum{ochoice}(oidxloc,:) - X{ochoice}(moved,:);
        C{nchoice}(nidxloc,:) = .5*sign(2*Xsum{nchoice}(nidxloc,:) - m(nidx)) + .5;
        C{ochoice}(oidxloc,:) = .5*sign(2*Xsum{ochoice}(oidxloc,:) - m(oidx)) + .5;
    end
    changed = sort([oidx nidx]);
  end % phase two
  
  if (~converged) & (display > 0)
    warning('stats:kmeans:FailedToConverge', ...
      'Failed to converge in %d iterations.', maxit);
  end
  
  % Calculate cluster-wise sums of distances
  nonempties = find(m(:)'>0);

  for iglob = nonempties,
    choice = kglobtochoice(iglob);
    iloc = kglobtoloc(iglob);
    D(:,iglob) = distfun(X{choice}, C{choice}(iloc,:), distance, iter);
  end

  d = D((idx-1)*n + (1:n)');
  sumD = zeros(ktotal,1);
  for i = 1:ktotal
    sumD(i) = sum(d(idx == i));
  end
  if display > 1 % 'final' or 'iter'
    disp(sprintf('%d iterations, total sum of distances = %g',iter,totsumD));
  end
  
  % Save the best solution so far
  if totsumD < totsumDBest
    totsumDBest = totsumD;
    idxBest = idx;
    Cbest = C;
    sumDBest = sumD;
    if nargout > 3
      Dbest = D;
    end
  end
end

% Return the best solution
idx = idxBest;
C = Cbest;
sumD = sumDBest;
if nargout > 3
  D = Dbest;
end


%------------------------------------------------------------------

function D = distfun(X, C, dist, iter)
%DISTFUN Calculate point to cluster centroid distances.
[n,p] = size(X);
D = zeros(n,size(C,1));
nclusts = size(C,1);

switch dist
case 'sqeuclidean'
    for i = 1:nclusts
        D(:,i) = sum((X - C(repmat(i,n,1),:)).^2, 2);
    end
case 'cityblock'
    for i = 1:nclusts
        D(:,i) = sum(abs(X - C(repmat(i,n,1),:)), 2);
    end
case {'cosine','correlation'}
    % The points are normalized, centroids are not, so normalize them
    normC = sqrt(sum(C.^2, 2));
    if any(normC < eps(class(normC))) % small relative to unit-length data points
        error('stats:kmeans:ZeroCentroid', ...
              'Zero cluster centroid created at iteration %d.',iter);
    end
    % This can be done without a loop, but the loop saves memory allocations
    for i = 1:nclusts
        D(:,i) = 1 - (X * C(i,:)') ./ normC(i);
    end
case 'hamming'
    for i = 1:nclusts
        D(:,i) = sum(abs(X - C(repmat(i,n,1),:)), 2) / p;
    end
end


%------------------------------------------------------------------

function [centroids, counts] = gcentroids(X, index, clusts, dist, Xsort, Xord)
%GCENTROIDS Centroids and counts stratified by group.
[n,p] = size(X);
num = length(clusts);
centroids = repmat(NaN, [num p]);
counts = zeros(num,1);
for i = 1:num
    members = find(index == clusts(i));
    if length(members) > 0
        counts(i) = length(members);
        switch dist
        case 'sqeuclidean'
            centroids(i,:) = sum(X(members,:),1) / counts(i);
        case 'cityblock'
            % Separate out sorted coords for points in i'th cluster,
            % and use to compute a fast median, component-wise
            Xsorted = reshape(Xsort(index(Xord)==clusts(i)), counts(i), p);
            nn = floor(.5*counts(i));
            if mod(counts(i),2) == 0
                centroids(i,:) = .5 * (Xsorted(nn,:) + Xsorted(nn+1,:));
            else
                centroids(i,:) = Xsorted(nn+1,:);
            end
        case {'cosine','correlation'}
            centroids(i,:) = sum(X(members,:),1) / counts(i); % unnormalized
        case 'hamming'
            % Compute a fast median for binary data, component-wise
            centroids(i,:) = .5*sign(2*sum(X(members,:), 1) - counts(i)) + .5;
        end
    end
end
