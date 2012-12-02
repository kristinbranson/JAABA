function [idx, C, sumD, D] = weightedkmeans(X, k, varargin)
% similar to MATLAB's kmeans function, but takes weights

if nargin < 2
    error('stats:kmeans:TooFewInputs','At least two input arguments required.');
end

if any(isnan(X(:)))
    warning('stats:kmeans:MissingDataRemoved','Removing rows of X with missing data.');
    X = X(~any(isnan(X),2),:);
end

% n points in p dimensional space
[n, p] = size(X);
Xsort = []; Xord = [];

pnames = {   'distance'  'start' 'replicates' 'maxiter' 'emptyaction' 'display' 'weights'};
dflts =  {'sqeuclidean' 'sample'          []       100        'error'  'notify' []};
[eid,errmsg,distance,start,reps,maxit,emptyact,display,weights] ...
                       = statgetargs(pnames, dflts, varargin{:});
if ~isempty(eid)
    error(sprintf('stats:kmeans:%s',eid),errmsg);
end

isweighted = ~isempty(weights);
if ~isweighted,
  weights = ones(n,1);
end

% remove data points with 0 weight
remove = isnan(weights) | weights <= 0;
X(remove,:) = [];
weights(remove) = [];
n = length(weights);
weights = weights / sum(weights);

if ischar(distance)
    distNames = {'sqeuclidean'};
    i = strmatch(lower(distance), distNames);
    if length(i) > 1
        error('stats:kmeans:AmbiguousDistance', ...
              'Ambiguous ''distance'' parameter value:  %s.', distance);
    elseif isempty(i)
        error('stats:kmeans:UnknownDistance', ...
              'Unknown ''distance'' parameter value:  %s.', distance);
    end
    distance = distNames{i};
    switch distance 
    case 'cityblock'
        [Xsort,Xord] = sort(X,1);        
    case 'cosine'
        Xnorm = sqrt(sum(X.^2, 2));
        if any(min(Xnorm) <= eps(max(Xnorm)))
            error('stats:kmeans:ZeroDataForCos', ...
                  ['Some points have small relative magnitudes, making them ', ...
                   'effectively zero.\nEither remove those points, or choose a ', ...
                   'distance other than ''cosine''.']);
        end
        X = X ./ Xnorm(:,ones(1,p));
    case 'correlation'
        X = X - repmat(mean(X,2),1,p);
        Xnorm = sqrt(sum(X.^2, 2));
        if any(min(Xnorm) <= eps(max(Xnorm)))
            error('stats:kmeans:ConstantDataForCorr', ...
                  ['Some points have small relative standard deviations, making them ', ...
                   'effectively constant.\nEither remove those points, or choose a ', ...
                   'distance other than ''correlation''.']);
        end
        X = X ./ Xnorm(:,ones(1,p));
    case 'hamming'
        if ~all(ismember(X(:),[0 1]))
            error('stats:kmeans:NonbinaryDataForHamm', ...
                  'Non-binary data cannot be clustered using Hamming distance.');
        end
    end
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
    elseif isempty(k)
        error('stats:kmeans:MissingK', ...
              'You must specify the number of clusters, K.');
    end
    start = startNames{i};
    if strcmp(start, 'uniform')
        if strcmp(distance, 'hamming')
            error('stats:kmeans:UniformStartForHamm', ...
                  'Hamming distance cannot be initialized with uniform random values.');
        end
        Xmins = min(X,[],1);
        Xmaxs = max(X,[],1);
    end
elseif isnumeric(start)
    CC = start;
    start = 'numeric';
    if isempty(k)
        k = size(CC,1);
    elseif k ~= size(CC,1);
        error('stats:kmeans:MisshapedStart', ...
              'The ''start'' matrix must have K rows.');
    elseif size(CC,2) ~= p
        error('stats:kmeans:MisshapedStart', ...
              'The ''start'' matrix must have the same number of columns as X.');
    end
    if isempty(reps)
        reps = size(CC,3);
    elseif reps ~= size(CC,3);
        error('stats:kmeans:MisshapedStart', ...
              'The third dimension of the ''start'' array must match the ''replicates'' parameter value.');
    end
    
    % Need to center explicit starting points for 'correlation'. (Re)normalization
    % for 'cosine'/'correlation' is done at each iteration.
    if isequal(distance, 'correlation')
        CC = CC - repmat(mean(CC,2),[1,p,1]);
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

if k == 1
    error('stats:kmeans:OneCluster', ...
          'The number of clusters must be greater than 1.');
elseif n < k
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
D = repmat(NaN,n,k);   % point-to-cluster distances
Del = repmat(NaN,n,k); % reassignment criterion
m = zeros(k,1);
z = zeros(k,1);

totsumDBest = Inf;
for rep = 1:reps
    switch start
    case 'uniform'
        C = unifrnd(Xmins(ones(k,1),:), Xmaxs(ones(k,1),:));
        % For 'cosine' and 'correlation', these are uniform inside a subset
        % of the unit hypersphere.  Still need to center them for
        % 'correlation'.  (Re)normalization for 'cosine'/'correlation' is
        % done at each iteration.
        if isequal(distance, 'correlation')
            C = C - repmat(mean(C,2),1,p);
        end
        if isa(X,'single')
            C = single(C);
        end
    case 'sample'
        if isweighted,
          C = X(samplepnoreplace(n,k,weights),:);
        else
          C = X(randsample(n,k),:);
        end
        if ~isfloat(C)      % X may be logical
            C = double(C);
        end
    case 'cluster'
      
        if isweighted,
          samples = samplepnoreplace(n,floor(.1*n),weights);
        else
          samples = randsample(n,floor(.1*n));
        end
        Xsubset = X(samples,:);
        weightssubset = weights(samples);
        [dum, C] = weightedkmeans(Xsubset, k, varargin{:}, 'start','sample', 'replicates',1,'weights',weightssubset);
    case 'numeric'
        C = CC(:,:,rep);
    end    
    changed = 1:k; % everything is newly assigned
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
        D(:,changed) = distfun(X, C(changed,:), distance, iter);
        
        % Compute the total sum of distances for the current configuration.
        % Can't do it first time through, there's no configuration yet.
        if iter > 0
            totsumD = sum(weights.*D((idx-1)*n + (1:n)'));
            % Test for a cycle: if objective is not decreased, back out
            % the last step and move on to the single update phase
            if prevtotsumD <= totsumD
                idx = previdx;
                [C(changed,:), z(changed), m(changed)] = gcentroids(X, idx, changed, distance, Xsort, Xord, weights);
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
            changed = 1:k;
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
        [C(changed,:), z(changed), m(changed)] = gcentroids(X, idx, changed, distance, Xsort, Xord, weights);
        iter = iter + 1;
        
        % Deal with clusters that have just lost all their members
        empties = changed(z(changed) == 0);
        if ~isempty(empties)
            switch emptyact
            case 'error'
                error('stats:kmeans:EmptyCluster', ...
                      'Empty cluster created at iteration %d.',iter);
            case 'drop'
                % Remove the empty cluster from any further processing
                D(:,empties) = NaN;
                changed = changed(z(changed) > 0);
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
                    % Find the point furthest away from its current cluster.
                    % Take that point out of its cluster and use it to create
                    % a new singleton cluster to replace the empty one.
                    [dlarge, lonely] = max(d);
                    from = idx(lonely); % taking from this cluster
                    C(i,:) = X(lonely,:);
                    m(i) = 1;
                    z(i) = weights(lonely);
                    idx(lonely) = i;
                    d(lonely) = 0;
                    
                    % Update clusters from which points are taken
                    [C(from,:), z(from), m(from)] = gcentroids(X, idx, from, distance, Xsort, Xord, weights);
                    changed = unique([changed from]);
                end
            end
        end
    end % phase one

    % Initialize some cluster information prior to phase two
    switch distance
    case 'cityblock'
        Xmid = zeros([k,p,2]);
        for i = 1:k
            if z(i) > 0
                % Separate out sorted coords for points in i'th cluster,
                % and save values above and below median, component-wise
                Xsorted = reshape(Xsort(idx(Xord)==i), m(i), p);
                nn = floor(.5*m(i));
                if mod(m(i),2) == 0
                    Xmid(i,:,1:2) = Xsorted([nn, nn+1],:)';
                elseif m(i) > 1
                    Xmid(i,:,1:2) = Xsorted([nn, nn+2],:)';
                else
                    Xmid(i,:,1:2) = Xsorted([1, 1],:)';
                end
            end
        end
    case 'hamming'
        Xsum = zeros(k,p);
        for i = 1:k
            if m(i) > 0
                % Sum coords for points in i'th cluster, component-wise
                Xsum(i,:) = sum(X(idx==i,:), 1);
            end
        end
    end
    
    %
    % Begin phase two:  single reassignments
    %
    changed = find(z' > 0);
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
                mbrs = (idx == i);
                sgn = weights.*(1 - 2*mbrs); % -w for members, w for nonmembers
                if m(i) == 1
                    sgn(mbrs) = 0; % prevent divide-by-zero for singleton mbrs
                end
                Del(:,i) = (z(i) ./ (z(i) + sgn)) .* sum((X - C(repmat(i,n,1),:)).^2, 2);
            end
        case 'cityblock'
            for i = changed
                if mod(m(i),2) == 0 % this will never catch singleton clusters
                    ldist = Xmid(repmat(i,n,1),:,1) - X;
                    rdist = X - Xmid(repmat(i,n,1),:,2);
                    mbrs = (idx == i);
                    sgn = repmat(1-2*mbrs, 1, p); % -1 for members, 1 for nonmembers
                    Del(:,i) = sum(max(0, max(sgn.*rdist, sgn.*ldist)), 2);
                else
                    Del(:,i) = sum(abs(X - C(repmat(i,n,1),:)), 2);
                end
            end
        case {'cosine','correlation'}
            % The points are normalized, centroids are not, so normalize them
            normC(changed) = sqrt(sum(C(changed,:).^2, 2));
            if any(normC < eps(class(normC))) % small relative to unit-length data points
                error('stats:kmeans:ZeroCentroid', ...
                      'Zero cluster centroid created at iteration %d.',iter);
            end
            % This can be done without a loop, but the loop saves memory allocations
            for i = changed
                XCi = X * C(i,:)';
                mbrs = (idx == i);
                sgn = 1 - 2*mbrs; % -1 for members, 1 for nonmembers
                Del(:,i) = 1 + sgn .*...
                      (m(i).*normC(i) - sqrt((m(i).*normC(i)).^2 + 2.*sgn.*m(i).*XCi + 1));
            end
        case 'hamming'
            for i = changed
                if mod(m(i),2) == 0 % this will never catch singleton clusters
                    % coords with an unequal number of 0s and 1s have a
                    % different contribution than coords with an equal
                    % number
                    unequal01 = find(2*Xsum(i,:) ~= m(i));
                    numequal01 = p - length(unequal01);
                    mbrs = (idx == i);
                    Di = abs(X(:,unequal01) - C(repmat(i,n,1),unequal01));
                    Del(:,i) = (sum(Di, 2) + mbrs*numequal01) / p;
                else
                    Del(:,i) = sum(abs(X - C(repmat(i,n,1),:)), 2) / p;
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
        
        % Update the cluster index vector, and rhe old and new cluster
        % counts and centroids
        idx(moved) = nidx;
        m(nidx) = m(nidx) + 1;
        m(oidx) = m(oidx) - 1;
        z(oidx) = z(oidx) - weights(moved);
        z(nidx) = z(nidx) + weights(moved);
        switch distance
        case 'sqeuclidean'
            C(nidx,:) = C(nidx,:) + weights(moved) * (X(moved,:) - C(nidx,:)) / z(nidx);
            C(oidx,:) = C(oidx,:) - weights(moved) * (X(moved,:) - C(oidx,:)) / z(oidx);
        case 'cityblock'
            for i = [oidx nidx]
                % Separate out sorted coords for points in each cluster.
                % New centroid is the coord median, save values above and
                % below median.  All done component-wise.
                Xsorted = reshape(Xsort(idx(Xord)==i), m(i), p);
                nn = floor(.5*m(i));
                if mod(m(i),2) == 0
                    C(i,:) = .5 * (Xsorted(nn,:) + Xsorted(nn+1,:));
                    Xmid(i,:,1:2) = Xsorted([nn, nn+1],:)';
                else
                    C(i,:) = Xsorted(nn+1,:);
                    if m(i) > 1
                        Xmid(i,:,1:2) = Xsorted([nn, nn+2],:)';
                    else
                        Xmid(i,:,1:2) = Xsorted([1, 1],:)';
                    end
                end
            end
        case {'cosine','correlation'}
            C(nidx,:) = C(nidx,:) + (X(moved,:) - C(nidx,:)) / m(nidx);
            C(oidx,:) = C(oidx,:) - (X(moved,:) - C(oidx,:)) / m(oidx);
        case 'hamming'
            % Update summed coords for points in each cluster.  New
            % centroid is the coord median.  All done component-wise.
            Xsum(nidx,:) = Xsum(nidx,:) + X(moved,:);
            Xsum(oidx,:) = Xsum(oidx,:) - X(moved,:);
            C(nidx,:) = .5*sign(2*Xsum(nidx,:) - m(nidx)) + .5;
            C(oidx,:) = .5*sign(2*Xsum(oidx,:) - m(oidx)) + .5;
        end
        changed = sort([oidx nidx]);
    end % phase two
    
    if (~converged) & (display > 0)
        warning('stats:kmeans:FailedToConverge', ...
                'Failed to converge in %d iterations.', maxit);
    end

    % Calculate cluster-wise sums of distances
    nonempties = find(m(:)'>0);
    D(:,nonempties) = distfun(X, C(nonempties,:), distance, iter);
    d = D((idx-1)*n + (1:n)');
    sumD = zeros(k,1);
    for i = 1:k
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

function [centroids, z, counts] = gcentroids(X, index, clusts, dist, Xsort, Xord, weights)
%GCENTROIDS Centroids and counts stratified by group.
[n,p] = size(X);
num = length(clusts);
centroids = repmat(NaN, [num p]);
counts = zeros(num,1);
z = zeros(num,1);
for i = 1:num
    members = find(index == clusts(i));
    counts(i) = length(members);
    z(i) = sum(weights(members));
    if z(i) > 0,
      switch dist
        case 'sqeuclidean'
          centroids(i,:) = sum(repmat(weights(members),[1,p]).*X(members,:),1) / z(i);
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
            centroids(i,:) = sum(repmat(weights(members),[1,p]).*X(members,:),1) / z(i);% unnormalized
        case 'hamming'
            % Compute a fast median for binary data, component-wise
            centroids(i,:) = .5*sign(2*sum(X(members,:).*repmat(weights(members),[1,p]), 1) - z(i)) + .5;
        end
    end
end
