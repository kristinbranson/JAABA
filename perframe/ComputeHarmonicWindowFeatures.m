function [y,feature_names,cache] = ComputeHarmonicWindowFeatures(x,varargin)

x = x(:)';
N = numel(x);
y = nan(0,N);
feature_names = {};

%% default parameters

[windows,window_offsets,...
  window_radii,min_window_radius,...
  max_window_radius,nwindow_radii] = ...
  SetDefaultWindowParameters();

% use all transformation types by default
trans_types = 'all';

% for debugging purposes
SANITY_CHECK = true;

% whether to cache results
DOCACHE = true;

% initialize empty cache
cache = InitializeCache();

% initialize feature_types already computed to empty
feature_types = {};

% harmonic features
num_harmonic = 1;

%% parse parameters

[...
  windows,...
  window_radii,window_offsets,...
  min_window_radius,max_window_radius,nwindow_radii,...
  feature_types,...
  trans_types,...
  SANITY_CHECK,...
  DOCACHE,...
  cache,...
  num_harmonic...
  ] = myparse(varargin,...
  'windows',windows,...
  'window_radii',window_radii,'window_offsets',window_offsets,...
  'min_window_radius',min_window_radius,'max_window_radius',max_window_radius,'nwindow_radii',nwindow_radii,...
  'feature_types',feature_types,...
  'trans_types',trans_types,...
  'sanitycheck',SANITY_CHECK,...
  'docache',DOCACHE,...
  'cache',cache,...
  'num_harmonic',num_harmonic); %#ok<ASGLU>

%% whether we've specified to use all trans types by default
if ischar(trans_types) && strcmpi(trans_types,'all'),
  trans_types = {'none','abs','flip'};
end

%% select default windows from various ways of specifying windows

[windows,window_radii,windowi2radiusi,nradii] = ...
  SetWindowParameters(...
  windows,window_offsets,...
  window_radii,min_window_radius,...
  max_window_radius,nwindow_radii);

%% main computation

for radiusi = 1:nradii,
  r = window_radii(radiusi);
  w = 2*r+1;
   
  for num_harmonic_curr = 1:num_harmonic,
    
    % can't split into more segments than frames
    if num_harmonic_curr >= w,
      continue;
    end
    
    fil = cos(linspace(0,pi*num_harmonic_curr,w))/w*(num_harmonic_curr+1);
    % res is nan from 1:r, N-r+1:N
    res = myconv(x,fil,nan,'corr','same');
    % boundary conditions: smallw used for frame
    for smallr = 1:r-1,
      smallw = 2*smallr+1;
      fil = cos(linspace(0,pi*num_harmonic_curr,smallw))/smallw*(num_harmonic_curr+1);
      % start
      res(smallr+1) = sum(fil.*x(1:smallw));
      % end
      res(end-smallr) = sum(fil.*x(end-smallw+1:end));
    end
    
    % all offsets for this radius
    windowis = find(windowi2radiusi == radiusi);
    for windowi = windowis',
      off = windows(windowi,2);
      res1 = padgrab(res,nan,1,1,1+off,N+off);
      
      if ismember('none',trans_types),
        
        y(end+1,:) = res1; %#ok<*AGROW>
        feature_names{end+1} = {'stat','harmonic','trans','none','radius',r,'offset',off,'num_harmonic',num_harmonic_curr};
        
        if SANITY_CHECK,
          
          res_dumb = nan(1,N);
          for n_dumb = 1:N,
            r_dumb = min([r,n_dumb+off-1,N-n_dumb-off]);
            if r_dumb < 1,
              continue;
            end
            w_dumb = 2*r_dumb+1;
            fil_dumb = cos(linspace(0,pi*num_harmonic_curr,w_dumb))/w_dumb*(num_harmonic_curr+1);
            res_dumb(n_dumb) = sum(fil_dumb.*x(n_dumb+off-r_dumb:n_dumb+off+r_dumb));
          end
          
          if any(isnan(y(end,:)) ~= isnan(res_dumb)),
            fprintf('SANITY CHECK: harmonic, trans = none, r = %d, off = %d, num_harmonic = %d, nan mismatch\n',r,off,num_harmonic_curr);
          else
            fprintf('SANITY CHECK: harmonic, trans = none, r = %d, off = %d, num_harmonic_curr = %d, max error = %f\n',r,off,num_harmonic_curr,max(abs(y(end,:)-res_dumb)));
          end
          
        end
        
      end
      
      if ismember('abs',trans_types),
        y(end+1,:) = abs(res1);
        feature_names{end+1} = {'stat','harmonic','trans','abs','radius',r,'offset',off,'num_harmonic',num_harmonic_curr};
        
        if SANITY_CHECK,
          
          res_dumb = nan(1,N);
          for n_dumb = 1:N,
            r_dumb = min([r,n_dumb+off-1,N-n_dumb-off]);
            if r_dumb < 1,
              continue;
            end
            w_dumb = 2*r_dumb+1;
            fil_dumb = cos(linspace(0,pi*num_harmonic_curr,w_dumb))/w_dumb*(num_harmonic_curr+1);
            res_dumb(n_dumb) = abs(sum(fil_dumb.*x(n_dumb+off-r_dumb:n_dumb+off+r_dumb)));
          end
          
          if any(isnan(y(end,:)) ~= isnan(res_dumb)),
            fprintf('SANITY CHECK: harmonic, trans = abs, r = %d, off = %d, num_harmonic = %d, nan mismatch\n',r,off,num_harmonic_curr);
          else
            fprintf('SANITY CHECK: harmonic, trans = abs, r = %d, off = %d, num_harmonic_curr = %d, max error = %f\n',r,off,num_harmonic_curr,max(abs(y(end,:)-res_dumb)));
          end
          
        end
      end
      
    end
  end
  
end