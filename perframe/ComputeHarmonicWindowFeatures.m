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
%trans_types = 'all';
trans_types = uint8(15);

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

% relative params
relativeParams = [];
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
  num_harmonic,...
  relativeParams...
  ] = myparse(varargin,...
  'windows',windows,...
  'window_radii',window_radii,'window_offsets',window_offsets,...
  'min_window_radius',min_window_radius,'max_window_radius',max_window_radius,'nwindow_radii',nwindow_radii,...
  'feature_types',feature_types,...
  'trans_types',trans_types,...
  'sanitycheck',SANITY_CHECK,...
  'docache',DOCACHE,...
  'cache',cache,...
  'num_harmonic',num_harmonic,...
  'relativeParams',relativeParams); %#ok<ASGLU>

%% whether we've specified to use all trans types by default
%if ischar(trans_types) && strcmpi(trans_types,'all'),
%  trans_types = {'none','abs','flip','relative'};
%end

%% select default windows from various ways of specifying windows

[windows,window_radii,windowi2radiusi,nradii] = ...
  SetWindowParameters(...
  windows,window_offsets,...
  window_radii,min_window_radius,...
  max_window_radius,nwindow_radii);

%% main computation

%if ismember('relative',trans_types)
if bitand(8,trans_types)
  if DOCACHE && ~isempty(cache.relX)
    modX = cache.relX;
  else
    modX = convertToRelative(x,relativeParams);
    cache.relX = modX;
  end
end

for radiusi = 1:nradii,
  r = window_radii(radiusi);
  w = 2*r+1;
   
  for num_harmonic_curr = 1:num_harmonic,
    
    % can't split into more segments than frames
    if num_harmonic_curr >= w,
      continue;
    end
    
    res = HarmonicWindowCore(x,w,num_harmonic_curr);
    
    %if ismember('relative',trans_types),
    if bitand(8,trans_types),
       resRel = HarmonicWindowCore(modX,w,num_harmonic_curr);
    end
    
    % all offsets for this radius
    windowis = find(windowi2radiusi == radiusi);
    for windowi = windowis',
      off = windows(windowi,2);
      res1 = padgrab(res,nan,1,1,1+off,N+off);
      
      %if ismember('none',trans_types),
      if bitand(1,trans_types),
        y(end+1,:) = res1; %#ok<*AGROW>
        feature_names{end+1} = {'stat','harmonic','trans','none','radius',r,'offset',off,'num_harmonic',num_harmonic_curr};
      end
   
      %if ismember('abs',trans_types),
      if bitand(2,trans_types),
        y(end+1,:) = abs(res1);
        feature_names{end+1} = {'stat','harmonic','trans','abs','radius',r,'offset',off,'num_harmonic',num_harmonic_curr};
      end
      
      %if ismember('relative',trans_types),
      if bitand(8,trans_types),
        resRel1 = padgrab(resRel,nan,1,1,1+off,N+off);
        y(end+1,:) = resRel1;
        feature_names{end+1} = {'stat','harmonic','trans','relative','radius',r,'offset',off,'num_harmonic',num_harmonic_curr};
      end
      
      if SANITY_CHECK,
        extraStr = sprintf('num_harmonic = %d ',num_harmonic_curr);
        funcType = 'harmonic';   
        
        %if ismember('none',trans_types),
        if bitand(1,trans_types),
          fastY = res1; %#ok<*AGROW>
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
          checkSanity(fastY,res1,r,off,funcType,'none',extraStr);
        end
        
        %if ismember('abs',trans_types),
        if bitand(2,trans_types),
          fastY = abs(res1);
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
          checkSanity(fastY,res1,r,off,funcType,'abs',extraStr);
        end
      end
      
    end
  end
  
end
