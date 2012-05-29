classdef Trx < handle

  % Trx class
  %
  % The main data consists of any number of experiments (each corresponding
  % to one video of multiple flies). 
  % ROOTREADDIR: There is a root directory ROOTREADDIR that is the root
  % directory for all data we only need read access to -- currently the
  % video and metadata file (someday, when Ctrax is part of the pipeline,
  % the trajectories output by Ctrax will be included here). 
  % ROOTWRITEDIR: There is a root directory ROOTWRITEDIR that contains all
  % data for which we need write access; data that can be created by the
  % FlyBowlAnalysis code. This is currently the registered trajectories, sex 
  % classification, per-frame derived behavior statistics. 
  % EXPDIR_BASES: Within each of these root directories, there will be a
  % subdirectory of the same name for each experiment (name stored in
  % EXPDIR_BASES). 
  % CTRAXFILESTR: For each experiment subdirectory of ROOTWRITEDIR
  % (currently), there will be a mat file containing the trajectories
  % output by Ctrax. The name of this file within the experiment
  % subdirectory is defined by CTRAXFILESTR. 
  % REGISTER: The trajectories can be registered (pixels converted to
  % milimeters, coordinate system translated and rotated so that they align
  % between experiments) with the function REGISTER. 
  % TRXFILESTR: The registered trajectories are stored in the file
  % TRXFILESTR within the experiment directory in ROOTWRITEDIR. 
  % DETECTREGISTRATIONPARAMS: The parameters for REGISTER are stored in the
  % struct DETECTREGISTRATIONPARAMS. 
  % PERFRAMEDIR: The derived, per-frame time series for a given experiment
  % are stored within the subdirectory PERFRAMEDIR of each experiment
  % subdirectory in ROOTWRITEDIR. Within PERFRAMEDIR, there is a mat file
  % for each type of per-frame derived measurement, e.g. speed, change in
  % orientation. These mat files each contain the time series for all flies
  % in the given experiment, and are created or loaded on demand. 
  % PERFRAMEHISTORY: For each per-frame derived measure, the last access
  % time for any experiment is stored in the NFNS x 2 PERFRAMEHISTORY. 
  % MAXDATACACHED: When the total number of doubles stored in data exceeds
  % MAXDATACACHED, the least recently used per-frame derived measurement is
  % cleared from the memory cache. 

  
  properties (SetAccess = protected, GetAccess = public)

    % number of flies
    nflies = 0;
    
    % number of experiment directories
    nexpdirs = 0;
    
    % number of flies per movie
    nfliespermovie = [];
    
    %% parameters of data locations
    
    trxfilestr = 'trx.mat';
    perframedir = 'perframe';
    moviefilestr = 'movie.ufmf';

    %% data locations
        
    % there may be two directories for the data
    rootwritedir = 0;
    
    % names of the experiment directories
    expdirs = {};
    outexpdirs = {};

    % names of mat files containing registered, processed trajectories,
    % with sex classified
    trxfiles = {};

    % movie locations
    movienames = {};
    
    %% landmark parameters

    landmark_params = [];
    
    %% per-frame parameters
    perframe_params = [];
    % default field of view for computing angle subtended
    %fov = pi;
    
    % smoothing orientation
    %thetafil = [1 4 6 4 1]/16;
    
    % units of per-frame properties
    units = struct;

    
    %% sex classifier parameters
    %sexclassifier_params = [];
        
    %% data caching
    
    % history of per-frame properties loaded & their last access time
    perframehistory = cell(0,2);
    
    % maximum number of doubles in data
    maxdatacached = 2^27;
    
    % number of doubles currently stored
    ndatacached = 0;
    
    % per experiment
    ndatacachedperexp = [];
    
    % data cached for each experiment
    nfnscached = [];
    datacached = {};
    fnscached = {};
    
    %% video info
    
    % number of frames in the video
    movie_nframes = [];
    
    % image sizes
    nrs = [];
    ncs = [];
    ncolors = [];
    
    % movie size in mm
    width_mms = [];
    height_mms = [];
    
    pxpermm = [];
    fps = [];
    
    %% trajectory frames
    
    firstframes = [];
    endframes = [];
    nframes = [];

    %% indexing stuff
    
    exp2flies = {};
    fly2exp = [];
    
  end
  
  properties (Hidden = true)
    
  end
  
  properties (Constant)
    
  end
  
  methods
    
    function obj = Trx(varargin)

      % all arguments are parameters
      for i = 1:2:nargin-1,
        obj.(varargin{i}) = varargin{i+1};
      end
      
      %obj.ReadAllParams();
      
    end
        
    % deconstructor
    function delete(obj)

      obj.nflies = 0;
      obj.expdirs = {};
      obj.perframehistory = cell(0,2);
      obj.ndatacached = 0;
      obj.datacached = {};
      obj.nfnscached = [];
      
    end
    
    function flyidx = getFlyIdx(obj,n,fly)
      flyidx = obj.exp2flies{n}(fly);
    end
    
    function [n,fly] = getExpFly(obj,flyidx)
      n = obj.fly2exp(flyidx);
      fly = nan(size(flyidx));
      for i = 1:numel(flyidx),
        flycurr = find(obj.exp2flies{n} == flyidx(i),1);
        if isempty(flycurr), 
          error('Sanity check: flyidx %d mapped to exp %d, but this is not in exp2flies{%d}',flyidx(i),n,n);
        end
        fly(i) = flycurr;
      end
    end
    
    function SetTrxFileStr(obj,trxfilestr)

      obj.trxfilestr = trxfilestr;
      
    end
    
    function SetPerFrameDir(obj,perframedir)

      obj.perframedir = perframedir;
      
    end
    
    function SetMoveFileStr(obj,moviefilestr)
      
      obj.moviefilestr = moviefilestr;
      
    end

    
    function SetLandmarkParams(obj,landmark_params)
      
      obj.landmark_params = landmark_params;
      % TODO: invalidate features computed with the old parameters
      
    end
    
    function SetPerFrameParams(obj,perframe_params)
      
      obj.perframe_params = perframe_params;
      % TODO: invalidate features computed with the old parameters

    end
    
%     function SetSexClassifierParams(obj,sexclassifier_params)
%       
%       obj.sexclassifier_params = sexclassifier_params;
%       % TODO: invalidate features computed with the old parameters
% 
%     end

%     function n = numel(obj)
%       
%       n = obj.nflies;
%       
%     end
% 
%     function varargout = size(obj)
% 
%       n = obj.nflies;
%       if nargout <= 1,
%         varargout = {n};
%       else
%         varargout = num2cell([1,n,ones(1,nargout-2)]);
%       end
%       
%     end

    function [data,units] = ComputePerFrameData(obj,fn,n)
      
      m_magveldiff = regexp(fn,'^magveldiff_(.*)$','tokens','once');
      m_veltoward = regexp(fn,'^veltoward_(.*)$','tokens','once');
      m_absthetadiff = regexp(fn,'^absthetadiff_(.*)$','tokens','once');
      m_absphidiff = regexp(fn,'^absphidiff_(.*)$','tokens','once');
      m_anglefrom1to2 = regexp(fn,'^anglefrom1to2_(.*)$','tokens','once');
      m_absanglefrom1to2 = regexp(fn,'^absanglefrom1to2_(.*)$','tokens','once');
      m_dnose2ellanglerange = regexp(fn,'^dnose2ell_angle_(\w+)to(\w+)$','tokens','once');
      m_closestfly_nose2ellanglerange = regexp(fn,'^closestfly_nose2ell_angle_(\w+)to(\w+)$','tokens','once');
      if ~isempty(m_magveldiff),
        [data,units] = compute_magveldiff(obj,n,m_magveldiff{1});
      elseif ~isempty(m_veltoward),
        [data,units] = compute_veltoward(obj,n,m_veltoward{1});
      elseif ~isempty(m_absthetadiff),
        [data,units] = compute_absthetadiff(obj,n,m_absthetadiff{1});
      elseif ~isempty(m_absphidiff),
        [data,units] = compute_absphidiff(obj,n,m_absphidiff{1});
      elseif ~isempty(m_anglefrom1to2),
        [data,units] = compute_anglefrom1to2(obj,n,m_anglefrom1to2{1});
      elseif ~isempty(m_absanglefrom1to2),
        [data,units] = compute_absanglefrom1to2(obj,n,m_absanglefrom1to2{1});
      elseif ~isempty(m_dnose2ellanglerange),
        m_dnose2ellanglerange = strrep(m_dnose2ellanglerange,'min','-');
        v = str2double(m_dnose2ellanglerange);
        if numel(v) ~= 2 || any(isnan(v)),
          error('anglerange must be something like min30to30');
        end
        v = v*pi/180;
        [data,units] = compute_dnose2ell_anglerange(obj,n,v);
      elseif ~isempty(m_closestfly_nose2ellanglerange),
        [data,units] = compute_closestfly_nose2ell_anglerange(obj,n,str2double(m_closestfly_nose2ellanglerange));
      else
        funname = sprintf('compute_%s',fn);
        [data,units] = feval(funname,obj,n);
      end
      filename = obj.GetPerFrameFile(fn,n);
      try
        save(filename,'data','units');
      catch ME
        warning('Could not save %s data to %s: %s',fn,filename,getReport(ME));
      end
      
    end

    %
    
    % function declarations

    % read parameters
%     ReadAllParams(obj,varargin)
%     ReadLandmarkParams(obj)
%     ReadPerFrameParams(obj)
%     ReadSexClassifierParams(obj)
        
    % AddExpDir(obj,expdir)
    % add a new experiment directory
    AddExpDir(obj,expdir,varargin)
    
    % RemoveExpDir(obj,expdir)
    % remove an experiment directory
    RemoveExpDir(obj,expdir)
        
    % StoreTrajectories(obj,n,traj)
    % store trajectories traj for experiment n
    StoreTrajectories(obj,n,traj,varargin)
    
    FreeDataCache(obj,ndataadd)
    
    [varargout] = subsref(obj,s)
    
    [varargout] = GetPerFrameData(obj,fn,varargin)
    
    SetPerFrameData(obj,fn,x,varargin)
    
    res = PerFrameDataExists(obj,fn,n)
    
    UpdatePerFrameAccessTime(obj,fn)
    
    x = LoadPerFrameData(obj,fn,n)
    
    labelidx = LoadLabelsFromFile(obj,labelfilestr,varargin)
    
    [scores,labelidx] = LoadScoresFromFile(obj,labelfilestr,varargin)
    
    
    CleanPerFrameData(obj,fns,varargin)
    
    ClearDataCache(obj,fn,ns)
    
    varargout = drawfly(obj,fly,t,varargin)
    
  end
  
  methods(Static)
    
    fns = PerFrameFieldNames()
    
    fns = TrajectoryFieldNames()
    
  end
  
end