classdef ProjectManager < handle
  
  properties (Access = public)
    projfile = '';
    curproj = [];
    projparams = [];
  end
  
  methods (Access = public)
    
    function obj = ProjectManager(projfile)
      
      if nargin < 1
        projfile = fullfile('params','BehaviorList.xml');
      end
      obj.projfile = projfile;
      
      if ~exist(obj.projfile,'file'),
        projs = {};
      else
        in =  ReadXMLParams(obj.projfile);
        projs = fieldnames(in);
      end
      
      for ndx = 1:numel(projs),
        obj.projparams(ndx).name = projs{ndx};
        obj.projparams(ndx).configfile = in.(projs{ndx}).configfile;
        obj.projparams(ndx).save = false;
        
        if ~exist(obj.projparams(ndx).configfile,'file'),
          uiwait(warndlg(sprintf('Config file %s does not exist for project %s\n Removing the project',...
            obj.projparams(ndx).configfile,obj.projparams(ndx).name)));
          obj.projparams(ndx) = [];
        end
        
        curparams = ReadXMLParams(obj.projparams(ndx).configfile);
        if isfield(curparams,'featureparamlist'),
          obj.projparams(ndx).pfList = fieldnames(curparams.featureparamlist);
          curparams = rmfield(curparams,'featureparamlist');
        else
          obj.projparams(ndx).pfList = [];
        end
        obj.projparams(ndx).config = curparams;
      end
      
      if numel(obj.projparams)>0,
        obj.curproj = 1;
      end
      
    end
    
    function projnum = FindProjFromConfigfile(obj,configfilename)
      projnum = [];
      for fndx = 1:numel(obj.projparams)
        if strcmp(configfilename,  ...
            obj.projparams(ndx).configFile),
          projnum = fndx;
          break;
        end
      end

    end
    
    function doesExist = checkExist(obj,projname)
      doesExist = any(strcmp(projname,{obj.projparams(:).name}));
    end
    
    function defaultConfig = DefaultConfig(obj,name)
      % TODO: 
      if nargin<3
        name = 'default';
      end
      defaultConfig.targets = struct('type','fly');
      defaultConfig.behaviors = struct('names',name,...
         'labelcolors',[0.7,0,0,0,0,0.7],...
         'unknowncolor',[0,0,0]);
      defaultConfig.file = struct('moviefilename','movie.ufmf',...
        'trxfilename','registered_trx.mat',...
        'labelfilename',sprintf('labeled%s.mat',name),...
        'gt_labelfilename',sprintf('gt_labeled%s.mat',name),...
        'scorefilename',sprintf('scores_%s.mat',name),...
        'perframedir','perframe',...
        'rootoutputdir','',...
        'featureparamfilename','',...
        'featureconfigfile',fullfile('params','featureConfig.xml'));
      defaultConfig.plot.trx = struct('colormap','jet',...
        'colormap_multiplier','.7',...
        'extra_linestyle','-',...
        'extra_marker','.',...
        'extra_markersize','12');
      defaultConfig.plot.labels = struct('colormap','line',...
        'linewidth','3');
      defaultConfig.perframe.params = struct('fov',3.1416,'nbodylengths_near',2.5,...
        'thetafil',[0.0625,0.25,0.375,0.25,0.0625],...
        'max_dnose2ell_anglerange',127);
      defaultConfig.perframe.landmark_params = struct(...
        'arena_center_mm_x',0,'arena_center_mm_y',0,'arena_radius_mm',60,'arena_type','circle');
      
    end
    
    function projlist = GetProjectList(obj)
      if numel(obj.projparams)>0  
        projlist = {obj.projparams(:).name};
      else
        projlist = {};
      end
    end
    
    function projname = GetCurProjName(obj)
      projname = obj.projparams(obj.curproj).name;
    end
    
    function configfile = GetCurProjConfigfile(obj)
      configfile = obj.projparams(obj.curproj).configfile;
    end
    
    function SetCurrentProject(obj,curproj)
      obj.curproj = curproj;
    end
    
    function curproj = GetCurrentProject(obj)
      curproj = obj.curproj;
    end
    
    function SetName(obj,projnum)
      name = obj.projparams(projnum).name;
      obj.projparams(projnum).config.behaviors.names = name;
      obj.projparams(projnum).config.file.labelfilename = sprintf('labeled%s.mat',name);
      obj.projparams(projnum).config.file.gt_labelfilename = sprintf('gt_labeled%s.mat',name);
      obj.projparams(projnum).config.file.scorefilename = sprintf('scores_%s.mat',name);
      [dd,~] = fileparts(obj.projparams(projnum).configfile);
      obj.projparams(projnum).config.file.featureparamfilename = ...
        fullfile(dd,sprintf('WindowFeatures_%s.xml',name));
    end
    
    function CheckClash(obj,projnum)
      str = '';
      for ndx = 1:numel(obj.projparams)
        if ndx == projnum; continue; end
        if isfield(obj.projparams(ndx).config.file,'labelfilename') && ...
            isfield(obj.projparams(projnum).config.file,'labelfilename') && ...
            strcmp(obj.projparams(ndx).config.file.labelfilename,...
            obj.projparams(projnum).config.file.labelfilename),
          str = sprintf('%s, Labelfilename is same as for project %s',str, obj.projparams(ndx).name);
        end
        if isfield(obj.projparams(ndx).config.file,'gt_labelfilename') && ...
            isfield(obj.projparams(projnum).config.file,'gt_labelfilename') && ...
        strcmp(obj.projparams(ndx).config.file.gt_labelfilename ,...
            obj.projparams(projnum).config.file.gt_labelfilename ),
          str = sprintf('%s, Labelfilename is same as for project %s',str, obj.projparams(ndx).name);
        end
        if isfield(obj.projparams(ndx).config.file,'scorefilename') && ...
            isfield(obj.projparams(projnum).config.file,'scorefilename') && ...
            strcmp(obj.projparams(ndx).config.file.scorefilename,...
            obj.projparams(projnum).config.file.scorefilename),
          str = sprintf('%s, Labelfilename is same as for project %s',str, obj.projparams(ndx).name);
        end
        if isfield(obj.projparams(ndx).config.file,'featureparamfilename') &&...
            isfield(obj.projparams(projnum).config.file,'featureparamfilename') && ...
            strcmp(obj.projparams(ndx).config.file.featureparamfilename,...
            obj.projparams(projnum).config.file.featureparamfilename),
          str = sprintf('%s, Labelfilename is same as for project %s',str, obj.projparams(ndx).name);
        end
      end
        
      if ~isempty(str),
        uiwait(warndlg(sprintf('%s. Please change the names of these file to avoid overwriting files.',str)));
      end
    end
    
    function AddProject(obj,name,configfile,copyconfigFile)
      
      obj.projparams(end+1).name = name;
      obj.projparams(end).configfile = configfile;
      obj.projparams(end).save = true;
      obj.curproj = numel(obj.projparams);
      
      fileToRead = '';
      
      if exist(configfile,'file')
        wstr = 'The configuration file %s for this project already exists';
        wstr = sprintf('%s\n, Overwrite or use exisiting settings for the file?',wstr);
        in = questdlg(wstr,'Overwrite','Overwrite','Keep and use existing settings');
        switch in,
          case 'Overwrite',
            delete(configfile);
            if ~isempty(copyconfigFile) && exist(copyconfigFile,'file')
              fileToRead = copyconfigFile;
            end
          case 'Keep and use existing settings',
            fileToRead = configfile;
            obj.projparams(end).save = false;
        end
      elseif ~isempty(copyconfigFile) && exist(copyconfigFile,'file')
        fileToRead = copyconfigFile;
      end
      
      if ~isempty(fileToRead)
        curparams = ReadXMLParams(fileToRead);
        if isfield(curparams,'featureparamlist'),
          obj.projparams(end).pfList = fieldnames(curparams.featureparamlist);
          curparams = rmfield(curparams,'featureparamlist');
        else
          obj.projparams(end).pfList = [];
        end
        obj.projparams(end).config = curparams;
      else
        obj.projparams(end).config = obj.DefaultConfig(name);
        obj.projparams(end).pfList = [];
        obj.projparams(end).save = true;
      end
      
      obj.SetName(numel(obj.projparams));
      obj.CheckClash(numel(obj.projparams));
    
    end
    
    function RemoveProject(obj,projnum)
      if nargin<2,
        projnum = obj.curproj;
      end
      
      obj.projparams(projnum) = [];
      if obj.curproj > numel(obj.projparams)
        obj.curproj = obj.curproj -1;
        if obj.curproj < 1
          obj.curproj = [];
        end
      end
    end
    
    function WriteProjectManager(obj)
      docNode = com.mathworks.xml.XMLUtils.createDocument('behaviors');
      toc = docNode.getDocumentElement;
      for ndx = 1:numel(obj.projparams)
        curN.configfile = obj.projparams(ndx).configfile;
        toc.appendChild(createXMLNode(docNode,obj.projparams(ndx).name,...
          curN));
      end
      xmlwrite(obj.projfile,docNode);
    end
    
    function config = GetProjConfig(obj,projnum)
      if nargin<2
        projnum = obj.curproj;
      end
      config = obj.projparams(projnum).config;
    end
    
    function [data,success] = GetConfigAsTable(obj)
      success = true;
      if isempty(obj.curproj);
        data = {}; 
        return;
      end
      configparams = obj.projparams(obj.curproj).config;
      data = obj.addToList(configparams,{},'');
      idx = cellfun(@iscell,data(:,2));
      if any(idx),
        for i = find(idx(:)'),
          if all(cellfun(@ischar,data{i,2})),
            data{i,2} = sprintf('%s,',data{i,2}{:});
            if data{i,2}(end) == ',',
              data{i,2} = data{i,2}(1:end-1);
            end
          end
        end
      end
        
      if any(cellfun(@iscell,data(:,2))),
        data = {}; success = false;
        return;
      end
      
    end
    
    function list = addToList(obj,curStruct,list,pathTillNow)
      if isempty(fieldnames(curStruct)), return; end
      fnames = fieldnames(curStruct);
      for ndx = 1:numel(fnames)
        if isstruct(curStruct.(fnames{ndx})),
          list = obj.addToList(curStruct.(fnames{ndx}),list,[pathTillNow fnames{ndx} '.']);
        else
          list{end+1,1} = [pathTillNow fnames{ndx}];
          param = curStruct.(fnames{ndx});
          if isnumeric(param)
            q = num2str(param(1));
            for jj = 2:numel(param)
              q = [q ',' num2str(param(jj))];
            end
            list{end,2} = q;
          else
            list{end,2} = param;
          end
        end
      end
    end
      
    function success  = AddConfig(obj,name,value)
      
      if ~isempty(str2num(value))
        value = str2num(value);
      end
      
      success = true;
      
      iname = fliplr(name);
      curstruct = obj.projparams(obj.curproj).config;
      while true
        [iname,lastfield] = splitext(iname);
        if isempty(lastfield)
          fexist = isfield(curstruct,iname);
          break;
        else
          fexist = isfield(curstruct,fliplr(iname(2:end)));
          if ~fexist, break; end
          curstruct = curstruct.(fliplir(iname(2:end)));
        end
        if fexist, success = false; return; end
      end
      
      eval(sprintf('obj.projparams(obj.curproj).config.%s = value;',name));
      obj.projparams(obj.curproj).save = true;
    end
    
    function RemoveConfig(obj,name)

      [fpath,lastfield] = splitext(name);
      if isempty(lastfield)
        obj.projparams(obj.curproj).config = ...
          rmfield(obj.projparams(obj.curproj).config,fpath);
      else
        evalStr = sprintf(...
          'obj.projparams(obj.curproj).config.%s = rmfield(obj.projparams(obj.curproj).config.%s,lastfield(2:end));',...
          fpath,fpath);
        eval(evalStr);
      end
      obj.projparams(obj.curproj).save = true;      
    end
    
    function EditConfigName(obj,oldName,newName)
      eval_str = sprintf(...
        'value = obj.projparams(obj.curproj).config.%s;',...
        oldName);
      eval(eval_str);
      obj.RemoveConfig(oldName);
      obj.AddConfig(newName,value);
      obj.projparams(obj.curproj).save = true;
    end
    
    function EditConfigValue(obj,name,value)
      eval_str = sprintf(...
        'obj.projparams(obj.curproj).config.%s = value;',...
        name);
      eval(eval_str);
      obj.projparams(obj.curproj).save = true;
    end
    
    function SetPerframeList(obj,pfList)
      obj.projparams(obj.curproj).pfList = pfList;
      obj.projparams(obj.curproj).save = true;
    end
    
    function pfList = GetPerframeList(obj)
      pfList = obj.projparams(obj.curproj).pfList;
    end
    
    function [allPfList selected missing]= GetAllPerframeList(obj)
      featureconfigfile = obj.projparams(obj.curproj).config.file.featureconfigfile;
      params = ReadXMLParams(featureconfigfile);
      allPfList = fieldnames(params.perframe);
      selected = false(numel(allPfList),1);
      missing = {};
      curpf = obj.projparams(obj.curproj).pfList;
      if ~isempty(curpf),
        for ndx = 1:numel(curpf),
          allndx = find(strcmp(curpf{ndx},allPfList));
          if isempty(allndx)
            missing{end+1} = curpf{ndx};
          else
            selected(allndx) = true;
          end
        end
      else
        missing = {};
        selected = true(numel(allPfList,1));
      end
      
    end
    
    function SaveConfig(obj,projnum)
      
      docNode = com.mathworks.xml.XMLUtils.createDocument('params');
      toc = docNode.getDocumentElement;
      topNode = obj.projparams(projnum).config;
      att = fieldnames(topNode);
      for ndx = 1:numel(att)
        toc.appendChild(createXMLNode(docNode,att{ndx},topNode.(att{ndx})));
      end
      if ~isempty(obj.projparams(projnum).pfList),
        pfStruct = struct;
        for ndx = 1:numel(obj.projparams(projnum).pfList)
          pfStruct.(obj.projparams(projnum).pfList{ndx}) = struct;
        end
        toc.appendChild(createXMLNode(docNode,'featureparamlist',pfStruct));
      end
      xmlwrite(obj.projparams(projnum).configfile,docNode);
      
    end
    
    function SaveAllConfig(obj)
      for pno = 1:numel(obj.projparams)
        if obj.projparams(pno).save
          obj.SaveConfig(pno);
        end
      end
    end
    
    function AddClassifier(obj,name)
    end
    
    function GetClassifier(obj,name)
    end
    
  end
  
end
