classdef ExperimentTags
% Experiments can be tagged for custom identification/grouping. 
% Tags are stored in a data structure 'expTags'.

  methods (Static)
    
    function expTags = expTags(expdirs)
      % Empty expDirTag constructor 
      % This signature is a little silly
    
      expTags = cellfun(@(x)cell(1,0),expdirs,'uni',0);
    end
    
    function verifyTags(expTags,nexp)
      assert(iscell(expTags),...
        'Experiment tags must be a cell array of row cellstrs.');
      assert(all(cellfun(@(x)iscellstr(x)&&isrow(x),expTags)),...
        'Experiment tags must be a cell array of row cellstrs.');
      assert(numel(expTags)==nexp,...
        'Unexpected number of experiment tag specifications.');
    end      
    
    function uniqueTags = allUniqueTags(expTags)
      alltags = cat(2,expTags{:});
      uniqueTags = unique(alltags(:));
    end
    
    function tf = findTag(expTags,tag)
      tf = cellfun(@(x)any(strcmp(tag,x)),expTags);
    end
    
    function expTags = rmTag(expTags,tag)
      % remove tag from all expTags
      for i = 1:numel(expTags)
        tmptags = expTags{i};
        tf = strcmp(tmptags,tag);
        tmptags(:,tf) = [];
        expTags{i} = tmptags;
      end
    end
    
    function expTags = cleanLegacyNotTags(expTags)
      unTags = ExperimentTags.allUniqueTags(expTags);
      Ntags = numel(unTags);
      PFIX = 'not_';
      for i = 1:Ntags
        tag = unTags{i};
        if strncmp(tag,PFIX,numel(PFIX))
          tag0 = tag(numel(PFIX)+1:end);
          if ismember(tag0,unTags)
            tf = ExperimentTags.findTag(expTags,tag);
            tf0 = ExperimentTags.findTag(expTags,tag0);
            if all(xor(tf,tf0))
              expTags = ExperimentTags.rmTag(expTags,tag); 
              warningNoTrace('ExperimentTags:notTags',...
                'Removing legacy not-tag ''%s''.',tag);
            else
              % Matching tags 'not_<tag>' and '<tag>' both exist; but they
              % are not perfect comlpements
              warningNoTrace('ExperimentTags:notTags',...
                'Experiment tag ''%s'' appears to be a legacy not-tag, but it is not the complement of tag ''%s''. Leaving as-is.',...
                tag,tag0);
            end
          end
        end
      end
    end
    
  end
  
end

