classdef Ethogram
  % TODO: to hold all ethogram/labelgram-related 
  
  methods (Static)    
    
    function boutExport(xlsfile,boutmat,expnames,behnames)      
      Nexp = numel(expnames);
      Nbeh = numel(behnames);
      assert(isequal(size(boutmat),[Nexp Nbeh 2]));
      
      nobehnames = Labels.noneOrNoBehaviorNames(behnames);
      
      fh = fopen(xlsfile,'w');
      assert(fh~=-1,'Cannot open file ''%s'' for writing.',xlsfile);
      
      for i = 1:Nexp
        fprintf(fh,'%s\n',expnames{i});
        for j = 1:Nbeh
          fprintf(fh,'\t%s\n',behnames{j});
          posbouts = boutmat{i,j,1};
          Ethogram.hlpExportBout(fh,posbouts);
          
          fprintf(fh,'\t%s\n',nobehnames{j});
          negbouts = boutmat{i,j,2};
          Ethogram.hlpExportBout(fh,negbouts);
        end
      end
      
      fclose(fh);      
    end
    function hlpExportBout(fh,bouts)
      Nbout = size(bouts,1);
      assert(size(bouts,2)==2);
      for k = 1:Nbout
        fprintf(fh,'\t\t%d\t%d\n',bouts(k,1),bouts(k,2));
      end
    end
    
  end
    
end