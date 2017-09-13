classdef JLabelInterf 
  % "External" interface to MouseJAABA
  
  methods (Static)

    % StartJAABA conceptually overlaps with this API
    
    function openJabExpFrame(jabfile,expdir,frameno) 
      h = findall(0,'type','figure','name','JAABA');
      if ~isempty(h)
        JLabelH = guidata(h);
        jlJabfile = JLabelH.data.everythingFileNameAbs;
        if strcmp(jlJabfile,jabfile)
          % JLabel is open and the current project is the right one
          JLabelInterf.setJLabelMovieFrame(h,expdir,frameno);
          return;
        else
          JLabel('figure_JLabel_CloseRequestFcn',h,[],guidata(h));
        end
      end
      
      JLabel('jabfile',jabfile);
      h = findall(0,'type','figure','name','JAABA');
      if isempty(h),
        error('Could not open JAABA');
      end
      JLabelInterf.setJLabelMovieFrame(h,expdir,frameno);
    end
    
    function setJLabelMovieFrame(hJL,expdir,frameno)
      JLabelH = guidata(hJL);
      jexpnum = find(strcmp(JLabelH.data.expdirs,expdir));
      if ~isempty(jexpnum)
        JLabelH = JLabel('SetCurrentMovie',JLabelH,jexpnum);
        guidata(hJL,JLabelH); % necessary?
        JLabel('SetCurrentFrame',JLabelH,1,frameno,hJL);
      else
        warning('Experiment ''%s'' does not exist in current project.',expdir);        
%         JLabelH.data.AddExpDir(expdir);
%         JLabelH = JLabel('SetCurrentMovie',JLabelH,JLabelH.data.nexps);
      end
    end
    
  end

end