try
  uiwait(PrepareJAABAData);
  
catch ME
  uiwait(errordlg(getReport(ME),'Error running PrepareJAABAData'));
end

delete(findall(0,'type','figure'));