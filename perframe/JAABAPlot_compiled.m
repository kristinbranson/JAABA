try
  uiwait(JAABAPlot);
  
catch ME
  errordlg(getReport(ME),'Error running JAABAPlot');
end

delete(findall(0,'type','figure'));