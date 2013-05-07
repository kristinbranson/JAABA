function successAll=testAll()

thisFileNameAbs=mfilename('fullpath');  % without .m, for some reason
thisDirNameAbs=fileparts(thisFileNameAbs);
thisFunctionName=fileNameRelFromAbs(thisFileNameAbs);
d=dir(fullfile(thisDirNameAbs,'*.m'));
testFileNames={d.name}';
% Get rid of backup and hidden files
testFileNames=cellFilter(@(str)(isempty(regexp(str,'^\.','once'))),testFileNames);
testFunctionNames=cellfun(@baseNameFromFileName,testFileNames,'UniformOutput',false);
testFunctionNames(strcmp(thisFunctionName,testFunctionNames))=[];  
  % eliminate the name of this file from the list
nTests=length(testFunctionNames);
success=false(size(testFunctionNames));
for i=1:nTests
  testFunctionName=testFunctionNames{i};
  fprintf('%s: \n',testFunctionName);
  try
    successThis=feval(testFunctionName);
    if successThis ,
      fprintf('%s succeded.\n\n',testFunctionName);
    else
      fprintf('%s Failed.\n\n',testFunctionName);
    end
  catch excp
      successThis=false;      
      fprintf('%s failed with uncaught exception: id: %s, message: %s.\n\n', ...
              testFunctionName, ...
              excp.identifier, ...
              excp.message);
  end
  success(i)=successThis;
end
nSucceeded=sum(success);
fprintf('Summary: %d of %d tests succeeded.\n',nSucceeded,nTests);
successAll=(nSucceeded==nTests);

end
