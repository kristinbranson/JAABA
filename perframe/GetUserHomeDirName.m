function homeDirName=GetUserHomeDirName()

if isunix() || ismac()
  homeDirName=getenv('HOME');
else
  % Windows
  homeDriveName=getenv('%HOMEDRIVE%');  
  homeDirName=getenv('%HOMEDIR%');  
  homeDirName=fullfile(homeDriveName,homeDirName);
end

end
