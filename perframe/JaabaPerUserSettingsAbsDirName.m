function jaabaPerUserDirName=JaabaPerUserSettingsAbsDirName()

if isunix() || ismac()
  settingsDirName=getenv('HOME');
  jaabaPerUserDirName=fullfile(settingsDirName,'.jaaba');
else
  % Windows
  settingsDirName=getenv('APPDATA');  
  jaabaPerUserDirName=fullfile(settingsDirName,'jaaba');
end
if ~exist(jaabaPerUserDirName,'dir');
  mkdir(jaabaPerUserDirName);
end

end
