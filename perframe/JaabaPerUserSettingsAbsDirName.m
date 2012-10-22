function jaabaPerUserDirName=JaabaPerUserSettingsAbsDirName()

if isunix() || ismac()
  settingsDirName=getenv('HOME');
  jaabatxt = '.jaaba.txt';
  jaabadir = '.jaaba';
else
  settingsDirName = getenv('APPDATA');
  jaabatxt = 'jaaba.txt';
  jaabadir = 'jaaba';
end
if isdeployed,
  if exist(fullfile(settingsDirName,jaabatxt),'file')
    ff = fopen(fullfile(settingsDirName,jaabatxt));
    jaabaPerUserDirName = fgetl(ff);
    fclose(ff);
  else
    jaabaPerUserDirName = uigetdir(settingsDirName,'Select a file to store JAABA configuration files');
    if ~exist(jaabaPerUserDirName)
      mkdir(jaabaPerUserDirName);
    end
    ff  = fopen(fullfile(settingsDirName,jaabatxt),'w');
    fprintf(ff,'%s',jaabaPerUserDirName);
    fclose(ff);
    if jaabaPerUserDirName == 0,
      error('Cannot continue without a settings dir');
    end
    dd = dir(fullfile(pwd,'params','*.xml'));
    for ndx = 1:numel(dd)
      copyfile(fullfile(pwd,'params',dd(ndx).name),....
        fullfile(jaabaPerUserDirName,dd(ndx).name));
    end
  end
else
  jaabaPerUserDirName=fullfile(settingsDirName,jaabadir);
end

if ~exist(jaabaPerUserDirName,'dir');
  mkdir(jaabaPerUserDirName);
end

end
