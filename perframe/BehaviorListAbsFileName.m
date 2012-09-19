function file_name=BehaviorListAbsFileName()

% figure out a good place to store BehaviorList.xml
jaabaPerUserDirName=JaabaPerUserSettingsAbsDirName();
file_name=fullfile(jaabaPerUserDirName,'BehaviorList.xml');

end
