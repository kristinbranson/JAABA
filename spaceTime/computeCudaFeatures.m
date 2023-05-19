function computeCudaFeatures(pathtocudascript, exp_list, varargin)

  [bash_path] = myparse(varargin,...
                'bash_path', '"C:\Program Files\Git\bin\bash.exe"');

  fullpathtoscript = pathtocudascript; %'C:\Users\27rut\BIAS\scripts\run_bias.sh';
  runcommandtoscript = [bash_path ' ' fullpathtoscript ' ' exp_list];
  disp(runcommandtoscript)
  system(runcommandtoscript);
  %cmdStr = ['"C:\Program Files\Git\bin\bash.exe" C:\Users\27rut\BIAS\scripts\run_bias.sh Y:\hantman_data\jab_experiments\STA14\STA14\20230503\hantman_modified.txt'];

  % convert cuda features to features.mat
  convertFeaturestoMatlab(exp_list)

end