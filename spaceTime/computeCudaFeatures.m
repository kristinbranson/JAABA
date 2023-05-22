function computeCudaFeatures(pathtocudascript, exp_list, varargin)

  [bash_path, compute_features, combine_movie] = myparse(varargin,...
                'bash_path', '"C:\Program Files\Git\bin\bash.exe"',...
                'compute_features', '1',...
                'combine_movie', '1');

  fullpathtoscript = pathtocudascript;
  runcommandtoscript = [bash_path ' ' fullpathtoscript ' ' exp_list ' ' compute_features ' ' combine_movie];
  disp(runcommandtoscript)
  try
     system(runcommandtoscript);
  catch
      fclose('all');
  end
  %cmdStr = ['"C:\Program Files\Git\bin\bash.exe" C:\Users\27rut\BIAS\scripts\run_bias.sh Y:\hantman_data\jab_experiments\STA14\STA14\20230503\hantman_modified.txt'];

  % convert cuda features to features.mat
  convertFeaturestoMatlab(exp_list)
  
end