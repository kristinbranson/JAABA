function DeleteBadPerFrameData(expdirs,perframefns)

if nargin < 2,
  perframefns = {'x_mm.mat','y_mm.mat','a_mm.mat','theta_mm.mat','b_mm.mat'};
end

for i = 1:numel(expdirs),
  if ~exist(expdirs{i},'dir'),
    continue;
  end
  if ~exist(fullfile(expdirs{i},'perframe'),'dir'),
    continue;
  end
  for j = 1:numel(perframefns),
    filename = fullfile(expdirs{i},'perframe',perframefns{j});
    if exist(filename,'file'),
      delete(filename);
    end
  end
  
end