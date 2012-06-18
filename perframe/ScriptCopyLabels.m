%% copy labels

% rootinputdir = 'B:\robiea\Projects_data\JAABA\GroundTruth';
% rootoutputdir = 'C:\Data\JAABA\groundtruth_pBDPGAL4U_data';
rootinputdir = 'B:\robiea\Projects_data\JAABA\Data_Righting';
rootoutputdir = 'C:\Data\JAABA\FlyBowl';

%labelfilestrs = {'gt_labeledWalks.mat','labeledWalks.mat'};
%labelfilestrs = {'labeledImpJumps_gt.mat','labeledImpJumps.mat'};
%labelfilestrs = {'labeledpivot_tail.mat','gt_labeledpivot_tail.mat','scores_pivot_tail.mat'};
%labelfilestrs = {'labeledStops.mat','gt_labeledStops.mat','scores_Stops.mat'};
labelfilestrs = {'labeledRighting.mat','labeledRighting_gt.mat','scores_Righting.mat'};
doforce = false;
doignore = false;

%%

exps = dir(fullfile(rootinputdir,'*201*T*'));
exps(~[exps.isdir]) = [];

for i = 1:numel(exps),
  inexpdir = fullfile(rootinputdir,exps(i).name);
  outexpdir = fullfile(rootoutputdir,exps(i).name);
  if ~exist(outexpdir,'dir'),
    fprintf('Output experiment directory %s does not exist, skipping.\n',outexpdir);
    continue;
  end
  for j = 1:numel(labelfilestrs),
    infile = fullfile(inexpdir,labelfilestrs{j});
    outfile = fullfile(outexpdir,labelfilestrs{j});
    if ~exist(infile,'file'),
      fprintf('File %s does not exist, skipping\n',infile);
      continue;
    end

    if exist(outfile,'file'),
      if doignore,
        fprintf('Skipping %s -> %s\n',infile,outfile);
        continue;
      end
      if ~doforce,
        res = questdlg(sprintf('File %s exists, overwrite with %s?',outfile,infile));
        if ismember(lower(res),{'no','cancel'}),
          fprintf('Skipping %s -> %s\n',infile,outfile);
          continue;
        end
      end
    end
    fprintf('Copying %s -> %s\n',infile,outfile);
    copyfile(infile,outfile);
  end
end
