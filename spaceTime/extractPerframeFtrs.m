function extractPerframeFtrs(outdir,ftrs,stationary,flowname)

if ~exist(outdir,'dir'),
  mkdir( outdir);
end

ff = fields(ftrs);

% Initialize the struct for features of all the frames
for fnum = 1:numel(ff)
  curf = ff{fnum};
  if strcmp(curf,'hogftrs') ,
    pfname = 'hf';
  elseif strcmp(curf,'flowftrs') && stationary,
    pfname = [flowname 's'];
  elseif strcmp(curf,'flowftrs') && ~stationary,
    pfname = flowname;
  else
    error('Unknown feature type');
  end
  for yy = 1:size(ftrs.(curf){1},1)
    for xx = 1:size(ftrs.(curf){1},2)
      for oo = 1:size(ftrs.(curf){1},3)
        perframe = struct;
        perframe.units.num = {};
        perframe.units.den = {'s'};
        for fly = 1:numel(ftrs.(ff{fnum}))
          tt = ftrs.(ff{fnum}){fly}(yy,xx,oo,:);
          perframe.data{fly} = tt(:);
        end
        perframe_name = sprintf('%s_%02d_%02d_%d',pfname,yy,xx,oo);
        outfilename = fullfile(outdir,perframe_name);
        save(outfilename,'-struct','perframe','-v7.3');
      end
    end
  end
end
