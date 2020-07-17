function extractSTPerframeFtrs(outdir,ftrs,stationary,flowname,params)

if ~exist(outdir,'dir')
  mkdir( outdir);
end

ff = fields(ftrs);

% Initialize the struct for features of all the frames
for fnum = 1:numel(ff)
  curf = ff{fnum};
  if strcmp(curf,'hogftrs') 
    pfname = 'hf';
  elseif strcmp(curf,'flowftrs') && stationary
    pfname = [flowname 's'];
  elseif strcmp(curf,'flowftrs') && ~stationary
    pfname = flowname;
  else
    error('Unknown feature type');
  end
  ny = size(ftrs.(curf){1},1);
  nx = size(ftrs.(curf){1},2);
  nb = size(ftrs.(curf){1},3);
  psz = params.psize;
  for yy = 1:ny
    for xx = 1:nx
      for oo = 1:nb
        perframe = struct;
        perframe.units.num = {};
        perframe.units.den = {'s'};
        for fly = 1:numel(ftrs.(ff{fnum}))
          tt = ftrs.(ff{fnum}){fly}(yy,xx,oo,:);
          perframe.data{fly} = tt(:);
        end
        perframe_name = sprintf('st_%s_%02d_%02d_%d_ny%d_nx%d_ns%d',pfname,yy,xx,oo,ny,nx,psz);
        perframe.params = params;
        outfilename = fullfile(outdir,perframe_name);
        save(outfilename,'-struct','perframe','-v7.3');
      end
    end
  end
end
