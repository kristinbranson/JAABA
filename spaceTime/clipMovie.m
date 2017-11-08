function clipMovie(indir,outdir,startFrame,endFrame,varargin)
% function clipMovie(indir,outdir,startFrame,endFrame,...)

[moviename,trxfilename] = myparse(varargin,...
  'moviename','movie.ufmf','trxfilename','trx.mat');

if ~exist(outdir,'dir'),
  mkdir(outdir,'perframe');
end

if ~exist(fullfile(outdir,moviename),'file')
  copyfile(fullfile(indir,moviename),fullfile(outdir,moviename));
end

if ~strcmp(trxfilename(end-3:end),'.mat')
  trxfilename = [trxfilename '.mat'];
end

Q = load(fullfile(indir,trxfilename));
selflds = {'x','y','theta','a','b','timestamps',...
  'x_mm','y_mm','theta_mm','a_mm','b_mm'};
pfstart = [];
pfend = [];
for ndx = 1:numel(Q.trx)
  if Q.trx(ndx).firstframe > endFrame || Q.trx(ndx).endframe < startFrame
    Q.trx(ndx).firstframe = startFrame;
    Q.trx(ndx).endframe = startFrame;
    for fnum = 1:numel(selflds)
      Q.trx(ndx).(selflds{fnum}) = Q.trx(ndx).(selflds{fnum})(1);
    end
    Q.trx(ndx).dt = 1;
    Q.trx(ndx).nframes = 1;
    pfstart(ndx) = 1; pfend(ndx) = 1;
  else
    curStart = max(Q.trx(ndx).firstframe,startFrame);
    curEnd = min(Q.trx(ndx).endframe,endFrame);
    for fnum = 1:numel(selflds)
      Q.trx(ndx).(selflds{fnum}) = Q.trx(ndx).(selflds{fnum})( (curStart:curEnd)-Q.trx(ndx).firstframe+1);
    end
    Q.trx(ndx).dt = Q.trx(ndx).dt( (curStart:curEnd-1)-Q.trx(ndx).firstframe+1);
    Q.trx(ndx).endframe = curEnd;
    Q.trx(ndx).nframes = curEnd-curStart+1;
    Q.trx(ndx).firstframe = curStart;
    Q.trx(ndx).off = -curStart+1;
    pfstart(ndx) = curStart-Q.trx(ndx).firstframe+1;
    pfend(ndx) = curEnd-Q.trx(ndx).firstframe+1;
  end
end
if ~isempty(Q.timestamps),
  Q.timestamps = Q.timestamps(startFrame:endFrame);
end
save(fullfile(outdir,trxfilename),'-struct','Q');

if ~exist(fullfile(outdir,'perframe'),'dir'),
  mkdir(fullfile(outdir,'perframe'));
end

dd = dir(fullfile(indir,'perframe','*.mat'));
for pf = 1:numel(dd)
  Q = load(fullfile(indir,'perframe',dd(pf).name));
  for ndx = 1:numel(Q.data)
    if isnan( pfstart(ndx))
      Q.data{ndx} = [];
    else
      Q.data{ndx} = Q.data{ndx}(pfstart(ndx):pfend(end));
    end
  end
  save(fullfile(outdir,'perframe',dd(pf).name),'-struct','Q');
  
end
