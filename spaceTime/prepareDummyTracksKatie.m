function prepareDummyTracksKatie(moviename)

expdir = fileparts(moviename);
[~,nframes,fid,~] = get_readframe_fcn(moviename);
if fid>0,
  fclose(fid);
end
T = struct;
T.timestamps = ones(1,nframes)*now;
T.trx.x = ones(1,nframes)*5;
T.trx.y = ones(1,nframes)*5;
T.trx.a = ones(1,nframes)*1;
T.trx.b = ones(1,nframes)*1;
T.trx.theta = zeros(1,nframes);
T.trx.id = 1;
T.trx.firstframe = 1;
T.trx.endframe = nframes;
T.trx.x_mm = T.trx.x;
T.trx.y_mm = T.trx.y;
T.trx.a_mm = T.trx.a;
T.trx.b_mm = T.trx.b;
T.trx.theta_mm = T.trx.theta;
T.trx.dt = 1/1300;
T.trx.fps = 1300;
T.trx.pxpermm = 1;
T.trx.sex = 'M';
T.trx.off = 0;
T.trx.nframes = nframes;

save(fullfile(expdir,'trx.mat'),'-struct','T');