function trx = apt2trx(apttrk,bodylandmarks)

if isempty(bodylandmarks),
  bodylandmarks = 1:(apttrk.npts);
end

nflies = apttrk.ntlts;
trx = [];
for fly = 1:nflies,
  ff = apttrk.startframes(fly);
  ef = apttrk.endframes(fly);
  xy = apttrk.getPTrkTgt(fly);
  trxcurr = struct;
  trxcurr.x = mean(xy(bodylandmarks,1,ff:ef),1);
  trxcurr.y = mean(xy(bodylandmarks,2,ff:ef),1);
  dx = xy(bodylandmarks(2),1,ff:ef)-xy(bodylandmarks(1),1,ff:ef);
  dy = xy(bodylandmarks(2),2,ff:ef)-xy(bodylandmarks(1),2,ff:ef);
  trxcurr.theta = atan2( dy,dx );
  trxcurr.a = sqrt(dx.^2+dy.^2)/4;
  trxcurr.b = trxcurr.a/3;
  trxcurr.a_mm = trxcurr.a;
  trxcurr.b_mm = trxcurr.b;
  trxcurr.theta_mm = trxcurr.theta;
  trxcurr.x_mm = trxcurr.x;
  trxcurr.y_mm = trxcurr.y;
  trxcurr.firstframe = ff;
  trxcurr.endframe = ef;
  trxcurr.nframes = ef-ff+1;
  trxcurr.off = 1-ff;
  trxcurr.kpts = reshape(xy(:,:,ff:ef),[size(xy,1)*size(xy,2),ef-ff+1]);
  trxcurr.trkInfo = apttrk.trkInfo;
  trxcurr.dt = ones(1,ef-ff);
  trxcurr.pxpermm = 1;
  trx = structappend(trx,trxcurr);
 
end