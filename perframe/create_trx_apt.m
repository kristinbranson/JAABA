function trx = create_trx_apt(trkfilename, aptInfo)

if ischar(trkfilename)
  trkfilename = {trkfilename};
end

trx = struct();
for view = 1:numel(trkfilename)
  trk = load(trkfilename{view},'-mat');
  frms = trk.pTrkFrm;
  dd = frms(2:end)-frms(1:end-1);
  assert(all(dd==1),'Frames in trkfile should be in order')
  assert(size(trk.pTrk,4)==1, 'Cannot create trx file for trk with multiple animals');
  if aptInfo.has_trx
    warning('APT project already has trx. This might overwrite it');
  end

  t = ones(1,numel(frms));
  switch aptInfo.apt_trx_type
    case 'crop'
      crop_loc = trk.trkInfo.crop_loc;
      x = mean(crop_loc(1:2))*t;
      y = mean(crop_loc(3:4))*t;
      theta = t*0;
    case 'frame'
      h_imsz = double(cell2mat(trk.trkInfo.params.imsz))/2;
      x = h_imsz(2)*t;
      y = h_imsz(1)*t;
      theta = t*0;
    case 'trk'
      sel_pts = aptInfo.apt_trx_centers;
      temp = trk.pTrk(sel_pts,:,:,1);
      temp = mean(temp,1);
      temp = shiftdim(temp,1);
      x = temp(1,:);
      y = temp(2,:);
      if aptInfo.apt_trx_orient < 1
        theta = t*0;
      else
        t_pt = squeeze(trk.pTrk(aptInfo.apt_trx_orient,:,:,1));
        theta = atan2(t_pt(1,:)-x,t_pt(2,:)-y);
      end
      
  end

  trx.(sprintf('x_view%d',view)) = x;
  trx.(sprintf('y_view%d',view)) = y;
  trx.(sprintf('theta_view%d',view)) = theta;

  if view == 1
    trx.pxpermm = 1;
    trx.nframes = numel(frms);
    trx.firstframe = min(frms);
    trx.endframe = max(frms);
    trx.id = 0;
    trx.trkfile = trkfilename;
    trx.from_apt = true;
    trx.aptInfo = aptInfo;
    trx.off = 1-trx.firstframe;
    trx.dt = ones(1,numel(frms)-1);
    trx.x = x;
    trx.y = y;
    trx.x_mm = trx.x;
    trx.y_mm = trx.y;
    trx.theta = theta;
    trx.theta_mm = theta;
    trx.pxpermm = 1;
    trx.a = t;
    trx.b = t;
    trx.a_mm = t;
    trx.b_mm = t;
  end
end