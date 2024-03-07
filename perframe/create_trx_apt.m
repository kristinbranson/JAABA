function [trx,success,msg] = create_trx_apt(trkfilename, aptInfo)

success = true; msg = '';
if ischar(trkfilename)
  trkfilename = {trkfilename};
end

if aptInfo.is_ma
  assert(exist(trkfilename{1},'file')>0,sprintf('Trkfile %s does not exist',trkfilename{1}));
  trk = TrkFile.load(trkfilename{1});
  trx= apt2trx(trk,aptInfo.ma_head_tail);
  return;
end

trx = struct();
for view = 1:numel(trkfilename)
  if ~exist(trkfilename{view},'file')
    success = false;
    msg = sprintf('Trkfile %s doesnt exist',trkfilename{view});
    return;
  end
  trk = TrkFile.load(trkfilename{view});
%   frms = trk.pTrkFrm;
%   dd = frms(2:end)-frms(1:end-1);
%   assert(all(dd==1),'Frames in trkfile should be in order')
  assert(size(trk.pTrk,4)==1, 'Cannot create trx file for trk with multiple animals');
  if aptInfo.has_trx
    warning('APT project already has trx. This might overwrite it');
  end

  temp_pts = trk.getPTrkTgt(1);
  aa = trk.isalive(1:size(temp_pts,3),1);
  ff = find(aa,1,'first');
  ef = find(aa,1,'last');
  temp_pts = temp_pts(:,:,ff:ef);
  t = ones(1,ef-ff+1);
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
      temp = temp_pts(sel_pts,:,:);
      temp = mean(temp,1);
      temp = shiftdim(temp,1);
      x = temp(1,:);
      y = temp(2,:);
      if aptInfo.apt_trx_orient < 1
        theta = t*0;
      else
        t_pt = squeeze(temp_pts(aptInfo.apt_trx_orient,:,:));
        theta = atan2(t_pt(1,:)-x,t_pt(2,:)-y);
      end

  end

  trx.(sprintf('x_view%d',view)) = x;
  trx.(sprintf('y_view%d',view)) = y;
  trx.(sprintf('theta_view%d',view)) = theta;
  trx.(sprintf('kpts_view%d',view)) = reshape(temp_pts,[size(temp_pts,1)*size(temp_pts,2),ef-ff+1]);

  if view == 1
    trx.pxpermm = 1;
    trx.nframes = ef-ff+1;
    trx.firstframe = ff;
    trx.endframe = ef;
    trx.id = 0;
    trx.trkfile = trkfilename;
    trx.from_apt = true;
    trx.aptInfo = aptInfo;
    trx.off = 1-trx.firstframe;
    trx.dt = ones(1,ef-ff);
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
    trx.kpts = reshape(temp_pts,[size(temp_pts,1)*size(temp_pts,2),ef-ff+1]);
    trx.trkInfo = trk.trkInfo;
  end
end