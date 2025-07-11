function trx = addAPTTrk2Trx(trx,trkfile,varargin)

[view,aptInfo,prev_width] = myparse(varargin,'view',1,'aptInfo',[],'prev_width',0);
if ~exist(trkfile, 'file')
  error('Trk file %s does not exist', trkfile) ;
end
trk = TrkFile.load(trkfile);

for fly = 1:numel(trx)
  %                     trx_ndx = find(trk.pTrkiTgt==i);
  %                     if isempty(trx_ndx)
  %                       msg = sprintf('APT trk %s does not tracking results for animal %d',trkfilename,i);
  %                       success = false;
  %                       return;
  %                     end
  cur_trk = trk.getPTrkTgt(fly);
  ff = trx(fly).firstframe;
  ef = trx(fly).endframe;
  cur_trk = cur_trk(:,:,ff:ef);
  cur_trk(:,1,:) = cur_trk(:,1,:) + prev_width;
  if view == 1
    trx(fly).kpts = reshape(cur_trk,[size(cur_trk,1)*size(cur_trk,2),ef-ff+1]);
    trx(fly).aptInfo = aptInfo;
    trx(fly).trkInfo = trk.trkInfo;
  end
  trx(fly).(sprintf('kpts_view%d',view)) = reshape(cur_trk,[size(cur_trk,1)*size(cur_trk,2),ef-ff+1]);
  %                     trx(fly).apt_trk{ndx} = cur_trk;
end

