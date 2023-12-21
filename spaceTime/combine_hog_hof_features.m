function feat = combine_hog_hof_features(feat_split,dotranspose)

if nargin > 1 && dotranspose,
  fns = fieldnames(feat_split);
  for i = 1:numel(fns),
    feat_split.(fns{i}) = feat_split.(fns{i})';
  end
end
featmult = 10;
nori = 8;

[nframes,nfeat_split_front] = size(feat_split.hog_front);
[nframes,nfeat_split_side] = size(feat_split.hog_side);
nlandmarks_front = nfeat_split_front/featmult^2/nori;
nlandmarks_side = nfeat_split_side/featmult^2/nori;
nlandmarks = nlandmarks_front+nlandmarks_side;

featsz = [featmult,featmult*nlandmarks,nori,2];
feat = zeros([nframes,featsz]);
idxside = 1:nlandmarks_side*featmult;
idxfront = nlandmarks_side*featmult+1:nlandmarks*featmult;

feat(:,:,idxfront,:,2) = reshape(feat_split.hog_front,[nframes,featmult,featmult*nlandmarks_front,nori]);
feat(:,:,idxside,:,2) = reshape(feat_split.hog_side,[nframes,featmult,featmult*nlandmarks_side,nori]);
feat(:,:,idxfront,:,1) = reshape(feat_split.hof_front,[nframes,featmult,featmult*nlandmarks_front,nori]);
feat(:,:,idxside,:,1) = reshape(feat_split.hof_side,[nframes,featmult,featmult*nlandmarks_side,nori]);

feat = feat(:,:);