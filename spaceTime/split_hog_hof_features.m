function matlabfeat_split = split_hog_hof_features(featin)

[nframes,nfeat] = size(featin);
nlandmarks = nfeat/100/8/2;
featsz = [10,10*nlandmarks,8,2];
feat = reshape(featin,[nframes,featsz]);
matlabfeat_split = struct;
matlabfeat_split.hog_front = reshape(feat(:,:,26:end,:,2),[nframes,nfeat/4]);
matlabfeat_split.hog_side = reshape(feat(:,:,1:25,:,2),[nframes,nfeat/4]);
matlabfeat_split.hof_front = reshape(feat(:,:,26:end,:,1),[nframes,nfeat/4]);
matlabfeat_split.hof_side = reshape(feat(:,:,1:25,:,1),[nframes,nfeat/4]);

