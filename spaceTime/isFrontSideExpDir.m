function tf = isFrontSideExpDir(edir)
moviefrt = fullfile(edir,'movie_frt.avi');
moviesde = fullfile(edir,'movie_sde.avi');
moviecmb = fullfile(edir,'movie_comb.avi');
tf = (exist(moviefrt,'file')==2 && exist(moviesde,'file')==2) || exist(moviecmb,'file')==2;