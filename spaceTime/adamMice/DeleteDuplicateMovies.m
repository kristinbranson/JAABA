function dupdirsdeleted = DeleteDuplicateMovies(rootdir)

dupdirsdeleted = {};


dd = dir(rootdir);
dd = dd([dd.isdir]);
for ndx =1:numel(dd)
  if strcmp(dd(ndx).name(1),'.'), continue; end
  isduplicate = false;
  if strcmp(dd(ndx).name,'movie'),
    firstmovie = fullfile(rootdir,'movie.avi');
    dupmovie = fullfile(rootdir,'movie','movie.avi');
    if exist(firstmovie,'file') && exist(dupmovie,'file'),
      [readframe1,nframes1] = get_readframe_fcn(firstmovie);
      [readframe2,nframes2] = get_readframe_fcn(dupmovie);
      if nframes1 == nframes2,
        im1 = readframe1(1);
        im2 = readframe2(1);
        if all(im1(:) == im2(:)),
          isduplicate = true;
        end
      end
      if isduplicate,
        fprintf('Deleting directory %s, as it is a duplicate of %s\n',fullfile(rootdir,'movie'),fullfile(rootdir,'movie.avi'));
        cmd = sprintf('rm -r %s',fullfile(rootdir,'movie'));
        if true,
          fprintf('%s\n',cmd);
          unix(cmd);
        else
          fprintf('%s\n',cmd);
        end
        dupdirsdeleted{end+1} = fullfile(rootdir,'movie');
        continue;
      else
        fprintf('Not deleting %s, movies do not match\n',fullfile(rootdir,'movie'));
      end
    else
      fprintf('Not deleting %s, both movies do not exist\n',fullfile(rootdir,'movie'));
    end
  end

  dupdirsdeleted = [dupdirsdeleted,DeleteDuplicateMovies(fullfile(rootdir,dd(ndx).name))];
  
end