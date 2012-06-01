%function convertmovie(filename)
%
%convert a .mov file to a .avi file as a work around to
%the slow load times of the former in matlab on linux
%
%if filename is 'movie.mov', output will be movie.mov.avi.

function convertmovie(filename)

tic;  readobj = VideoReader(filename);  toc;
writeobj = VideoWriter(filename,'Motion JPEG AVI');
%writeobj = VideoWriter(filename,'Uncompressed AVI');
%writeobj = VideoWriter(filename,'MPEG-4');
open(writeobj);

nFrames = readobj.NumberOfFrames;

for k = 1 : nFrames
  writeVideo(writeobj,read(readobj,k));
  if(rem(k,nFrames/100)<1)
    disp([num2str(100*k/nFrames,3) '% done']);
  end
end

close(writeobj);
