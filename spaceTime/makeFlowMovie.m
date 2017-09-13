function success = makeFlowMovie(inMovie,outFile)

clipat = 25;

success =false;
try
  [readframefcn,~,fid,headerinfo] = get_readframe_fcn(inMovie);
catch exception,
  warning('Could not open %s for computing the features %s\n',vidPath,exception.message); %#ok<*WNTAG>
  return;
end

writerObj = VideoWriter(sprintf('%s.avi',outFile));
writerObj.FrameRate = 10;  
writerObj.Quality = 95;    
open(writerObj);           

h = waitbar(0,'Computing Flow 0%% done');

Iprev = readframefcn(1);
t1 = tic;
for ndx = 2:headerinfo.nframes
  Icur = readframefcn(ndx);
  
  [Vx,Vy,~] = optFlowLk(Iprev,Icur,[],3);
  Iwrite = cat(3,Vx,Vy);
  Iwrite(:,:,3) = 0;
  
  Iwrite(Iwrite>clipat) = clipat;
  Iwrite(Iwrite<-clipat) = -clipat;
  Iwrite = Iwrite+clipat;
  Iwrite = uint8(Iwrite/2/clipat*255);
  writeVideo(writerObj,Iwrite);
  Iprev = Icur;
  if mod(ndx,20)
    t2 = toc(t1);
    wmsg{1} = sprintf('Computing Flow %.2f%% done',ndx/headerinfo.nframes*100);
    wmsg{2} = sprintf('Elapsed Time:%.2f',t2);
    wmsg{3} = sprintf('Time Remaining:%.2f',(headerinfo.nframes-ndx)/ndx*t2);
    waitbar(ndx/headerinfo.nframes,h,wmsg);
  end
end

close(h);
close(writerObj);

if fid>0,
  fclose(fid);
end

success = true;