function catmov(mov1,mov2,outmov)
% catmov(mov1,mov2,outmov)
% Concatenate movies (frame by frame)
%
% mov1, mov2: movie filenames. Movies must have the same number of frames
% and have compatible dimensions (currently concat along x-dir).
% outmov: output filename
%

assert(exist(mov1,'file')==2,'Cannot find ''%s''',mov1);
assert(exist(mov2,'file')==2,'Cannot find ''%s''',mov1);

if exist(outmov,'file')
  warning('catmov:overwrite','Overwriting %s.',outmov);
  delete(outmov);
end

didopen = false;
for tryi = 1:5,
  try
    [vr1,nframes1,fid1,headerinfo1] = get_readframe_fcn(mov1);
    im1 = vr1(1);
    [nr1,nc1,~] = size(im1);
    [vr2,nframes2,fid2,headerinfo2] = get_readframe_fcn(mov2);
    im2 = vr2(1);
    [nr2,nc2,~] = size(im2);
 
%     vr1 = VideoReader(mov1);
%     vr2 = VideoReader(mov2);
  catch ME,
    warning('Failed to open input movies on try %d: %s',tryi,getReport(ME));
    clear vr1 vr2;
    if tryi < 5,
      fprintf('Pausing 5 seconds before trying again!\n');
      pause(5);
    end
    continue;
  end
  didopen = true;
  break;
end

if ~didopen,
  error('Could not open input movie files %s and %s for reading',mov1,mov2);
end

Nfrm = min(nframes1,nframes2);
off1 = 0;
off2 = 0;
%assert(vr1.NumberOfFrames==vr2.NumberOfFrames,'Number of frames do not match.');
assert(nr1==nr2,'Image heights do not match.');

DISP_MOD = 500;

didopen = false;
for tryi = 1:5,
  try
    vw = VideoWriter(outmov); %#ok<TNMLP>
    vw.open();
  catch ME,
    warning('Failed to open movie %s on try %d: %s',outmov,tryi,getReport(ME));
    try
      clear vw;
      if exist(outmov,'file'),
        fprintf('Deleting failed movie file %s.\n',outmov);
        delete(outmov);
      end
    catch ME2,
      warning('Could not delete failed movie file %s: %s',outmov,getReport(ME2));
    end
    if tryi < 5,
      fprintf('Pausing 5 seconds before trying again!\n');
      pause(5);
    end
    continue;
  end
  didopen = true;
  break;
end

if ~didopen,
  error('Could not open file %s',outmov);
end

try
  for i = 1:Nfrm
    im1 = vr1(i+off1);
    im2 = vr2(i+off2);
    im = [im1 im2];
    vw.writeVideo(im);
    
    if mod(i,DISP_MOD)==0 || i==Nfrm
      fprintf(1,'Wrote concatenated frame %d.\n',i);
    end
  end
  
  vw.close();  
catch ME
  error('catmov:write','Error while concatenating movies: %s\n',ME.message);  
end