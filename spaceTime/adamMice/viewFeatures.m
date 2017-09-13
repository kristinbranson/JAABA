function viewFeatures(jabfile,movname)

V = load(jabfile,'-mat');

Fpos = [];
Hpos = [];
Fneg = [];
Hneg = [];
imMatpos = [];
imMatneg = [];
patchMatpos = [];
patchMatneg = [];
width = 640;
height = 480;

if nargin<2
  savemov = false;
else
  savemov = true;
end
  

for expi = 1:numel(V.x.expDirNames)
  vidfile = fullfile(V.x.expDirNames{expi},V.x.file.moviefilename);
  [readframe,nframes,fid,headerinfo] = get_readframe_fcn(vidfile);
  
  readFrameFcns.readframe = readframe;
  readFrameFcns.nframes = nframes;
  readFrameFcns.headerinfo = headerinfo;

  for fly = 1:numel(V.x.labels(expi).t0s)
    for bnum = 1:numel(V.x.labels(expi).names{fly})
      t0 = V.x.labels(expi).t0s{fly}(bnum);
      t1 = V.x.labels(expi).t1s{fly}(bnum);

      curimmat = zeros(headerinfo.Height,headerinfo.Width,t1-t0);
      count = 1;
      for ndx = t0:t1-1
        curimmat(:,:,count) = rgb2gray(readframe(ndx));
        count = count+1;
      end
      
      if size(curimmat,1)>height
        curimmat(height+1:end,:,:) = [];
      elseif size(curimmat,1)<height
        curimmat(end+1:height,:,:) = 0;
      end
      
      if size(curimmat,2)>width
        curimmat(:,width+1:end,:) = [];
      elseif size(curimmat,2)<width
        curimmat(:,end+1:width,:) = 0;
      end

      
      [F,H,I] = genFeatures('ts',t0:t1-1,'readframeFcns',readFrameFcns,...
        'trxfile',fullfile(V.x.expDirNames{expi},V.x.file.trxfilename));

      if strcmpi(V.x.labels(expi).names{fly}{bnum},'none')
        Fneg = cat(4,Fneg,F);
        Hneg = cat(4,Hneg,H);
        imMatneg = cat(3,imMatneg,curimmat);
        patchMatneg = cat(3,patchMatneg,I);
      else
        Fpos = cat(4,Fpos,F);
        Hpos = cat(4,Hpos,H);
        imMatpos = cat(3,imMatpos,curimmat);
        patchMatpos = cat(3,patchMatpos,I);
      end
      
    end
    
    
  end
  
  if ~isempty(vidfile) && fid>0,
    fclose(fid);
  end

  
end
Fposframe(2) = struct('cdata',[],'colormap',[]);
Fnegframe(2) = struct('cdata',[],'colormap',[]);

prm = struct('nn',2);
fclims = prctile(Fneg(:),[0.5 99.5]);
hclims = prctile(Hneg(:),[0.5 99.5]);

figpos = figure('units','normalized','outerposition',[0 0 1 1]);
for ndx = 1:size(Fpos,4)
  figure(figpos);
  curH = Hpos(:,:,:,ndx);
  subplot(2,3,[1 4]); montage2(Fpos(:,:,:,ndx),prm);
  set(gca,'CLim',fclims);
  subplot(2,3,[2 5]) ; montage2(curH,prm);
  set(gca,'CLim',hclims);
  subplot(2,3,3);
  imshow(imMatpos(:,:,ndx)/255);
  subplot(2,3,6);
  imshow(patchMatpos(:,:,ndx)/255);
  Fposframe(ndx) = getframe;
  pause;
end
figneg = figure('units','normalized','outerposition',[0 0 1 1]);
for ndx = 1:size(Fneg,4)
  figure(figneg);
  curH = Hneg(:,:,:,ndx);
  subplot(2,3,[1 4]); montage2(Fneg(:,:,:,ndx),prm);
  set(gca,'CLim',fclims);
  subplot(2,3,[2 5]); montage2(curH,prm);
  set(gca,'CLim',hclims);
  subplot(2,3,3);
  imshow(imMatneg(:,:,ndx)/255);
  subplot(2,3,6);
  imshow(patchMatneg(:,:,ndx)/255);
  Fnegframe(ndx) = getframe;
  pause;
end

if savemov,
  
  writerObj = VideoWriter([movname 'pos.avi']);
  writerObj.FrameRate = 4;
  open(writerObj);
  for ndx = 1:numel(Fposframe)
    writeVideo(writerObj,Fposframe(ndx));
  end
  close(writerObj);
  
  
  writerObj = VideoWriter([movname 'neg.avi']);
  writerObj.FrameRate = 4;
  open(writerObj);
  for ndx = 1:numel(Fnegframe)
    writeVideo(writerObj,Fnegframe(ndx));
  end
  close(writerObj);
  
  
end