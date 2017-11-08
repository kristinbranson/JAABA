function [hogftrs flowftrs framendx framehdndx] = genFeatures(dirstruct,viddet)
%% Temporary code to add more interval windows.

%% Parameters

wsz = 5; % Window size.
nwin = 4; % Number of time windows.
intsz = wsz*nwin;


npatches = 4;
psize = 40; % patch size for hog/hof
nbins = 8; % number of bins in hog/hof

readsz = 120;

debug = false;

%% Video Details

tracks = getTracks(viddet,dirstruct);
vidstats = load(viddet.vidStatsFile);
vidstats = vidstats.vidStats;
hdframes = vidstats.supplement_frame_reference{viddet.vidName};


%%
try
  [readframefcn,~,fid,headerinfo] = get_readframe_fcn(viddet.hdvid);
catch exception,
  warning('Could not open %s for computing the features %s\n',vidPath,exception.message); %#ok<*WNTAG>
  return;
end

%% If number of hdframes differs from headerinfo.nframes by 1, 
% then still continue.. its a bug with the codecs on mayank's linux 
% workstation.

if (numel(hdframes)-headerinfo.nframes) == 1,
  hdframes(end) = [];
end


%% Find intervals of hd

supphd = [-1; hdframes];
dd = supphd(2:end)-supphd(1:end-1);
int_start = find(dd>1);
int_end = int_start(2:end)-1;
int_end = [int_end;numel(hdframes)];
int_size = int_end-int_start+1;
int_remove = int_size<(2*intsz+1); 
% remove intervals smaller than the window size.
int_start(int_remove) = [];
int_end(int_remove) = [];

if isempty(int_start), 
  warning('Not big enough intervals in video %s for computing features\n',...
    viddet.vidPath);
    hogftrs = []; flowftrs = []; framendx = []; framehdndx = [];
  return;
end


%% Count the number of frames on which we will compute the features..
frcount = 0;
for fno = 1:numel(hdframes)
  if fno-intsz< 1, continue; end
  if fno+intsz> numel(hdframes), continue; end
  
  jj = hdframes( (fno-intsz):(fno+intsz));
  if any( (jj(2:end)-jj(1:end-1))>1), continue; end
  
  frcount = frcount+1;
end

% And the size of the features
patchsz = 2*npatches*psize;
fsize = size(gradientHist(single(zeros(patchsz)),single(zeros(patchsz)),psize,nbins,1));

%%

% cache to store the frames.
imMat = zeros(headerinfo.Height,headerinfo.Width,readsz);
imMatHDidx = zeros(1,readsz);
imMatActualidx = zeros(1,readsz);
gradMat = struct('M',zeros(size(imMat)),'O',zeros(size(imMat)));
flowMat = struct('M',zeros(size(imMat)),'O',zeros(size(imMat)));
flowftrs = zeros(fsize(1),fsize(2),fsize(3),2,2,frcount);
hogftrs = zeros(fsize(1),fsize(2),fsize(3),2,2,frcount);
framendx = zeros(1,frcount);
framehdndx = zeros(1,frcount);

count = 1;
if debug,
  timestart = tic; %#ok<*UNRCH>
end

for intno = 1:numel(int_start)
  
  
  tsFrame = int_start(intno):int_end(intno);
  
  imMat(:) = 0;
  imMatHDidx(:) = 0;
  imMatActualidx(:) = 0;
  
  if any(tsFrame>headerinfo.nframes),
    warning('Tracking has more frames than HD supplement video for %s (%s)\n',viddet.hdvid,viddet.vidPath);
    hogftrs = []; flowftrs = []; framendx = []; framehdndx = [];
    return;
  end
  

  %% Loop over the ts that have sufficient frames around it.
  for t = (intsz+1):(numel(tsFrame)-intsz)

    
    %% Fill up the cache..
    if ~ismember( tsFrame(t+intsz),imMatHDidx)
      imMat(:) = 0;
      imMatHDidx(:) = 0;
      imMatActualidx(:) = 0;
      gradMat.M(:) = 0; gradMat.O(:) = 0;
      flowMat.M(:) = 0; flowMat.O(:) = 0;
      
      startframe = max(tsFrame(1),tsFrame(t-intsz));
      endframe = min(tsFrame(t-intsz)+readsz-1,tsFrame(end));
      

      tt = readframesinrange(readframefcn,startframe,endframe);
      if size(tt,4)>readsz,
        fprintf('Read more than %d frames\n',readsz);
      end
      sizet = size(tt,4);
      framesreadHD = startframe:endframe;
      
      for ndx = 1:sizet
        imMat(:,:,ndx) = single(rgb2gray(tt(:,:,:,ndx)));
      end
      
      imMatP = imMat(:,:,2:end);
      GM = gradMat.M; GO = gradMat.O;
      FM = flowMat.M; FO = flowMat.O;
      
      for ndx = 1:sizet
        
        imMatHDidx(ndx) = framesreadHD(ndx);
        [M,O] = gradientMag(single(imMat(:,:,ndx))/255); 
        GM(:,:,ndx) = M;
        GO(:,:,ndx) = min(O,pi-1e-6);
        
        if ndx<sizet
          if ndx>size(imMatP,3),
            fprintf('ouuouo\n');
          end
          I1 = imMat(:,:,ndx);
          I2 = imMatP(:,:,ndx);
          [Vx,Vy,~] = optFlowLk(I1,I2,[],3);
          M = sqrt(Vx.^2 + Vy.^2);
          O = mod(atan2(Vy,Vx)/2,pi);
          O = min(O,pi-1e-6);
          FM(:,:,ndx) = M;
          FO(:,:,ndx) = O;
        end

      end
      
      gradMat.M = GM; gradMat.O = GO;
      flowMat.M = FM; flowMat.O = FO;
      
    end % Done with cache part
    
    %% Adjust for the tracking..
    
    matndx = find(tsFrame(t)==imMatHDidx);
    trackNdx = hdframes(tsFrame(t));

    if size(tracks.topTracks,1)<trackNdx,
      break;
    end
    
    locs = {round(tracks.topTracks(trackNdx,:)),...
            round(tracks.botTracks(trackNdx,:))};
          
    if any(isnan(locs{1})) || any(isnan(locs{2})),
      continue;
    end
    
    patchsz = 2*npatches*psize;
    gradFrames = (matndx-intsz):(matndx+intsz-1);
    flowFrames = (matndx-intsz):(matndx+intsz-1);
    
    curGrad.M = extractPatches(gradMat.M(:,:,gradFrames),locs,patchsz,0);
    curGrad.O = extractPatches(gradMat.O(:,:,gradFrames),locs,patchsz,pi/2);
    curFlow.M = extractPatches(flowMat.M(:,:,flowFrames),locs,patchsz,0);
    curFlow.O = extractPatches(flowMat.O(:,:,flowFrames),locs,patchsz,0);

    
    %% Feature Computation

    Hpatches = {}; Fpatches = {};
    for pno = 1:numel(curGrad.M)
      for ndx = 1:size(curGrad.M{1},3)
        Hpatches{pno}(:,:,:,ndx) = gradientHist(...
          single(curGrad.M{pno}(:,:,ndx)),...
          single(curGrad.O{pno}(:,:,ndx)),...
          psize,nbins,1); %#ok<AGROW>
      end
      
      for ndx = 1:size(curFlow.M{1},3)
        Fpatches{pno}(:,:,:,ndx) = gradientHist(...
          single(curFlow.M{pno}(:,:,ndx)),...
          single(curFlow.O{pno}(:,:,ndx)),...
          psize,nbins,1); %#ok<AGROW>
      end
    end
    
    % if current frame is t, then 
    % hogftrs(:,:,:,:,1,:) is grad hist of t-wsz+1,t-wsz+1,...,t 
    %  Before t and including t.
    % hogftrs(:,:,:,:,2,:) is grad hist of t,t+1,...,t+wsz-1   
    %  After t and including t.
    % flowftrs(:,:,:,:,1,:) is grad hist of pairs (t-wsz,t-wsz+1),...,(t-1,t)  
    %  This is strictly before t.
    % flowftrs(:,:,:,:,2,:) is grad hist of pairs (t,t+1),...,(t+wsz-1,t+wsz)      
    %  This is strictly after t.
    
    % Indexing into hogftrs and flowftrs..
    % 1st dim: Y loc of the hog cell
    % 2nd dim: X loc of the hog cell
    % 3rd dim: gradient bin
    % 4th dim: Interest point or patch no around which hog/hof is computed.
    % 5th dim: Before or after the current frame.
    % 6th dim: Frame number.
    
    for idx = 1:2*nwin
      
      curint = ((idx-1)*wsz) + (1:wsz);
      
      for pno = 1:numel(curGrad.M)
        hogftrs(:,:,:,pno,idx,count) = sum(Hpatches{pno}(:,:,:,curint),4);
        flowftrs(:,:,:,pno,idx,count) = sum(Fpatches{pno}(:,:,:,curint),4);
      end
      
    end
    framendx(count) = hdframes(tsFrame(t));
    framehdndx(count) = tsFrame(t);
    
    count = count+1;
    
  end
  
end

if debug,
  fprintf('Time without debug');
  toc(timestart);
end

%%

if debug,
  tic;
  count = 1;
  imMat = zeros(headerinfo.Height,headerinfo.Width,intsz);
  flowftrsd = zeros(fsize(1),fsize(2),fsize(3),2,2,frcount);
  hogftrsd = zeros(fsize(1),fsize(2),fsize(3),2,2,frcount);
  framendxd = zeros(1,frcount);
  framehdndxd = zeros(1,frcount);

  ii = readframefcn([1 headerinfo.nframes]);
  for fno = 1:numel(hdframes)
    fprintf('.');
    if mod(fno,30)==0, fprintf('\n'); end
    if fno-intsz< 1, continue; end
    if fno+intsz> numel(hdframes), continue; end
    
    jj = hdframes( (fno-intsz):(fno+intsz));
    if any( (jj(2:end)-jj(1:end-1))>1), continue; end
    
    trackNdx = hdframes(fno);
    locs = {round(tracks.topTracks(trackNdx,:)),...
            round(tracks.botTracks(trackNdx,:))};
    patchsz = 2*npatches*psize;

    % Before t features.
    gradMat = struct('M',zeros(size(imMat)),'O',zeros(size(imMat)));
    flowMat = struct('M',zeros(size(imMat)),'O',zeros(size(imMat)));

    for ndx = (fno-intsz):(fno-1)
      I = single(rgb2gray(ii(:,:,:,ndx)));
      I1 = single(rgb2gray(ii(:,:,:,ndx+1)));
      [Vx,Vy,~] = optFlowLk(I,I1,[],3);
      M = sqrt(Vx.^2 + Vy.^2);
      O = mod(atan2(Vy,Vx)/2,pi);
      flowMat.M(:,:,ndx-fno+intsz+1) = M;
      flowMat.O(:,:,ndx-fno+intsz+1) = min(O,pi-1e-6);
      
      [M,O] = gradientMag(single(I)/255);
      gradMat.M(:,:,ndx-fno+intsz+1) = M;
      gradMat.O(:,:,ndx-fno+intsz+1) = min(O,pi-1e-6);

    end
    

    curGrad.M = extractPatches(gradMat.M,locs,patchsz,0);
    curGrad.O = extractPatches(gradMat.O,locs,patchsz,pi/2);
    curFlow.M = extractPatches(flowMat.M,locs,patchsz,0);
    curFlow.O = extractPatches(flowMat.O,locs,patchsz,0);

    Hpatches = {}; Fpatches = {};
    for pno = 1:numel(curGrad.M)
      for ndx = 1:size(curGrad.M{1},3)
        Hpatches{pno}(:,:,:,ndx) = gradientHist(...
          single(curGrad.M{pno}(:,:,ndx)),...
          single(curGrad.O{pno}(:,:,ndx)),...
          psize,nbins,1); 
      end
      
      for ndx = 1:size(curFlow.M{1},3)
        Fpatches{pno}(:,:,:,ndx) = gradientHist(...
          single(curFlow.M{pno}(:,:,ndx)),...
          single(curFlow.O{pno}(:,:,ndx)),...
          psize,nbins,1); 
      end
    end

    for idx = 1:nwin
      curint = ((idx-1)*wsz) + (1:wsz);

      for pno = 1:numel(curGrad.M)
        hogftrsd(:,:,:,pno,idx,count) = sum(Hpatches{pno}(:,:,:,curint),4);
        flowftrsd(:,:,:,pno,idx,count) = sum(Fpatches{pno}(:,:,:,curint),4);
      end
    end
    
    
    
    
    % After t features.
    gradMat = struct('M',zeros(size(imMat)),'O',zeros(size(imMat)));
    flowMat = struct('M',zeros(size(imMat)),'O',zeros(size(imMat)));

    for ndx = fno:(fno+intsz-1)
      I = single(rgb2gray(ii(:,:,:,ndx)));
      I1 = single(rgb2gray(ii(:,:,:,ndx+1)));
      [Vx,Vy,~] = optFlowLk(I,I1,[],3);
      M = sqrt(Vx.^2 + Vy.^2);
      O = mod(atan2(Vy,Vx)/2,pi);
      flowMat.M(:,:,ndx-fno+1) = M;
      flowMat.O(:,:,ndx-fno+1) = min(O,pi-1e-6);
      
      [M,O] = gradientMag(single(I)/255);
      gradMat.M(:,:,ndx-fno+1) = M;
      gradMat.O(:,:,ndx-fno+1) = min(O,pi-1e-6);

    end
    

    curGrad.M = extractPatches(gradMat.M,locs,patchsz,0);
    curGrad.O = extractPatches(gradMat.O,locs,patchsz,pi/2);
    curFlow.M = extractPatches(flowMat.M,locs,patchsz,0);
    curFlow.O = extractPatches(flowMat.O,locs,patchsz,0);

    Hpatches = {}; Fpatches = {};
    for pno = 1:numel(curGrad.M)
      for ndx = 1:size(curGrad.M{1},3)
        Hpatches{pno}(:,:,:,ndx) = gradientHist(...
          single(curGrad.M{pno}(:,:,ndx)),...
          single(curGrad.O{pno}(:,:,ndx)),...
          psize,nbins,1); 
      end
      
      for ndx = 1:size(curFlow.M{1},3)
        Fpatches{pno}(:,:,:,ndx) = gradientHist(...
          single(curFlow.M{pno}(:,:,ndx)),...
          single(curFlow.O{pno}(:,:,ndx)),...
          psize,nbins,1); 
      end
    end
    
    for idx = 1:nwin
      curint = ((idx-1)*wsz) + (1:wsz);
      
      for pno = 1:numel(curGrad.M)
        hogftrsd(:,:,:,pno,idx+nwin,count) = sum(Hpatches{pno}(:,:,:,curint),4);
        flowftrsd(:,:,:,pno,idx+nwin,count) = sum(Fpatches{pno}(:,:,:,curint),4);
      end
    end


    framendxd(count) = hdframes(fno);
    framehdndxd(count) = fno;
    
    count = count+1;
  end
  
  fprintf('Time with debug\n');
  toc;

  if ~isequal(hogftrs,hogftrsd)  || ...
    ~isequal(flowftrs,flowftrsd) || ...
    ~isequal(framendx,framendxd) || ...
    ~isequal(framehdndx,framehdndxd)
    fprintf('Debug failed. Something doesn'' match\n');
  end

end


if fid>0,
  fclose(fid);
end



  
function tt = readframesinrange(readframe,ts2readstart,ts2readend)
s = warning('error','MATLAB:audiovideo:VideoReader:incompleteRead'); %#ok<CTPCT>
try
  tt = readframe([ts2readstart ts2readend]);
catch ME,
  if strcmp(ME.identifier,'MATLAB:audiovideo:VideoReader:incompleteRead')
    fprintf('Reading Frames individually\n');
    tt = [];
    for ndx = ts2readstart:ts2readend
      tt(:,:,:,end+1) = readframe(ndx); %#ok<AGROW>
    end
  else
    rethrow(ME);
  end
end
fprintf('Read frames:%d to %d\n',ts2readstart,ts2readend);
warning(s);


function patches = extractPatches(img,locs,psz,padval)

patches = cell(1,numel(locs));
for ndx = 1:numel(locs)
  patches{ndx} = padgrab(img,padval,locs{ndx}(2)-psz/2,locs{ndx}(2)+psz/2-1,...
                             locs{ndx}(1)-psz/2,locs{ndx}(1)+psz/2-1,...
                             1,size(img,3));
  
end
