function [flowftrs,hogftrs,info] = genFeatures(varargin)

%%

wsize = 5;
psize = 40;
nbins = 8;
halfpsize = round(psize*10/2);
hsz = round( (wsize - 1)/2);
debug = false;

% * Video: either readframeFcns or vidfile should be specified.
% * Interest pts: If interest points are desired, specify either ips or
% trxfile. If interest points are not desired (side only), leave both ips
% and trxfile empty.
[ts,readframeFcns,vidfile,trxfile,ips,face,framechunksize,frontside,...
  dopatchprint] = myparse(varargin,...
  'ts',[],...
  'readframeFcns',[],... 
  'vidfile','',...
  'trxfile','',...
  'ips',[],...
  'face','',...
  'framechunksize',[],...
  'frontside',false,...
  'dopatchprint',true);

if isempty(framechunksize),
  
  framechunksize = 100;

end

if isempty(face)
  error('Face direction of the mice is not specified');
end

if strcmpi(face,'left')
  do_fliplr = true;
else
  do_fliplr = false;
end

if isempty(readframeFcns) && isempty(vidfile)
  error('Both arguments are empty');
end

tfvidfile = false;
if isempty(readframeFcns)
  if frontside
    assert(isstruct(vidfile) && all(isfield(vidfile,{'front' 'side'})),...
      'vidfile must be a struct with fields .front, .side.');
    [readFrameFcns.front.readframe,readFrameFcns.front.nframes,readFrameFcns.front.fid,readFrameFcns.front.headerinfo] = get_readframe_fcn(vidfile.front);
    [readFrameFcns.side.readframe,readFrameFcns.side.nframes,readFrameFcns.side.fid,readFrameFcns.side.headerinfo] = get_readframe_fcn(vidfile.side);       
  else
    [readFrameFcns.side.readframe,readFrameFcns.side.nframes,readFrameFcns.side.fid,readFrameFcns.side.headerinfo] = get_readframe_fcn(vidfile);
  end
  tfvidfile = true;
else
  if ~frontside
    readFrameFcns.side = readframeFcns;
  end
end

% offside = 0;
% offfront = 0;
if frontside
  if readFrameFcns.front.nframes~=readFrameFcns.side.nframes
    warning('genFeatures:nframe',...
      'Number of frames in front video (%d) and side video (%d) do not match. Using minimum.',...
      readFrameFcns.front.nframes,readFrameFcns.side.nframes);
  end
  nframes = min(readFrameFcns.front.nframes,readFrameFcns.side.nframes);  
  % turns out we should cut from the end
  %offfront = readFrameFcns.front.nframes-nframes;
  %offside = readFrameFcns.side.nframes-nframes;
else
  nframes = readFrameFcns.side.nframes;
end

if isempty(ts)
  ts = 1:nframes;
  if any(diff(ts)~=1),
    error('ts must be starframe:endframe');
  end   
end

if ~isempty(ips) && ~isempty(trxfile)
  error('Both ips and trxfile specified.');  
elseif ~isempty(ips) && isempty(trxfile)
  if frontside
    assert(all(isfield(ips,{'front' 'side'})),'Expected fields missing from ips.');
  elseif ~isstruct(ips)
    ips = struct('side',ips);
  end
elseif isempty(ips) && ~isempty(trxfile)
  T = load(trxfile);
  ips.side(1,:) = T.trx(1).arena.food;
  ips.side(2,:) = T.trx(1).arena.mouth;
  ips.side(3,:) = T.trx(1).arena.perch;
  ips.side = round(ips.side);
  if frontside
    ips.front(1,:) = T.trx(1).arena.foodfront;
    ips.front(2,:) = T.trx(1).arena.mouthfront;
    ips.front = round(ips.front);
  end
  assert(strcmp(face,T.trx(1).arena.face),'Specified ''face'' does not match ''face'' in specified trxfile');
else
  % isempty(ips) && isempty(trxfile)
  % no-op; handled below
end

if do_fliplr
  imw = readframeFcns.side.headerinfo.Width;
  ips.side(:,1) = imw - ips.side(:,1);  
  if frontside
    imw = readframeFcns.front.headerinfo.Width;
    ips.front(:,1) = imw - ips.front(:,1);    
  end
end

if isempty(ips),
  assert(~frontside,'Frontside currently requires interest point specification');
  patchsz = [headerinfo.Height,headerinfo.Width];
else
  Npatch = size(ips.side,1);
  if frontside
    Npatch = Npatch + size(ips.front,1);
  end  
  patchsz1 = halfpsize*2+1;
  patchsz = [patchsz1,patchsz1*Npatch];
end

% initialize
flowftrs = [];
hogftrs = [];

% buffer of the past wsize frames
Fall = {};
Hall = {};
buffert0 = -1;
buffert1 = -2;
first = true;
% Fallslow = {};
% Hallslow = {};

% previous and current patch
patchPrev = zeros(patchsz,'single');
patchCurr = zeros(patchsz,'single');

% which frames to read
ts2readstart = max(ts(1)-hsz,1);
ts2readend = min(ts(end)+hsz,nframes);

s = warning('error','MATLAB:audiovideo:VideoReader:incompleteRead');

% initialize buffer

if dopatchprint && tfvidfile && ~isdeployed && numel(ts) >= 2+hsz,
  randomframe = randsample(ts(2+hsz:end),1);
else
  randomframe = nan;
end

tic;
for ndx = 1:numel(ts),

  t = ts(ndx);

  % window considered during this frame
  wstart = max(1,t-hsz);
  wend = min(nframes,t+hsz);
  curn = wend - wstart +1;  

  % reuse the previously computed features.
  if numel(Hall)>=curn,
    Fall(1:end-1) = Fall(2:end);
    Hall(1:end-1) = Hall(2:end);
    buffert0 = buffert0+1;
  end

  % initialize HOG, HOF buffer
  if first,
       
    for wndx = 1:curn,
      
      patchPrev = patchCurr;
      patchCurr = lclGrabPatch(readFrameFcns,wstart+wndx-1,do_fliplr,tfvidfile,frontside,ips,halfpsize,patchCurr);
      
      % HOG
      % step size is two for gradient computation
      % approx equiv Matlab code, ignoring borders
      % gx = patchCurr(2:end-1,3:end)-patchCurr(2:end-1,1:end-2); 
      % gy = patchCurr(3:end,2:end-1)-patchCurr(1:end-2,2:end-1);
      % M = sqrt(gx.^2+gy.^2)/2;
      % O = modrange(atan2(gy,gx),0,pi);
      % will match M(2:end-1,2:end-1) and O(2:end-1,2:end-1)
      [M,O] = gradientMag(patchCurr/255); 
      H=gradientHist(M,O,psize,nbins,1);
      Hall{wndx} = H;

      if wndx > 1,
      
        % HOF
        [Vx,Vy,~] = optFlowLk(patchPrev,patchCurr,[],3);
        M = sqrt(Vx.^2 + Vy.^2);
        O = mod(atan2(Vy,Vx)/2,pi);
        O = min(O,pi-1e-6);
        H = gradientHist(single(M),single(O),psize,nbins,1);
        Fall{wndx-1} = H;
        
      end

    end
    buffert0 = wstart;
    buffert1 = wend;

      
  else

    % add to the buffer
    if buffert1 < wend,

      patchPrev = patchCurr;
      patchCurr = lclGrabPatch(readFrameFcns,wend,do_fliplr,tfvidfile,frontside,ips,halfpsize,patchCurr);
      
      % HOG feature for frame wend
      [M,O] = gradientMag(patchCurr/255);
      H=gradientHist(M,O,psize,nbins,1);
      Hall{curn} = H; %#ok<AGROW>
    
      [Vx,Vy,~] = optFlowLk(patchPrev,patchCurr,[],3);
      % Gradient hist requires orientation to be between 0 and pi. Makes sense for
      % image gradients. But for flow pi/2 is different than -pi/2. So adjust for
      % that.
      
      M = sqrt(Vx.^2 + Vy.^2);
      O = mod(atan2(Vy,Vx)/2,pi);
      O = min(O,pi-1e-6);
      H = gradientHist(single(M),single(O),psize,nbins,1);
      Fall{curn-1} = H; %#ok<AGROW>
      buffert1 = buffert1 + 1;
    end
    
  end
  
  curFall = Fall{1};
  curHall = Hall{1};
  for ii = 2:curn-1
    curFall = curFall + Fall{ii};
    curHall = curHall + Hall{ii};
  end
  curHall = curHall + Hall{curn};
  
  
  curFall = curFall/(curn-1);
  curHall = curHall/curn;

  if first,
    flowftrs = zeros([size(curFall),numel(ts)],'single');
    hogftrs = zeros([size(curFall),numel(ts)],'single');
  end
  
  flowftrs(:,:,:,ndx) = curFall;
  hogftrs(:,:,:,ndx) = curHall;

  if mod(ndx,10)==0,
    fprintf('Processed frame %d / %d\n',ndx,numel(ts));
  end
  first = false;  
  
  if buffert1 == randomframe,
    hFig = figure;
    imshow(patchCurr,[0 256]);
    fname = sprintf('patch_%d.jpg',randomframe);
    if ischar(vidfile)
      pth = fileparts(vidfile);
    else
      pth = fileparts(vidfile.side);
    end
    fname = fullfile(pth,fname);
    print(hFig,'-djpeg',fname);
    delete(hFig);
  end    
  
end
warning(s);

fprintf('\n');
fprintf('Generating features for %d frames took %f seconds.\n',numel(ts),toc);

info = struct();
info.nframes = nframes;

function imMat = lclReadFrame(rff,i0,i1,do_fliplr,tfvidfile)
try
  tt = rff.readframe([i0 i1]);
catch ME,
  if strcmp(ME.identifier,'MATLAB:audiovideo:VideoReader:incompleteRead')
    fprintf('Reading Frames individually\n');
    tt = [];
    for ndx = i0:i1
      tt(:,:,:,end+1) = rff.readframe(ndx); %#ok<AGROW>
    end
  else
    rethrow(ME);
  end
end

iscolor = size(tt,3) > 1;
  
nfrm = i1-i0+1;
imMat = zeros(rff.headerinfo.Height,rff.headerinfo.Width,nfrm,'single');
  
for ndx = 1:nfrm,
  if iscolor,
    ttgraycurr = rgb2gray(tt(:,:,:,ndx));
  else
    ttgraycurr = tt(:,:,:,ndx);
  end
  if do_fliplr,
    qq = fliplr(ttgraycurr);
  else
    qq = ttgraycurr;
  end
  imMat(:,:,ndx) = single(qq);
end

if tfvidfile && rff.fid>0
  fclose(rff.fid);
end

function patchMat = lclGrabPatch(readFrameFcns,f,do_fliplr,tfvidfile,frontside,ips,halfpsize,patchMat)

imMat.side = lclReadFrame(readFrameFcns.side,f,f,do_fliplr,tfvidfile);
if frontside
  imMat.front = lclReadFrame(readFrameFcns.front,f,f,do_fliplr,tfvidfile);
end

patchsz1 = size(patchMat,1);
if ~isempty(ips)
  Nipsside = size(ips.side,1);
  for pndx = 1:Nipsside
    pp = padgrab(imMat.side,0,...
      ips.side(pndx,2)-halfpsize,ips.side(pndx,2)+halfpsize,...
      ips.side(pndx,1)-halfpsize,ips.side(pndx,1)+halfpsize,...
      1,size(imMat.side,3));
    patchMat(:,(pndx-1)*patchsz1+1:pndx*patchsz1) = pp;
  end
  if frontside
    for pndx = 1:size(ips.front,1)
      pp = padgrab(imMat.front,0,...
        ips.front(pndx,2)-halfpsize,ips.front(pndx,2)+halfpsize,...
        ips.front(pndx,1)-halfpsize,ips.front(pndx,1)+halfpsize,...
        1,size(imMat.front,3));
      patchMat(:,(pndx-1+Nipsside)*patchsz1+1:(pndx+Nipsside)*patchsz1) = pp;
    end
  end
else
  assert(~frontside);
  patchMat = imMat.side;
end