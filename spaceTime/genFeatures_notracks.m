function ftrs = genFeatures_notracks(readfcn,headerinfo,fstart,fend,method,params)

%% parameters

% params = getParams;
npatchesy = params.npatchesy;
npatchesx = params.npatchesx;
psize = params.psize;
nbins = params.nbins; 

optflowwinsig = params.optflowwinsig ;
optflowsig = params.optflowsig ;
optreliability = params.optreliability ;


%% Read the images into a buffer.

tic;

nframes = numel(fstart:fend+1);
im = uint8(zeros(headerinfo.nr,headerinfo.nc,nframes));
count = 1;
for ndx = fstart:fend
  im(:,:,count) = readfcn(ndx);
  count = count + 1;
end
if fend==headerinfo.nframes
  im(:,:,count) = im(:,:,count-1);
else
  im(:,:,count) = readfcn(fend+1);
end


%% preallocate hogftrs and flowftrs

nflies = 1;
ftrsz = [npatchesy npatchesx nbins]; 

hogftrs = cell(1,1);
flowftrs = cell(1,1);
curframes =fstart:fend;
frames_fly = numel(curframes);
hogftrs{1} = zeros([ftrsz frames_fly],'single');
flowftrs{1} = zeros([ftrsz frames_fly],'single');


%% Compute the features.

flycount = ones(1,nflies);
if strcmp(method,'LK') % OLD 
  methodC = 1;
elseif strcmp(method,'hs-brightness')
  methodC = 2;
elseif strcmp(method,'deep-sup')
  methodC = 3;
else
  error('Unknown Method')
end

for ndx = fstart:fend
  
  if params.presmoothing
    curpatch = gaussSmooth(im(:,:, ndx-fstart + 1),params.smoothing_sigma,'same');
    curpatch2 = gaussSmooth(im(:,:, ndx-fstart + 2),params.smoothing_sigma,'same');
  else
    curpatch = im(:,:, ndx-fstart + 1);
    curpatch2 = im(:,:, ndx-fstart + 2);
  end
  
  curpatch = single(curpatch);
  curpatch2 = single(curpatch2);
  [M,O] = gradientMag(curpatch(:,:,1)/255);
  O = min(O,pi-1e-6);
  hogftrs{1}(:,:,:,ndx-fstart+1) = single(gradientHist(M,O,psize,nbins,false));
  if methodC ==1
    [Vx,Vy,~] = optFlowLk(curpatch,curpatch2,[],optflowwinsig,optflowsig,optreliability); % Using soft windows
  elseif methodC ==2
    curpatch = uint8(curpatch);
    curpatch2 = uint8(curpatch2);
    uv = estimate_flow_interface(curpatch,curpatch2,method,...
      {'max_warping_iters',params.warping_iters});
    Vx = uv(:,:,1);
    Vy = uv(:,:,2);
  else
      error('Undefined method');
  end
  M = single(sqrt(Vx.^2 + Vy.^2));
  O = mod(atan2(Vy,Vx)/2,pi);
  O = single(min(O,pi-1e-6));
  flowftrs{1}(:,:,:,ndx-fstart+1) = single(gradientHist(M,O,psize,nbins,true));
    
end
  
ftrs.hogftrs = hogftrs;
ftrs.flowftrs = flowftrs;
  
