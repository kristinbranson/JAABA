function ftrs = genSTFeatures(readfcn,headerinfo,fstart,fend,tracks,stationary,method,params)
% stationary = true will center the animals before
% computing the flow

if nargin< 6
    if isfield(params,'is_stationary')
        stationary=params.is_stationary;
    else
        stationary = false;
    end
end

%% parameters

% params = getParams;
npatches_x = params.npatches_x;
npatches_y = params.npatches_y;
psize = params.psize;
nbins = params.nbins; 
patchsz_x = npatches_x*psize;
patchsz_y = npatches_y*psize;

optflowwinsig = params.optflowwinsig ;
optflowsig = params.optflowsig ;
optreliability = params.optreliability ;


%% Read the images into a buffer.

tic;

nframes = numel(fstart:fend+1);
im = uint8(zeros(headerinfo.nr,headerinfo.nc,nframes));
count = 1;

fr = readfcn(fstart);
n_colors = 1;
if ndims(fr) > 2,
  if size(fr,3)==3,
    n_colors = 3;
  end
end

for ndx = fstart:fend
  if n_colors > 1,
    im(:,:,count) = rgb2gray(readfcn(ndx));
  else
    im(:,:,count) = readfcn(ndx);
  end
  count = count + 1;
end
if fend==headerinfo.nframes
  im(:,:,count) = im(:,:,count-1);
else
  if n_colors > 1,
    im(:,:,count) = rgb2gray(readfcn(fend+1));
  else
    im(:,:,count) = readfcn(fend+1);
  end
end

tread =toc;
% fprintf('Read %d frames in %.2fs\n',nframes,tread);

%% preallocate hogftrs and flowftrs

nflies = numel(tracks);
ftrsz = [npatches_y npatches_x nbins]; 

hogftrs = cell(1,nflies);
flowftrs = cell(1,nflies);
for fly = 1:nflies
  curframes =fstart:fend;
  frames_fly = nnz( curframes<=tracks(fly).endframe & ...
                    curframes>=tracks(fly).firstframe);
  hogftrs{fly} = zeros([ftrsz frames_fly],'single');
  flowftrs{fly} = zeros([ftrsz frames_fly],'single');
  
end

%% Compute the features.

flycount = ones(1,nflies);
methodC = 0;
if strcmp(method,'LK') % OLD 
  methodC = 1;
elseif strcmp(method,'hs-brightness')
  methodC = 2;
elseif strcmp(method,'hs-sup') 
  % HS with flow in background and fly suppressed
  methodC = 3;
elseif strcmp(method,'deep-sup')
  methodC = 4;
else
  error('Unknown Method')
end

% For background suppression
fr_params = struct;
% bwimg = zeros(patchsz,patchsz);
% ctr = [ceil( (patchsz+1)/2),ceil( (patchsz+1)/2)];
% bwimg(ctr(1),ctr(2))=1;
% dimg = bwdist(bwimg,'euclidean');
% [xx,yy]= meshgrid(1:patchsz,1:patchsz);
% aimg = atan2(-(yy-ctr(1)),-(xx-ctr(2)));
% params.dimg = dimg; params.aimg = aimg;

for ndx = fstart:fend
  
  for fly = 1:nflies
    if tracks(fly).firstframe > ndx || tracks(fly).endframe < ndx,
      continue;
    end
    trackndx = ndx - tracks(fly).firstframe + 1;
%     locy = double(round(tracks(fly).y(trackndx)));
%     locx = double(round(tracks(fly).x(trackndx)));
%    curpatch = extractPatch(im(:,:, ndx-fstart + 1),...
%      locy,locx,tracks(fly).theta(trackndx),patchsz);
    locy = tracks(fly).y(trackndx);
    locx = tracks(fly).x(trackndx);
    curpatch = CropImAroundTrx(im(:,:, ndx-fstart + 1),...
      locx,locy,tracks(fly).theta(trackndx)-pi/2,(patchsz_x-1)/2,(patchsz_y-1)/2);

    if stationary && ndx<tracks(fly).endframe
%       locy = double(round(tracks(fly).y(trackndx+1)));
%       locx = double(round(tracks(fly).x(trackndx+1)));
      locy = tracks(fly).y(trackndx+1);
      locx = tracks(fly).x(trackndx+1);
    end
    
    if ~stationary || ndx==tracks(fly).endframe
%       curpatch2 = extractPatch(im(:,:, ndx-fstart + 2),...
%         locy,locx,tracks(fly).theta(trackndx),patchsz);
      curpatch2 = CropImAroundTrx(im(:,:, ndx-fstart + 2),...
        locx,locy,tracks(fly).theta(trackndx)-pi/2,(patchsz_x-1)/2,(patchsz_y-1)/2);
    else
%       curpatch2 = extractPatch(im(:,:, ndx-fstart + 2),...
%         locy,locx,tracks(fly).theta(trackndx+1),patchsz);
      curpatch2 = CropImAroundTrx(im(:,:, ndx-fstart + 2),...
        locx,locy,tracks(fly).theta(trackndx+1)-pi/2,(patchsz_x-1)/2,(patchsz_y-1)/2);
    end
    
    
    curpatch = single(curpatch);
    curpatch2 = single(curpatch2);
    [M,O] = gradientMag(curpatch(:,:,1)/255);
    O = min(O,pi-1e-6);
    hogftrs{fly}(:,:,:,flycount(fly)) = single(gradientHist(M,O,psize,nbins,false));
    if methodC ==1,
      [Vx,Vy,~] = optFlowLk(curpatch,curpatch2,[],optflowwinsig,optflowsig,optreliability); % Using soft windows
    elseif methodC ==2
      curpatch = uint8(curpatch);
      curpatch2 = uint8(curpatch2);
      uv = estimate_flow_interface(curpatch,curpatch2,method,{'max_warping_iters',2});
      Vx = uv(:,:,1);
      Vy = uv(:,:,2);
    elseif methodC >=3
      curpatch = uint8(curpatch);
      curpatch2 = uint8(curpatch2);
      if ndx<tracks(fly).endframe
%         fr_params.dx = round(tracks(fly).x(trackndx+1)) - round(tracks(fly).x(trackndx));
%         fr_params.dy = round(tracks(fly).y(trackndx+1)) - round(tracks(fly).y(trackndx));
%         fr_params.dtheta = tracks(fly).theta(trackndx+1) - tracks(fly).theta(trackndx);
%         fr_params.rdx = tracks(fly).x(trackndx+1) - tracks(fly).x(trackndx)-fr_params.dx;
%         fr_params.rdy = tracks(fly).y(trackndx+1) - tracks(fly).y(trackndx)-fr_params.dy;
        fr_params.dx = tracks(fly).x(trackndx+1) - tracks(fly).x(trackndx);
        fr_params.dy = tracks(fly).y(trackndx+1) - tracks(fly).y(trackndx);
        fr_params.dtheta = tracks(fly).theta(trackndx+1) - tracks(fly).theta(trackndx);
        fr_params.rdx = tracks(fly).x(trackndx+1) - tracks(fly).x(trackndx)-fr_params.dx;
        fr_params.rdy = tracks(fly).y(trackndx+1) - tracks(fly).y(trackndx)-fr_params.dy;
      else
        fr_params.dx = 0;
        fr_params.dy = 0;
        fr_params.dtheta = 0;
        fr_params.rdx = 0;
        fr_params.rdy = 0;
      end
      fr_params.theta = tracks(fly).theta(trackndx);
      fr_params.stationary = stationary;
      if methodC==3
        [Vx,Vy] = computeFlowBkgSup(curpatch,curpatch2,fr_params,params);
      elseif methodC==4
        [Vx,Vy] = computeDeepFlowBkgSup(curpatch,curpatch2,fr_params,params);        
      else
        error('Undefined method');
      end
    end
%    [Vx,Vy,~] = optFlowLk(curpatch,curpatch2,4,[],optflowsig); % Using hard windows
    M = single(sqrt(Vx.^2 + Vy.^2));
    O = mod(atan2(Vy,Vx)/2,pi);
    O = single(min(O,pi-1e-6));
    flowftrs{fly}(:,:,:,flycount(fly)) = single(gradientHist(M,O,psize,nbins,true));
    
    
    flycount(fly) = flycount(fly) + 1;
  end
end
  
ftrs.hogftrs = hogftrs;
ftrs.flowftrs = flowftrs;
  
