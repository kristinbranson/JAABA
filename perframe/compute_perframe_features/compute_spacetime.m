% absolute azimuthal angle to from fly to closest fly according to type
function [data,units] = compute_spacetime(trx,n,theta_in,r_in,type_in)

try
  load(fullfile(trx.expdirs{n},trx.perframedir,'spacetime.mat'));
  
catch
  flies = trx.exp2flies{n};
  nflies = numel(flies);
  data = cell(1,nflies);

  % parameters

  meana = prctile(cellfun(@nanmean,trx(:).a),90);
  meanb = prctile(cellfun(@nanmean,trx(:).b),90);

  [binidx, nbins, featurenames, featureboundaries, featurecenters] = compute_spacetime_mask(meana, meanb);
  
%   featurenames{1}={featurenames{1}{:} ...
%       'theta28_r1' 'theta28_r2' 'theta28_r3' ...
%       'theta37_r1' 'theta37_r2' 'theta37_r3' ...
%       'theta46_r1' 'theta46_r2' 'theta46_r3' ...
%       'theta42_r1' 'theta42_r2' 'theta42_r3' ...
%       'theta51_r1' 'theta51_r2' 'theta51_r3' ...
%       'theta68_r1' 'theta68_r2' 'theta68_r3' ...
%       'theta1_r0' 'theta2_r0' 'theta3_r0' 'theta4_r0' ...
%       'theta5_r0' 'theta6_r0' 'theta7_r0' 'theta8_r0'};
%   
  % generate a random image
  % 
  % nr = 1024;
  % nc = 1024;
  % pos = struct;
  % pos.x = 234;
  % pos.y = 345;
  % pos.a = meana+3;
  % pos.b = meanb-1;
  % pos.theta = pi/6;
  % 
  % clf;
  % 
  % tmp = ellipsepixels([pos.x,pos.y,pos.a*4,pos.b*4,pos.theta],[1,nr,1,nc]);
  % im = zeros(size(tmp));
  % im(tmp) = rand([nnz(tmp),1]);
  % 
  % subplot(1,4,1);
  % imagesc(im,[0,1]); axis image;

  % create the bin template image

  chunk_size=100;
  for i1=1:nflies
    x0{i1}=trx(flies(i1)).x;
    y0{i1}=trx(flies(i1)).y;
    theta0{i1}=trx(flies(i1)).theta;
    a0{i1}=trx(flies(i1)).a;
    b0{i1}=trx(flies(i1)).b;
    dt0{i1}=trx(flies(i1)).dt;
    firstframes=trx.firstframes;
    endframes=trx.endframes;
  end
  parfor_tmp=cell(1,ceil((max(endframes)-min(firstframes)+1)/chunk_size));

  parfor chunk=1:ceil((max(endframes)-min(firstframes)+1)/chunk_size)
    frame_from = min(firstframes) + (chunk-1)*chunk_size;
    frame_to = min(max(endframes), frame_from+chunk_size);
    for i1=1:nflies
      parfor_tmp{chunk}{i1}=cell(1, frame_to-frame_from);
    end
    imnorm=nan([nflies size(binidx{1})]);
    imnorm_last=[];
    [readframe,nframes,fid,headerinfo] = get_readframe_fcn(fullfile(trx.expdirs{n},trx.moviefilestr));
    
    for framei = frame_from : frame_to
      disp(['frame ' num2str(framei) ', ' num2str(100*(framei-min(firstframes))/(max(endframes)-min(firstframes)),3) '%']);
      im=readframe(framei);

      for i1 = 1:nflies,
        if ((framei<firstframes(i1)) || (framei>endframes(i1)))  continue;  end

        fly1 = flies(i1);

        % extract out image region

        x=x0{fly1}(framei-firstframes(i1)+1);
        y=y0{fly1}(framei-firstframes(i1)+1);
        theta=theta0{fly1}(framei-firstframes(i1)+1);
        a=a0{fly1}(framei-firstframes(i1)+1);
        b=b0{fly1}(framei-firstframes(i1)+1);

        imnorm(i1,:,:) = compute_spacetime_transform(im,x,y,theta,a,b,meana,meanb);

        % example: average value within each box

        if ((framei==frame_from) || (framei==firstframes(i1)))  continue;  end

        parfor_tmp{chunk}{i1}{framei-frame_from} = ...
            compute_spacetime_gradient(imnorm(i1,:,:),imnorm_last(i1,:,:),binidx,nbins,featurenames,...
            dt0{fly1}(framei-firstframes(i1)));
      end

      imnorm_last=imnorm;
    end
    try  % could also set fid=-1 for .seq in get_readframe_fcn and then if(fid~=-1) here...
      fclose(fid);
    catch
    end
  end
  
  tmp=cell(1,nflies);
  for flyi=1:nflies
    cellfun(@(x) x{flyi}, parfor_tmp, 'uniformoutput',false);
    tmp{flyi}=[ans{:}];
  end

  data=cell(1,nflies);
  for flyi=1:nflies
    data{flyi}=cell(1,length(featurenames));
    for maski=1:length(featurenames)
      cellfun(@(x) x{maski}, tmp{flyi}, 'uniformoutput',false);
      data{flyi}{maski}=[ans{:}];
    end
  end

  units = parseunits('??/s');
  
  save(fullfile(trx.expdirs{n},trx.perframedir,'spacetime.mat'),'data','units','featurenames');
end

if ~strcmp(type_in,'_difference')
  idx=[];  idx2=0;
  while isempty(idx) && (idx2<length(featurenames))
    idx2 = idx2+1;
    idx=find(strcmp(featurenames{idx2},['t' num2str(theta_in) '_r' num2str(r_in) type_in]));
  end
  for i=1:length(data)
    data{i}=data{i}{idx2}(idx,:);
  end
else
  theta_in=num2str(theta_in);
  idx=[];  idx2=0;
  while isempty(idx) && (idx2<length(featurenames))
    idx2 = idx2+1;
    idx=find(strcmp(featurenames{idx2},['t' theta_in(1) '_r' num2str(r_in)]));
  end
  idxb=[];  idx2b=0;
  while isempty(idxb) && (idx2b<length(featurenames))
    idx2b = idx2b+1;
    idxb=find(strcmp(featurenames{idx2b},['t' theta_in(2) '_r' num2str(r_in)]));
  end
  for i=1:length(data)
    data{i} = abs(data{i}{idx2}(idx,:) - data{i}{idx2b}(idxb,:));
  end
end
