% absolute azimuthal angle to from fly to closest fly according to type
function [data,units] = compute_spacetime(trx,n,theta_in,r_in)

try
  load(fullfile(trx.expdirs{n},trx.perframedir,'spacetime.mat'));
  
catch
  flies = trx.exp2flies{n};
  nflies = numel(flies);
  data = cell(1,nflies);

  % parameters

  meana = prctile(cellfun(@mean,trx(:).a),90);
  meanb = prctile(cellfun(@mean,trx(:).b),90);

  [binidx, nbins, featurenames, featureboundaries, featurecenters] = compute_spacetime_mask(meana, meanb);
  
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
    for i1=1:nflies
      parfor_tmp{chunk}{i1}=nan(min(chunk_size,max(endframes)-min(firstframes)+1-chunk*chunk_size),nbins);
    end
    imnorm=nan([nflies size(binidx)]);
    imnorm_last=[];
    [readframe,nframes,fid,headerinfo] = get_readframe_fcn(fullfile(trx.expdirs{n},trx.moviefilestr));
    
    pooh=min(firstframes)+(chunk-1)*chunk_size;
    for framei=pooh:min(max(endframes),pooh+chunk_size)
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

        if ((framei==pooh) || (framei==firstframes(i1)))  continue;  end

        parfor_tmp{chunk}{i1}(framei-pooh,:) = ...
            compute_spacetime_gradient(imnorm(i1,:,:),imnorm_last(i1,:,:),binidx,nbins,...
            dt0{fly1}(framei-firstframes(i1)));
      end

      imnorm_last=imnorm;
    end
    fclose(fid);
  end
  
  data=cell(1,nflies);
  for i=1:nflies
    cellfun(@(x) transpose(x{i}), parfor_tmp,'uniformoutput',false);
    data{i}=[ans{:}]';
  end

  units = parseunits('??/s');
  
  save(fullfile(trx.expdirs{n},trx.perframedir,'spacetime.mat'),'data','units','featurenames');
end

idx=find(strcmp(featurenames,['theta' num2str(theta_in) '_r' num2str(r_in)]));
for i=1:length(data)
  data{i}=data{i}(:,idx);
end