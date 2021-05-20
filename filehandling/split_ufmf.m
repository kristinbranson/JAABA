function outfiles = split_ufmf(infile,outfilestr,chunksize,varargin)

[startframe,endframe,dooverwrite] = ...
  myparse(varargin,'startframe',1,'endframe',inf,'dooverwrite',false);

% read in the header
input_header = ufmf_read_header(infile);
endframe = min(endframe,input_header.nframes);
noutfiles = ceil((endframe-startframe+1)/chunksize);

% find the end of the header
headerstart = 0;
% loc of something after the header
headerend = min(min(input_header.mean2file),min(input_header.frame2file));

% read in the header
fseek(input_header.fid,headerstart,'bof');
headerbytes = fread(input_header.fid,headerend-headerstart);

outfiles = cell(1,noutfiles);
for outi = 1:noutfiles,
  
  t0 = startframe + (outi-1)*chunksize;
  t1 = min(endframe,t0+chunksize-1);
  
  outfiles{outi} = sprintf('%s_%02d.ufmf',outfilestr,outi);
  if exist(outfiles{outi},'file') && ~dooverwrite,
    res = input(sprintf('File %s exists. Overwrite? (y/N): ',outfiles{outi}),'s');
    if ~strcmpi(res,'y'),
      fclose(input_header.fid);
      return;
    end
  end
  
  output_fid = fopen( outfiles{outi}, 'wb' , 'ieee-le');
  if output_fid < 0 ,
    error('Unable to open output file %s', output_file_name) ;
  end
  cleaner = onCleanup(@()(fclose(output_fid))) ;
  
  % write the header
  fwrite(output_fid,headerbytes);
  
  lastmean = nan;
  
  index = struct;
  index.frame = struct;
  index.frame.loc = int64([]);
  index.frame.timestamp = input_header.timestamps(t0:t1);
  index.keyframe = struct;
  index.keyframe.mean = struct;
  index.keyframe.mean.loc = int64([]);
  index.keyframe.mean.timestamp = [];
  

  % write the frames
  for t = t0:t1,

    % write a mean
    if input_header.frame2mean(t) ~= lastmean,
      meanstart = input_header.mean2file(input_header.frame2mean(t));
      meanend = min([min(input_header.mean2file(input_header.mean2file>meanstart)),...
        min(input_header.frame2file(input_header.frame2file>meanstart)),...
        input_header.indexloc]);
      assert(~isempty(meanend) && ~isnan(meanend));
      index.keyframe.mean.loc(end+1) = ftell(output_fid);
      index.keyframe.mean.timestamp(end+1) = input_header.mean_timestamps(input_header.frame2mean(t));
      fseek(input_header.fid,meanstart,'bof');
      fwrite(output_fid,fread(input_header.fid,meanend-meanstart));
      lastmean = input_header.frame2mean(t);
    end
    
    framestart = input_header.frame2file(t);
    frameend = min([min(input_header.mean2file(input_header.mean2file>framestart)),...
        min(input_header.frame2file(input_header.frame2file>framestart)),...
        input_header.indexloc]);
    assert(~isempty(frameend) && ~isnan(frameend));
    fseek(input_header.fid,framestart,'bof');
    index.frame.loc(end+1) = ftell(output_fid);
    fwrite(output_fid,fread(input_header.fid,frameend-framestart));
  end
  
%   % write the index
%   wrapupUFMF(fidout,index,header.indexlocloc);
%   fclose(fidout);
  
  % write the index
  index_offset = ftell(output_fid) ;
  ufmf_write_struct(output_fid, index) ; 
  
  % overwrite the index location in the header
  offset = input_header.indexlocloc ;  % The offset where we should write the index offset
  fseek(output_fid, offset, 'bof') ;
  fwrite(output_fid, index_offset, 'uint64') ;  
end
