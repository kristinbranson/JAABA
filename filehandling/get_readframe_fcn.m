% [readframe,nframes,fid,headerinfo] = get_readframe_fcn(filename)
% get_readframe_fcn inputs the name of a video, opens this video, and creates a
% function that can be used for random access to frames in the video. it can
% input files of type fmf, sbfmf, and avi. if videoio is installed on your
% machine, then get_readframe_fcn will use videoio to read all types of avi
% supported by videoio. otherwise, it uses VideoReader to open any type of avi
% supported by VideoReader. Run "help VideoReader" to determine what types of avi
% files are readable through VideoReader on your OS. 
%
% Input:
% filename is the name of the video to read. The type of video is determined
% based on the extension, so be sure the extension of your movie is correct. 
%
% Outputs:
% readframe is the handle to a function which inputs a frame number, reads this
% frame from filename, and outputs the corresponding frame. 
% nframes is the number of frames in the video. For avis, this may be
% approximate, as it is computed from the duration and framerate. 
% fid is the file identifier for filename. get_readframe_fcn will open the file
% filename. you must close fid yourself when you are done using the returned
% readframe. 
% headerinfo is a struct with fields depending on the type of movie

function [readframe,nframes,fid,headerinfo] = get_readframe_fcn(filename,varargin)

% allow videoio library to be used if it is installed and on the path
%CTRAX_ISVIDEOIO = exist('videoReader','file');
CTRAX_ISVIDEOIO = false;

[~,ext] = splitext(filename);
if strcmpi(ext,'.fmf'),
  [header_size,version,nr,nc,bytes_per_chunk,nframes,data_format] = fmf_read_header(filename);
  fid = fopen(filename);
  readframe = @(f) fmfreadframe(fid,header_size+(f-1)*bytes_per_chunk,nr,nc,bytes_per_chunk,data_format);
  headerinfo = struct('header_size',header_size,'version',version,'nr',nr,'nc',nc,...
    'bytes_per_chunk',bytes_per_chunk,'nframes',nframes,'data_format',data_format,'type','fmf');
elseif strcmpi(ext,'.sbfmf'),
  [nr,nc,nframes,bgcenter,bgstd,frame2file] = sbfmf_read_header(filename);
  fid = fopen(filename);
  readframe = @(f) sbfmfreadframe(f,fid,frame2file,bgcenter);
  headerinfo = struct('nr',nr,'nc',nc,'nframes',nframes,'bgcenter',bgcenter,...
    'bgstd',bgstd,'frame2file',frame2file,'type','sbfmf');
elseif strcmpi(ext,'.ufmf'),
  headerinfo = ufmf_read_header(filename);
  readframe = ufmf_get_readframe_fcn(headerinfo,varargin{:});
  nframes = headerinfo.nframes;
  fid = headerinfo.fid;
elseif strcmpi(ext,'.mmf'),
  headerinfo = mmf_read_header(filename);
  readframe = mmf_get_readframe_fcn(headerinfo,varargin{:});
  nframes = headerinfo.nframes;
  fid = headerinfo.fid;
elseif strcmpi(ext,'.seq'),
  [indexfilename] = myparse(varargin,'indexfilename',0);
  if ischar(indexfilename),
    headerinfo = r_readseqinfo(filename,indexfilename);
  else
    headerinfo = r_readseqinfo(filename);
  end
  if numel(headerinfo.m_aiSeekPos) < headerinfo.m_iNumFrames,
    headerinfo.m_iNumFrames = numel(headerinfo.m_aiSeekPos);
  end
  headerinfo.nr = headerinfo.m_iHeight;
  headerinfo.nc = headerinfo.m_iWidth;
  headerinfo.nframes = headerinfo.m_iNumFrames;
  headerinfo.type = 'seq';
  fid = 0;
  nframes = headerinfo.m_iNumFrames;
  readframe = @(f) read_seq_frame(headerinfo,f);
else
  fid = 0;
  if CTRAX_ISVIDEOIO,
    readerobj = videoReader(filename,'preciseFrames',30,'frameTimeoutMS',5000);
    info = getinfo(readerobj);
    nframes = info.numFrames;
    seek(readerobj,0);
    seek(readerobj,1);
    readframe = @(f) videoioreadframe(readerobj,f);
    headerinfo = info;
    headerinfo.type = 'avi';
  else
    readerobj = VideoReader(filename);
    nframes = get(readerobj,'NumberOfFrames');
    if isempty(nframes),
      % approximate nframes from duration
      nframes = get(readerobj,'Duration')*get(readerobj,'FrameRate');
    end
    readframe = @(f) read(readerobj,f);
    %readframe = @(f) flipdim(read(readerobj,f),1);
    headerinfo = get(readerobj);
    headerinfo.type = 'avi';
  end
end
