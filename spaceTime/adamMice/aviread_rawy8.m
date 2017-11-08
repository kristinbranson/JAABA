function mov = aviread_rawy8(filename,varargin)
%AVIREAD Read AVI file.
%   AVIREAD will be removed in a future release. Use VIDEOREADER instead.
%
%   MOV = AVIREAD(FILENAME) reads the AVI movie FILENAME into the
%   MATLAB movie structure MOV. If FILENAME does not include an
%   extension, then '.avi' will be used. MOV has two fields, "cdata"
%   and "colormap". If the frames are truecolor images, MOV.cdata will
%   be Height-by-Width-by-3 and MOV.colormap will be empty. If the
%   frames are indexed images, then MOV.cdata field will be
%   Height-by-Width and MOV.colormap will be M-by-3. On UNIX, FILENAME
%   must be an uncompressed AVI file.
%
%   MOV = AVIREAD(...,INDEX) reads only the frame(s) specified by
%   INDEX.  INDEX can be a single index or an array of indices into the
%   video stream, where the first frame is number one.
%
%   The supported frame types are 8-bit (indexed or grayscale),
%   16-bit grayscale or 24-bit (TrueColor) frames.
%
%   See also VIDEOREADER, VIDEOWRITER, MOVIE.

%   Copyright 1984-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.18.2.1 $  $Date: 2011/12/13 18:36:41 $

%warning(message('MATLAB:audiovideo:aviread:FunctionToBeRemoved'));
    

% Initialization
map = [];
index = -10000;

% Validate input/output.
error(nargoutchk(1,1,nargout, 'struct'));
error(nargchk(1,2,nargin, 'struct'));

if nargin == 2
    if isnumeric(varargin{1})
        index = varargin{1};
    else
        error(message('MATLAB:audiovideo:aviread:nonNumericIndex'));
    end
end

if ~ischar(filename)
    error(message('MATLAB:audiovideo:aviread:invalidFilename'));
end

[path,name,ext] = fileparts(filename);
if isempty(ext)
    filename = strcat(filename,'.avi');
end

info = getAviInfo( filename );

if ismember(lower(info.VideoFrameHeader.CompressionType),{'y800','y8'}),
  info.VideoFrameHeader.CompressionType = 'RAW ';
end

if (index(1) ~= -10000)
    if ( isfield(info.MainHeader,'HasIndex') == 0 )
        error(message('MATLAB:audiovideo:aviread:indexParameterNotSupported', filename));
    elseif any(index > info.MainHeader.TotalFrames)
        error(message('MATLAB:audiovideo:aviread:indexOutOfRange', info.MainHeader.TotalFrames));
    elseif any(index <= 0)
        error(message('MATLAB:audiovideo:aviread:invalidZeroIndex'));
    end
end

% if ispc
%     try
%         if index == -10000
%             X = readavi(info.Filename,-1);
%         else
%             X = readavi(info.Filename,index-1);
%         end
%     catch exception
%         if strcmpi(exception.identifier,'MATLAB:audiovideo:aviread:invalidcodec') == 1 && ...
%             strcmpi(info.VideoFrameHeader.CompressionType,'indeo5')
%             % if we get an invalid codec error and are asking for indeo5
%             % inform the user of the problem and link to mathworks.com
%             error(message('MATLAB:audiovideo:aviread:indeo5NotFound'));
%         else
%             throw(exception);
%         end
%     end
% 
%     % The width of the frames may be padded to be on 4 byte boundaries because
%     % this is how bitmaps are typically written.  However, AVI files do not
%     % restrict the frames to lie on 4 byte boundries.  Because of
%     % inconsistencies with the header information in AVI files, the READAVI
%     % MEX-function determines if the frame is padded and reads the correct
%     % number of bytes for each frame.  If the returned data is more than
%     % expected (different for indexed and RGB) then the data is padded.  The
%     % padding is removed before frames are returned by AVIREAD.
% 
%     width = info.VideoFrameHeader.Width;
%     if ((info.VideoFrameHeader.Width*info.VideoFrameHeader.Height ~= length(X(1).cdata)) && (info.VideoFrameHeader.BitDepth == 8)) || ((info.VideoFrameHeader.Width*info.VideoFrameHeader.Height*3 ~= length(X(1).cdata)) && (info.VideoFrameHeader.BitDepth == 24) )
%         paddedWidth = 4*ceil(width/4);
%     else
%         paddedWidth = width;
%     end
%     height = info.VideoFrameHeader.Height;
% 
%     if info.VideoFrameHeader.BitDepth == 8 || info.VideoFrameHeader.BitDepth == 16
%         % Indexed frames
%         for i=1:length(X)
%             if height<0
%                 % Movie frames are stored top down (height is negative).
%                 mov(i).cdata = reshape(X(i).cdata, paddedWidth, abs(height))';
%             else
%                 mov(i).cdata = rot90(reshape(X(i).cdata, paddedWidth, height));
%             end
%             if paddedWidth ~= width
%                 mov(i).cdata = mov(i).cdata(:,1:width);
%             end
%             map = reshape(X(i).colormap, 4, ...
%                 info.VideoFrameHeader.NumColormapEntries);
%             map = double(flipud(map(1:3,:))')/255;
%             mov(i).colormap = map;
%         end
%     elseif info.VideoFrameHeader.BitDepth == 24
%         paddedWidth = 4*ceil(width*3/4);
%         % Truecolor frames
%         for i=1:length(X)
%             f = X(i).cdata;
%             if height<0
%                 % Movie frames are stored top down (height is negative).
%                 f = permute(reshape(f, paddedWidth, abs(height)),[2 1 3]);
%             else
%                 f = rot90(reshape(f, paddedWidth,height));
%             end
%             if paddedWidth ~= width
%                 f = f(:,1:width*3);
%             end
%             % Movie frames are stored top down (height is negative).
%             RGB(1:abs(height), 1:width,3) = f(:,1:3:end);
%             RGB(:, :, 2) = f(:,2:3:end);
%             RGB(:, :, 1) = f(:,3:3:end);
%             mov(i).cdata = RGB;
%             mov(i).colormap = [];
%         end
%     else
%         error(message('MATLAB:audiovideo:aviread:invalidColorBitDepth'));
%     end
% end %End of PC specific code

% if isunix

    if isempty(strmatch(lower(info.VideoFrameHeader.CompressionType),...
            {'dib ', 'raw ','none','raw ',char([0 0 0 0])}))
        error(message('MATLAB:audiovideo:aviread:compressedFileInvalid'));
    end

    if strcmpi(info.VideoFrameHeader.CompressionType,char([0 0 0 0]))
        info.VideoFrameHeader.CompressionType = 'none';
    end

    fid = fopen(filename,'r','l');
    if fid == -1
        error(message('MATLAB:audiovideo:aviread:unableToOpenFile', filename));
    end

    % Find RIFF chunk
    [chunk, msg,msgID] = findchunk(fid,'RIFF');
    errorWithFileClose(msgID,msg,fid);

    % Read AVI chunk
    [rifftype,msg,msgID] = readfourcc(fid);
    errorWithFileClose(msgID,msg,fid);
    if ( strcmp(rifftype,'AVI ') == 0 )
        errorWithFileClose('MATLAB:audiovideo:aviread:MissingAVIChunk',xlate('Not a valid AVI file. Missing ''AVI '' chunk.'), fid);
    end

    % Find hdrl LIST chunk
    [hdrlsize,msg,msgID] = findlist(fid,'hdrl');
    errorWithFileClose(msgID,msg,fid);

    % Find and skip avih chunk
    [chunk,msg,msgID] = findchunk(fid,'avih');
    errorWithFileClose(msgID,msg,fid);
    [msg,msgID] = skipchunk(fid,chunk);
    errorWithFileClose(msgID,msg,fid);

    % Find the video stream
    for  i = 1:info.MainHeader.NumStreams
        % Find strl LIST chunk
        [strlsize,msg,msgID] = findlist(fid,'strl');
        errorWithFileClose(msgID,msg,fid);
        % Read strh chunk
        [strhchunk,msg,msgID] = findchunk(fid,'strh');
        errorWithFileClose(msgID,msg,fid);
        % Determine stream type
        streamType = readfourcc(fid);
        % Break if it is a video stream
        if(strcmp(streamType,'vids'))
            found = 1;
            break;
        else
            found  = 0;
            % Seek to end of strl list minus the amount read
            if ( fseek(fid,listsize - 16,0) == -1 )
                errorWithFileClose('MATLAB:audiovideo:aviread:invalidChunkSize',getString(message('MATLAB:audiovideo:aviread:invalidChunkSizeAVI')), fid);
            end
        end
    end

    if (found == 0)
        errorWithFileClose('MATLAB:audiovideo:aviread:noVideoStreamFound',getString(message('MATLAB:audiovideo:aviread:noVideoStreamFound')), fid);
    end

    % Skip the strh chunk minus the fourcc (4 bytes) already read.
    strhchunk.cksize = strhchunk.cksize - 4;
    [msg,msgID] = skipchunk(fid,strhchunk);
    errorWithFileClose(msgID,msg,fid);

    % Read strf chunk
    [strfchunk,msg,msgID] = findchunk(fid,'strf');
    errorWithFileClose(msgID,msg,fid);

    if info.VideoFrameHeader.BitDepth == 24
        % For TrueColor images, skip the Stream Format chunk
        [msg,msgID] = skipchunk(fid,strfchunk);
        errorWithFileClose(msgID,msg,fid);
    elseif  info.VideoFrameHeader.BitDepth == 8
        % If bitmap has a palette Seek past the BITMAPINFOHEADER to put the
        % file pointer at the beginning of the colormap
        if  fseek(fid,info.VideoFrameHeader.BitmapHeaderSize,0) == -1
            errorWithFileClose('MATLAB:audiovideo:aviread:invalidBitmapInfoHeaderSize',getString(message('MATLAB:audiovideo:aviread:invalidBitmapInfoHeaderSize')),fid);
        end
        map = readColormap(fid,info.VideoFrameHeader.NumColorsUsed);
    else
        errorWithFileClose('MATLAB:audiovideo:aviread:invalidColorBitDepth',getString(message('MATLAB:audiovideo:aviread:invalidColorBitDepthBitmap')), fid);
    end

    % Search for the movi LIST
    [movisize,msg,msgID] = findlist(fid,'movi');
    errorWithFileClose(msgID,msg,fid);
    % movioffset will be used when using idx1. The offsets stored in idx1 are
    % with respect to just after the 'movi' LIST (not including 'movi')
    movioffset = ftell(fid) -4;

    %Determine method of reading movie
    if ( index ~= -10000 )
        method = 'UseIndex';
    elseif ( isfield(info,'IsInterleaved') == 1 )
        method = 'Interleaved';
    else
        method = 'Normal';
    end

    totalReadFrames = 1;
    switch method
        case 'Interleaved'
            % Read movies that are interleaved (contain rec lists)
            while(totalReadFrames < info.MainHeader.TotalFrames )
                [recsize, msg] = findlist(fid,'rec ');
                currentPos = ftell(fid);
                % The recListPos is the current position plus the rec list size minus
                % 4 bytes because we read the four character LIST name in findlist
                recListEndPos = currentPos + recsize - 4;
                chunk = readchunk(fid);
                while(~strcmpi(chunk.ckid,'00db') && ~strcmpi(chunk.ckid,'00dc'))
                    [msg,msgID] = skipchunk(fid,chunk);
                    errorWithFileClose(msgID,msg,fid);
                    if (ftell(fid)==recListEndPos)
                        [recsize,msg,msgID] = findlist(fid,'rec ');
                        errorWithFileClose(msgID,msg,fid);
                        currentPos = ftell(fid);
                        recListEndPos = currentPos + recsize - 4; % minus 4 because we read
                    end
                    chunk = readchunk(fid);
                end

                % Prepare input for readbmpdata
                tempinfo = info;
                tempinfo.ImageDataOffset = ftell(fid);
                tempinfo.CompressionType = info.VideoFrameHeader.CompressionType;
                tempinfo.BitDepth = info.VideoFrameHeader.BitDepth;
                % readbmpdata opens and closes the file so aviread must also to
                % have a valid fid.
                status = fclose(fid);

                RGB = readbmpdata(tempinfo);

                fid = fopen(filename,'r','l');
                fseek(fid,tempinfo.ImageDataOffset,'bof');
                % readbmpdata does not move the file pointer
                [msg,msgID] = skipchunk(fid,chunk);
                errorWithFileClose(msgID,msg,fid);
                % Assign to MATLAB movie structure
                mov(totalReadFrames).cdata = RGB;
                totalReadFrames = totalReadFrames+1;
            end
        case 'Normal'
            % Find each of the 00db or 00dc chunks and read the frames
            for i = 1:info.MainHeader.TotalFrames
                chunk = readchunk(fid);
                while (~strcmpi(chunk.ckid,'00db') && ~strcmpi(chunk.ckid,'00dc') );
                    [msg,msgID] = skipchunk(fid,chunk);
                    errorWithFileClose(msgID,msg,fid);
                    chunk = readchunk(fid);
                end

                % Prepare data to be sent to readbmpdata
                tempinfo.Filename = info.Filename;
                tempinfo.Width = info.VideoFrameHeader.Width;
                tempinfo.Height = info.VideoFrameHeader.Height;
                tempinfo.ImageDataOffset = ftell(fid);
                tempinfo.CompressionType = 'none';
                tempinfo.BitDepth = info.VideoFrameHeader.BitDepth;
                % readbmpdata opens and closes the file so aviread must also to
                % have a valid fid
                status = fclose(fid);
                % Read RGB frame
                RGB = readbmpdata(tempinfo);

                fid = fopen(filename,'r','l');
                fseek(fid,tempinfo.ImageDataOffset,'bof');
                [msg,msgID] = skipchunk(fid,chunk);
                errorWithFileClose(msgID,msg,fid);
                % Assign to MATLAB movie structure
                mov(i).cdata = RGB;
            end
        case 'UseIndex'
            % Skip the movi LIST (minus 4 because 'movi' was read) and use idx1.
            if ( fseek(fid,movisize-4,0) == -1 )
                errorWithFileClose('MATLAB:audiovideo:aviread:InvalidChunkSize',getString(message('MATLAB:audiovideo:aviread:InvalidChunkSizeWAV')), fid);
            end
            % Find idx1 chunk
            [idx1chunk,msg,msgID] = findchunk(fid,'idx1');
            errorWithFileClose(msgID,msg,fid);
            idx1ChunkPos = ftell(fid);

            for j = 1:length(index)
                fseek(fid,idx1ChunkPos,'bof');
                for i = 1:index(j)
                    found = 0;
                    while(found == 0)
                        id = readfourcc(fid);
                        if (strcmpi(id,'00db') || strcmpi(id,'00dc'))
                            found = 1;
                        end
                        [idx1data,msg,msgID] = readIDX1(fid);
                        errorWithFileClose(msgID,msg,fid);
                    end

                    % If the very first index inside the 'idx1' chunk contains a
                    % data offset that is larger than the 'movi' list offset, then
                    % it is an absolute offset.  Otherwise, we will calculate the
                    % offset relative to the beginning of the 'movi' list.
                    if (i == 1)
                        if (idx1data.offset > movioffset)
                            isIdx1OffsetAbsolute = true;
                        else
                            isIdx1OffsetAbsolute = false;
                        end
                    end
                end

                % Prepare data to be sent to readbmpdata
                tempinfo.Filename = info.Filename;
                tempinfo.Width = info.VideoFrameHeader.Width;
                tempinfo.Height = info.VideoFrameHeader.Height;
                tempinfo.BitDepth = info.VideoFrameHeader.BitDepth;
                % 8 is the riffheadersize

                % Add the 'movi' offset only if the data offset is not absolute.
                if isIdx1OffsetAbsolute
                    tempinfo.ImageDataOffset = idx1data.offset + 8;
                else
                    tempinfo.ImageDataOffset = movioffset + idx1data.offset + 8;
                end
                %tempinfo.CompressionType = info.VideoFrameHeader.CompressionType;
                tempinfo.CompressionType = 'none';
                % readbmpdata opens and closes the file so aviread must also to
                % have a valid fid
                currPos = ftell(fid);
                status = fclose(fid);

                frame = readbmpdata(tempinfo);

                fid = fopen(filename,'r','l');
                fseek(fid,currPos,'bof');
                mov(j).cdata = frame;
            end
    end
    fclose(fid);
    % Formulate outputs
    [mov.colormap] = deal(map);
    varargout{1} = mov;
%end
return

function map = readColormap(fid,numColors)
% Read colormap for 8-bit indexed images
map = fread(fid,numColors*4,'*uint8');
map = reshape(map,4,numColors);
map = double(flipud(map(1:3,:))')/255;
return;


function [idx1data,msg,msgID] = readIDX1(fid)
% Read the data in the idx1 chunk.
msg = '';
msgID='';
[idx1data.flags, count] = fread(fid,1,'uint32');
if ( count ~= 1 )
    msgID = 'MATLAB:audiovideo:aviread:IncorrectIDX1ChunkSize';
    msg = getString(message(msgID));
end
[idx1data.offset, count] = fread(fid,1,'uint32');
if ( count ~= 1 )
    msgID = 'MATLAB:audiovideo:aviread:IncorrectIDX1ChunkSize';
    msg = getString(message(msgID));
end
[idx1data.length, count] = fread(fid,1,'uint32');
if ( count ~= 1 )
    msgID = 'MATLAB:audiovideo:aviread:IncorrectIDX1ChunkSize';
    msg = getString(message(msgID));    
end
return;

function errorWithFileClose(msgID,msg,fid)
%Close open file the error
if ~isempty(msgID)
    fclose(fid);
    error(msgID,msg);
end
return;

function info = getAviInfo( filename )
    aviinfoWarningID = 'MATLAB:audiovideo:aviinfo:FunctionToBeRemoved';
    S = warning('OFF', aviinfoWarningID);
    cleaner = onCleanup(@()warning(S));
    
    % Save the state of lastwarn
    [preWarningMsg, preWarningID] = lastwarn(); 
    
    info = aviinfo(filename,'Robust');
    
    [~, lastWarningID] = lastwarn();
    
    % If aviinfo warns about deprecation restore 
    % the penultimate last warning
    if (strcmp(lastWarningID, aviinfoWarningID))
        lastwarn( preWarningMsg, preWarningID );
    end
    

