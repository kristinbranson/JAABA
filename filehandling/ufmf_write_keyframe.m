function frameloc = ufmf_write_keyframe(fp,stamp,im,keyframe_type)

w = size(im,1);
h = size(im,2);
switch dtype(im),
 case 'uint8'
  dtype_char = 'B';
 case 'float32',
  dtype_char = 'f';
 otherwise
  error('Unsupported type %s',dtype(im));
end

structtype.s = 'char';
structtype.c = '??';
structtype.I = 'uint32';
structtype.d = 'double';
structtype.Q = '??';
structtype.H = '??';
structtype.B = '??';

%CHUNKID = '<B', # 0 = keyframe, 1 = points
% KEYFRAME2 = '<cHHd', # (dtype, width,height,timestamp)

% this is a keyframe
fwrite(fp,0,structtype.B);

% write keyframe_type
fwrite(fp,len(keyframe_type),'uint8');
fwrite(fp,keyframe_type,structtype.s);

% write dtype (one of 'B' (uint8), 'f' (float32))
fwrite(fp,dtype,structtype.c);
% width
fwrite(fp,w,structtype.H);
% height
fwrite(fp,h,structtype.H);
% stamp
fwrite(fp,stamp,structtype.H);

% write the data
fwrite(fp,im,dtype(im));