function frameloc = ufmf_write_frame(fp,stamp,x0,y0,w,h,val)

structtype.s = 'char';
structtype.I = 'uint32';
structtype.d = 'double';
structtype.Q = '??';
structtype.H = '??';
structtype.B = '??';

%CHUNKID = '<B', # 0 = keyframe, 1 = points
%KEYFRAME2 = '<cHHd', # (dtype, width,height,timestamp)
%POINTS1 = '<dH', # timestamp, n_pts
%POINTS2 = '<HHHH', # x0, y0, w, h

frameloc = ftell(fp);
% write the chunk id (points=1)
fwrite(fp,1,structtype.B);
% write timestamp
fwrite(fp,stamp,structtype.d);
% write point position
fwrite(fp,[x0,y0,w,h],structtype.H);
% write the region intensities
fwrite(fp,val(:),'uint8');
