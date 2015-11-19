function headerinfo = ReadIndexedMJPGHeader(filename,indexfilename)

headerinfo = struct('filename',filename,'indexfilename',indexfilename,...
  'nr',0,'nc',0,'nframes',0,'frame2file',[],'timestamp',[],'type','mjpg');

fid = fopen(indexfilename,'r');
if fid < 0,
  error('Could not open file %s for reading',indexfilename);
end

t = 0;
while true,
  s = fgetl(fid);
  if ~ischar(s),
    break;
  end
  s = regexp(s,'\s','split');
  s = strtrim(s);
  assert(str2double(s{1})==t);
  t = t+1;
  headerinfo.timestamp(t) = str2double(s{2});
  headerinfo.frame2file(t) = str2double(s{3});
end
fclose(fid);

headerinfo.nframes = t;

im = parsejpg8(filename,headerinfo.frame2file(1));
[headerinfo.nr,headerinfo.nc,~] = size(im);


