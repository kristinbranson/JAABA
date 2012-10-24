function trx = readTrx(trxname)
% format:
% version (double)
% moviename_length (int)
% moviename (char*)
% matname_length (int)
% matname (char*)
% nflies (int)
% fly ids (int*)
% firstframe (int)
% lastframe (int)
% sex (M/F) (char)
% fps (double)
% nfields (int)
% fieldname1_length (int)
% fieldname1 (char*)
% units_numerator1_length (int)
% units_numerator1 (char*)
% units_denominator1_length (int)
% units_denominator1 (char*)
% trx.fieldname1(firstframe) (double)
% ...
% trx.fieldname1(endframe) (double)
% ...
% fieldnamen_length (int)
% fieldnamen (char*)
% units_numeratorn_length (int)
% units_numeratorn (char*)
% units_denominatorn_length (int)
% units_denominatorn (char*)
% trx.fieldnamen(firstframe) (double)
% ...
% trx.fieldnamen(endframe) (double)
              
fid = fopen(trxname,'rb');
if fid < 0,
  error('Could not open file %s for reading',trxname);
end

version = fread(fid,1,'double');
moviename = freadstring(fid);
matname = freadstring(fid);
nflies = fread(fid,1,'int');
flies = fread(fid,nflies,'int');
t0 = fread(fid,1,'int');
t1 = fread(fid,1,'int');
nframes = t1-t0+1;
sex = fread(fid,1,'*char');
fps = fread(fid,1,'double');
nfields = fread(fid,1,'int');
fields = zeros(nframes,nfields);
fields = double(fields);
for fn=1:nfields
    fieldnames{fn} = freadstring(fid);
    units.num = freadstring(fid);
    untis.den = freadstring(fid);
    for frm=1:nframes
        fields(frm,fn) = fread(fid,1,'double');
    end
end
  
fclose(fid);

trx.version = version;
trx.moviename = moviename;
trx.matname = matname; 
trx.flies = flies;
trx.t0 = t0;
trx.t1 = t1;
trx.sex = sex;
trx.fps = fps;
trx.nfields = nfields;
trx.fieldnames = fieldnames;
trx.fields = fields;

function s = freadstring(fid)

l = double(fread(fid,1,'int'));
s = char(fread(fid,l,'char')');
