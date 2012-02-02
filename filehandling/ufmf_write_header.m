function [indexlocpr] = ufmf_write_header(fp,version,boxw,boxh,coding)

structtype.s = 'char';
structtype.I = 'uint32';
structtype.Q = '??';
structtype.H = '??';
structtype.B = '??';

fwrite(fp,'ufmf',structtype.s);
fwrite(fp,version,structtype.I);
% don't know indexloc yet
indexlocpr = ftell(fp);
fwrite(fp,0,structtype.Q); 
fwrite(fp,boxw,structtype.H);
fwrite(fp,boxh,structtype.H);
% coding is the encoding of each bit, e.g. MONO8
coding_str_len = length(coding);
fwrite(fp,coding_str_len,structtype.B);
fwrite(fp,coding,structtype.s);