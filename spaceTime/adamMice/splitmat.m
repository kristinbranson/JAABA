function splitmat(filename)
f=fopen([filename '.mat'], 'rb');
header=fread(f,128);
i=1;
while true
    h2=fread(f,2,'int32');
    if length(h2) < 2
        disp 'finished reading file';
        break;
    end
    if h2(2) == 0
        disp(sprintf('Found bad 0-byte size at variable #%d.', i));
        break;
    end
    fout=fopen(sprintf('%s_%d.mat', filename, i), 'wb');
    fwrite(fout, header);
    fwrite(fout, h2, 'int32');
    data = fread(f, h2(2));
    fwrite(fout, data);
    fclose(fout);
    i = i+1;
end
fclose(f);