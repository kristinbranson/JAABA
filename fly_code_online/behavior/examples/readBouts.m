
function bouts = readBouts(filename)
    fid = fopen(filename);

    % read number of features
    line = fgetl(fid);
    num_features = textscan(line,'num_features = %d'); 
    num_features = num_features{1};
    % read number of bouts
    line = fgetl(fid);
    num_bouts = textscan(line,'num_bouts = %d'); 
    num_bouts = num_bouts{1};
    % load all bout features
    bouts_vec = fscanf(fid,'%f');
    % [sequence_id, f_start, f_end, class_id, <bout_features>]
    bouts = reshape(bouts_vec, num_features+4, num_bouts)';

    fclose(fid);
end
