function [nanimals, tracklens] = trxStats(trxfile)
% function [nanimals, tracklens] = trxStats(trxfile)
% Reads a trx file and outputs the number of trajectories in the file and
% track length of each trajectory.

Q = load_tracks(trxfile);
nanimals = numel(Q);
tracklens = zeros(1,nanimals);

for ndx = 1:nanimals
  tracklens(ndx) = Q(ndx).nframes;
end
