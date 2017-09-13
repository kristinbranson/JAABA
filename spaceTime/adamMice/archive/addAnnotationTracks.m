function addAnnotationTracks(infile)

Q = load(infile);
Q.trx(2) = Q.trx(1);
Q.trx(3) = Q.trx(1);
Q.trx(1).x(:) = Q.trx(1).arena.food(1);
Q.trx(1).y(:) = Q.trx(1).arena.food(2);
Q.trx(2).x(:) = Q.trx(1).arena.mouth(1);
Q.trx(2).y(:) = Q.trx(1).arena.mouth(2);
Q.trx(3).x(:) = Q.trx(1).arena.perch(1);
Q.trx(3).y(:) = Q.trx(1).arena.perch(2);

save(infile,'-struct','Q');