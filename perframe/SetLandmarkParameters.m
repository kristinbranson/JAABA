function trx = SetLandmarkParameters(trx,arenatype,arenacenterx,arenacentery,...
  arenaradius,arenawidth,arenaheight)

switch lower(arenatype),
  case 'none',
    return;
  case 'circle',
    for i = 1:numel(trx),
      trx(i).arena.arena_radius_mm = arenaradius;
      trx(i).arena.arena_center_mm_x = arenacenterx;
      trx(i).arena.arena_center_mm_y = arenacentery;
    end
  case 'rectangle',
    for i = 1:numel(trx),
      trx(i).arena.tl = [arenacenterx - arenawidth/2,arenacentery - arenaheight/2];
      trx(i).arena.tr = [arenacenterx + arenawidth/2,arenacentery - arenaheight/2];
      trx(i).arena.bl = [arenacenterx - arenawidth/2,arenacentery + arenaheight/2];
      trx(i).arena.br = [arenacenterx + arenawidth/2,arenacentery + arenaheight/2];
    end
end