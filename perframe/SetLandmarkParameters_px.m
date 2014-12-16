function trx = SetLandmarkParameters_px(trx,arenatype,arenacenterx,arenacentery,...
  arenaradius)

switch lower(arenatype),
  case 'none',
    return;
  case 'circle',
    for i = 1:numel(trx),
      trx(i).arena.r = arenaradius;
      trx(i).arena.x = arenacenterx;
      trx(i).arena.y = arenacentery;
    end
  case 'rectangle',
end
