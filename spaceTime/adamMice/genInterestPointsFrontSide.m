function [success ips] = genInterestPointsFrontSide(edir)

assert(isFrontSideExpDir(edir));

moviefrt = fullfile(edir,'movie_frt.avi');
moviesde = fullfile(edir,'movie_sde.avi');

[tf,ipssde] = genInterestPoints(moviesde,'ipnames',{'food','mouth','perch'},...
  'doaskleftright',true);
if ~tf
  success = false;
  ips = [];
  return;
end
[tf,ipsfrt] = genInterestPoints(moviefrt,'ipnames',{'food','mouth'},...
  'doaskleftright',false);
if ~tf
  success = false;
  ips = [];
  return;
end

success = true;
ips.food = ipssde.food;
ips.mouth = ipssde.mouth;
ips.perch = ipssde.perch;
ips.foodfront = ipsfrt.food;
ips.mouthfront = ipsfrt.mouth;
ips.face = ipssde.face;
