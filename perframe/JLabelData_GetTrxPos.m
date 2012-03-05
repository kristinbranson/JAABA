% [x,y,theta,a,b] = GetTrxPos1(obj,expi,fly,ts)
% Returns the position for the input experiment, SINGLE fly, and
% frames. If ts is not input, then all frames are returned.
function pos = GetTrxPos1(obj,expi,fly,ts)

pos = struct;

if all(expi ~= obj.expi),
  % TODO: generalize to multiple flies
  [success,msg] = obj.PreLoad(expi,fly);
  if ~success,
    error('Error loading trx for experiment %d: %s',expi,msg);
  end
end

switch obj.targettype,
  
  case 'fly',
    
    if nargin < 4,
      pos.x = obj.trx(fly).x;
      pos.y = obj.trx(fly).y;
      pos.theta = obj.trx(fly).theta;
      pos.a = obj.trx(fly).a;
      pos.b = obj.trx(fly).b;
      return;
    end
    
    pos.x = obj.trx(fly).x(ts + obj.trx(fly).off);
    pos.y = obj.trx(fly).y(ts + obj.trx(fly).off);
    pos.theta = obj.trx(fly).theta(ts + obj.trx(fly).off);
    pos.a = obj.trx(fly).a(ts + obj.trx(fly).off);
    pos.b = obj.trx(fly).b(ts + obj.trx(fly).off);
    
  case 'larva',
    
    if nargin < 4,
      pos.x = obj.trx(fly).x;
      pos.y = obj.trx(fly).y;
      pos.skeletonx = obj.trx(fly).skeletonx;
      pos.skeletony = obj.trx(fly).skeletony;
      return;
    end
    
    pos.x = obj.trx(fly).x(ts + obj.trx(fly).off);
    pos.y = obj.trx(fly).y(ts + obj.trx(fly).off);
    pos.skeletonx = obj.trx(fly).skeletonx(:,ts + obj.trx(fly).off);
    pos.skeletony = obj.trx(fly).skeletony(:,ts + obj.trx(fly).off);

  case 'wingedfly',
    
    if nargin < 4,
      pos.x = obj.trx(fly).x;
      pos.y = obj.trx(fly).y;
      pos.theta = obj.trx(fly).theta;
      pos.a = obj.trx(fly).a;
      pos.b = obj.trx(fly).b;
      pos.xwingl = obj.trx(fly).xwingl;
      pos.ywingl = obj.trx(fly).ywingl;
      pos.xwingr = obj.trx(fly).xwingr;
      pos.ywingr = obj.trx(fly).ywingr;
      return;
    end
    
    pos.x = obj.trx(fly).x(ts + obj.trx(fly).off);
    pos.y = obj.trx(fly).y(ts + obj.trx(fly).off);
    pos.theta = obj.trx(fly).theta(ts + obj.trx(fly).off);
    pos.a = obj.trx(fly).a(ts + obj.trx(fly).off);
    pos.b = obj.trx(fly).b(ts + obj.trx(fly).off);
    pos.xwingl = obj.trx(fly).xwingl(ts + obj.trx(fly).off);
    pos.ywingl = obj.trx(fly).ywingl(ts + obj.trx(fly).off);
    pos.xwingr = obj.trx(fly).xwingr(ts + obj.trx(fly).off);
    pos.ywingr = obj.trx(fly).ywingr(ts + obj.trx(fly).off);

end