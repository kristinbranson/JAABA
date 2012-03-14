% [x,y,theta,a,b] = GetTrxPos1(obj,expi,fly,ts)
% Returns the position for the input experiment, SINGLE fly, and
% frames. If ts is not input, then all frames are returned.
function pos = JLabelData_GetTrxPos(obj,expi,fly,ts)

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
    
  case 'larvacontour',
    
    
    if nargin < 4,
      pos.x = obj.trx(fly).x;
      pos.y = obj.trx(fly).y;
      pos.theta = obj.trx(fly).theta;
      pos.a = obj.trx(fly).a;
      pos.b = obj.trx(fly).b;
      pos.xcontour = obj.trx(fly).xcontour;
      pos.ycontour = obj.trx(fly).ycontour;
      return;
    end
    
    pos.x = obj.trx(fly).x(ts + obj.trx(fly).off);
    pos.y = obj.trx(fly).y(ts + obj.trx(fly).off);
    pos.theta = obj.trx(fly).theta(ts + obj.trx(fly).off);
    pos.a = obj.trx(fly).a(ts + obj.trx(fly).off);
    pos.b = obj.trx(fly).b(ts + obj.trx(fly).off);
    pos.xcontour = obj.trx(fly).xcontour{ts + obj.trx(fly).off};
    pos.ycontour = obj.trx(fly).ycontour{ts + obj.trx(fly).off};
    
  case 'larvasamuel',
    
    
    if nargin < 4,
      pos.x = obj.trx(fly).x;
      pos.y = obj.trx(fly).y;
      pos.theta = obj.trx(fly).theta;
      pos.a = obj.trx(fly).a;
      pos.b = obj.trx(fly).b;
      pos.xcontour = obj.trx(fly).xcontour;
      pos.ycontour = obj.trx(fly).ycontour;
      pos.xspine = obj.trx(fly).xspine;
      pos.yspine = obj.trx(fly).yspine;
      pos.xhead = obj.trx(fly).xhead;
      pos.yhead = obj.trx(fly).yhead;
      pos.xmid = obj.trx(fly).xmid;
      pos.ymid = obj.trx(fly).ymid;
      pos.xtail = obj.trx(fly).xtail;
      pos.ytail = obj.trx(fly).ytail;
      return;
    end
    
    is = ts + obj.trx(fly).off;
    pos.x = obj.trx(fly).x(is);
    pos.y = obj.trx(fly).y(is);
    pos.theta = obj.trx(fly).theta(is);
    pos.a = obj.trx(fly).a(is);
    pos.b = obj.trx(fly).b(is);
    pos.xcontour = obj.trx(fly).xcontour{is};
    pos.ycontour = obj.trx(fly).ycontour{is};
    pos.xspine = obj.trx(fly).xspine(:,is);
    pos.yspine = obj.trx(fly).yspine(:,is);
    pos.xhead = obj.trx(fly).xhead(is);
    pos.yhead = obj.trx(fly).yhead(is);
    pos.xmid = obj.trx(fly).xmid(is);
    pos.ymid = obj.trx(fly).ymid(is);
    pos.xtail = obj.trx(fly).xtail(is);
    pos.ytail = obj.trx(fly).ytail(is);
    
  case 'larvaspecies',
    
    if nargin < 4,
      pos.x = obj.trx(fly).x;
      pos.y = obj.trx(fly).y;
      pos.theta = obj.trx(fly).theta;
      pos.a = obj.trx(fly).a;
      pos.b = obj.trx(fly).b;
      pos.xcontour = obj.trx(fly).xcontour;
      pos.ycontour = obj.trx(fly).ycontour;
      pos.xspine = obj.trx(fly).xspine;
      pos.yspine = obj.trx(fly).yspine;
      return;
    end
    
    is = ts + obj.trx(fly).off;
    pos.x = obj.trx(fly).x(is);
    pos.y = obj.trx(fly).y(is);
    pos.theta = obj.trx(fly).theta(is);
    pos.a = obj.trx(fly).a(is);
    pos.b = obj.trx(fly).b(is);
    pos.xcontour = obj.trx(fly).xcontour{is};
    pos.ycontour = obj.trx(fly).ycontour{is};
    pos.xspine = obj.trx(fly).xspine(:,is);
    pos.yspine = obj.trx(fly).yspine(:,is);
    
end