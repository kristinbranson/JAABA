classdef ExpPPLaserEvent
% ExpPPLaserEvent
% Laser-pulse event for behavioral analysis. Laser interrupts a behavior
% and can cause behavioral progression.

  properties
    laserOn = nan; % frame laser is turned on
    
    intrptNbeh = 0; % number of behavior bouts real-interrupted
    
    intrptBeh0Name = ''; % name of interrupted behavior bout 
    intrptBeh0t0 = nan; % start of interrupted behavior bout
    intrptBeh0IsPhantom = false; % if true, interrupted behavior is phantom-interrupted
    
    intrptTFProgress = false; % if true, progression occurred
    intrptBehProgName = ''; % name of progression behavior
    intrptBehProgt0 = nan; % start of progression behavior
  end
  properties (Dependent)  
    intrptBeh0tdel % elapsed time from start of interrupted behavior bout to laserOn
  end
  
  methods 
    function v = get.intrptBeh0tdel(obj)
      v = obj.laserOn - obj.intrptBeh0t0;
    end
  end
    
end

