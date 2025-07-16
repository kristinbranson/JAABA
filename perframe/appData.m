% create a dummy class so that we can pass it to APTProject_app
classdef appData < handle
  properties
    aptStruct = struct();
  end
  methods
    function obj = appData()
      return;
    end
  end
end