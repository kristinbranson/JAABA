function synchStructuredBehaviors( jdata )
  behaviors = {};
  for i=1:numel(jdata.labelnames), 
      behavior = struct('name', jdata.labelnames(i), 'color', randi(hex2dec('FFFFFF'))); %jdata.labelcolors(i));
      if i==1, behaviors = behavior; else behaviors(i) = behavior; end
  end
  query = struct('method', 'set_behaviors', 'jsonrpc', '2.0', 'behaviors', behaviors);
  response = jsonrpc_request(jdata.h.ip_address, jdata.h.port, query);
  if isfield(response, 'error'), warning(response.error); end
end

