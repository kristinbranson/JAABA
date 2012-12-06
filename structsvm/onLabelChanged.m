function response = onLabelChanged(jdata, expi, fly_id)
    params = struct('method', 'relabel_example',...
        'x', getBoutData(jdata, expi, fly_id), ...
        'y', getBoutLabel(jdata, expi, fly_id));
    response = jsonrpc_request(jdata.h.ip_address, jdata.h.port, params);
end

