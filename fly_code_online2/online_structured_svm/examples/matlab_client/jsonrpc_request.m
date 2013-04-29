function response = jsonrpc_request(host, port, query)
    query.jsonrpc = '2.0';
    str = tcp_request(host, port, savejson(query));
    response = loadjson(str);
end