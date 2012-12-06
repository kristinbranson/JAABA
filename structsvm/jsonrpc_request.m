function response = jsonrpc_request(host, port, query)
    query.jsonrpc = '2.0';
    str = tcp_request(host, port, savejson(query));
    display(sprintf('Query: %s\nResponse: %s', savejson(query), str));
    response = loadjson(str);
end