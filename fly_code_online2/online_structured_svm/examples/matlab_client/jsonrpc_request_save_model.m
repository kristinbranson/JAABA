function response = jsonrpc_request_save_model(host, port, filename, cachename, saveCacheFull)
%  response=jsonrpc_request_plot_statistics(host, port, filename, plot_type)   
%     Signal an online structured learner over the network to save its
%     current model parameters
%  PARAMETERS:
%     host: the host IP or domain name (usually set this to 'localhost'
%     port: the socket the server is listening on
%     filename: the name of the file to save the model to
	
    query = {};
    query.method = 'save';
    query.filename = filename;
    query.saveCacheFull = false;
    if nargin > 3, query.saveCache = cachename; end
    if nargin > 4, query.saveCacheFull = saveCacheFull; end
    response = jsonrpc_request(host, port, query);
end