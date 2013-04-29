function response = jsonrpc_request_plot_statistics(host, port, filename, plot_type)
%  response=jsonrpc_request_plot_statistics(host, port, filename, plot_type)   
%     Generate a training/test error decomposition plot, querying
%     statistics from an online structured learner over the network
%  PARAMETERS:
%     host: the host IP or domain name (usually set this to 'localhost'
%     port: the socket the server is listening on
%     filename: the name of the file to save the plot to.  The server will
%       export this to a file filename.m (so don't include an extension on
%       the filename)
%     plot_type: should either be 'time', 'example', or 'iteration'.  This
%       defines what is plotted on the x-axis of the plot
	
    query = {};
    query.method = 'plot_stats';
    query.plot_name = filename;
    query.plot_by = plot_type;
    response = jsonrpc_request(host, port, query);
    run(filename);
end