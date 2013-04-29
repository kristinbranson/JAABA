function response = jsonrpc_request_classify_example(host, port, x, imagename)
%  response=jsonrpc_classify_example(host, port, filename, plot_type)   
%     Generate a training/test error decomposition plot, querying
%     statistics from an online structured learner over the network
%  PARAMETERS:
%     host: the host IP or domain name (usually set this to 'localhost'
%     port: the socket the server is listening on
%     x: JSON string encoding the example data
%     imagename: optional name of the image file to save a visualization of
%       the prediction to
    query = {};
    query.method = 'classify_example';
    query.x = x;
    if nargin >= 4, query.visualization = imagename; end
    response = jsonrpc_request(host, port, query);
    if nargin >= 4, imshow(imagename); end
end