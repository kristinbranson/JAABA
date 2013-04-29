function response = jsonrpc_request_visualize(host, port, htmlDir, visualize_type, maxExamples)
%  response=jsonrpc_request_visualize(host, port, htmlDir, visualize_type, maxExamples)   
%     Generate an html visualization of the training set, or of the hardest
%     examplsin the training set
%  PARAMETERS:
%     host: the host IP or domain name (usually set this to 'localhost'
%     port: the socket the server is listening on
%     htmlDir: the directory where generated html is stored
%     visualize_type: 'train' (visualize the whole training set), 'slack' 
%         (order examples based on average slack during training), 'alpha' 
%         (order examples by dual parameter weight)
%     maxExamples: the maximum number of training examples to visualize
	
    query = {};
    query.method = 'visualize';
    query.htmlDir = htmlDir;
    if nargin >= 4, query.visualize = visualize_type; end
    if nargin >= 5, query.maxExamples = maxExamples; end
    response = jsonrpc_request(host, port, query);
end