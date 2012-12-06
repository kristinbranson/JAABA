%
% This is an example demonstrating how to access an online structured
% learning server from a matlab client over the network.  This example is
% intended to demonstrate how to invoke different commands, as opposed to
% being directly useful for any application.  
% 
% Before running this program, you need to start the server using a command
% like:
%   bin/debug_static/svm_fly_behavior_sequence.out -P 8001 -B Data/midres_flies/Params/BoyMeetsBoySVMBehaviorParams.txt -F Data/midres_flies/Params/BoyMeetsBoySVMFeatureParams.txt -d Data/midres_flies/TrainingData/traindata_1.txt
% The server will start training a behavior detector for the behaviors
% defined with in BoyMeetsBoySVMBehaviorParams.txt, features defined in
% BoyMeetsBoySVMFeatureParams.txt, and with the initial training set
% defined in traindata_1.txt.  It will listen on port 8001 for incoming
% requests over the network.
%
% This script demonstrates how to make different kinds of requests over the
% network.  It reads in a file (traindata_2.txt) with a list of new training
% examples.  For each training example, it 
%    1) Makes a network request to ask the server to predict a behavior
%       segmentation using the current learned model
%    2) Makes a network request to ask the server to add a new training
%       example
% It also demonstrates how to ask the server to save its current model to
% disk.
%


% Add functions for a matlab JSON parser to the path
addpath('../../online_structured_svm/examples/matlab_client/');

host = '127.0.0.1';  % The server's IP address (localhost)
port = 8001;    % The port the server is running on
trainset = 'Data/midres_flies/TrainingData/traindata_2.txt';   % training set file
pause_time=.5;  % time in seconds to pause before adding a new training example

train_files = textread(trainset, '%s');
dir = fileparts(trainset);

for i=1:size(train_files,1)
    [ignore fname] = fileparts(train_files{i});
    fname = strcat(dir, '/', fname);
    
    % Use the current model to classify a new training example
    response = jsonrpc_request(host, port, loadjson(sprintf('{"method":"classify_example","x":{"fname":"%s.trx"}}', fname)))
    
    % Add the new training example
    % WARNING: when the server adds a new training example, it will append
    % an entry to its training set file on disk (e.g., traindata_1.txt)
    response = jsonrpc_request(host, port, loadjson(sprintf('{"method":"add_example","x":{"fname":"%s.trx"},"y":{"fname":"%s.label"}}', fname, fname)))
    
    pause(pause_time); 
end

% Plot a decomposition of the test error (in terms of the structured SVM
% error function), to get a sense of whether the test error is coming from
% insufficient training data, bad feature space, or insufficient
% computation time spent 
% Red Curve: average structured svm error on a new unseen training example
%            (before it has been processed)
% Blue Curve: average error on each training example (after it has been
%             iterated over multiple times)
% Green Curve: lower bound on the minimal achievable model error
% 
% If the gap between the red and blue curve is high, it suggests one would
% benefit from adding more training examples
% If the gap between the blue and green curve is high, it suggests one
% would benefit from spending more computation time
% If the raw value of the green curve is high, it suggests the feature
% space has saturated and one cannot get better training error without
% augmenting the feature space (or potentially reducing the regularization
% parameter)
response = jsonrpc_request(host, port, loadjson('{"method":"plot_stats","error_decomp_by_iteration_plot_name":"error_decomp_plot"}'))
error_decomp_plot


% Save the state of the learning algorithm, such that the learned model can
% be loaded from disk to resume training using a command like
%   bin/debug_static/svm_fly_behavior_sequence.out -P 8001 -B Data/midres_flies/Params/BoyMeetsBoySVMBehaviorParams.txt -F Data/midres_flies/Params/BoyMeetsBoySVMFeatureParams.txt -i learned_model.txt
% One can also use the learned model to classify a testset using a command
% like:
%   bin/debug_static/svm_fly_behavior_sequence.out -B Data/midres_flies/Params/BoyMeetsBoySVMBehaviorParams.txt -F Data/midres_flies/Params/BoyMeetsBoySVMFeatureParams.txt -i learned_model.txt -t Data/midres_flies/TrainingData/testdata_1.txt
response = jsonrpc_request(host, port, loadjson(sprintf('{"method":"save","filename":"learned_model.txt","savefull":true}')))



