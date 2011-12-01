host = '127.0.0.1';  % The server's IP address (localhost)
port = 8001;    % The port the server is running on: ./svm_fly_behavior_sequence.out -P 8000
trainset = 'Data/midres_flies/TrainingData/traindata_2.txt';   % training set file
pause_time=.5;  % time in seconds to pause

train_files = textread(trainset, '%s');
dir = fileparts(trainset);

for i=1:size(train_files,1)
    [ignore fname] = fileparts(train_files{i});
    fname = strcat(dir, '/', fname);
    
    % Use the current model to classify a new training example
    response = jsonrpc_request(host, port, loadjson(sprintf('{"method":"classify_example","x":{"fname":"%s.trx"}}', fname)))
    
    % Add the new training example
    response = jsonrpc_request(host, port, loadjson(sprintf('{"method":"add_example","x":{"fname":"%s.trx"},""y":{"fname":"%s.label"}}', fname, fname)))
    
    pause(pause_time); 
end

% Save the state of the learning algorithm, such that the learned model can
% be loaded from disk to resume training
response = jsonrpc_request(host, port, loadjson(sprintf('{"method":"save","filename":"Data/midres_flies/train.tmp","savefull":true}')))