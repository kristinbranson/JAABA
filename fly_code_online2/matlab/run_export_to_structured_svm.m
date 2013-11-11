
function run_export_to_structured_svm(dataset_id)
if dataset_id == 1
    datadir = '/scratch/Datasets/FlyPairs/Jon_scored_movies/';
    export_dir = 'data_ssvm_jon';
    train = {'movie1','movie2','movie3'};
    test = {'movie4','movie5','movie6'};
    FPS = 200;
elseif dataset_id == 2
    datadir = '/scratch/Datasets/FlyPairs/Eric_scored_movies/';
    export_dir = 'data_ssvm_eric';
    train = {'movie3','movie9','movie10','movie13','movie28'};
    test = {'movie1','movie5','movie11','movie15','movie16'};
    FPS = 30;
else
    datadir = '/scratch/Datasets/FlyPairs/Eric_courtship_movies/';
    export_dir = 'data_ssvm_eric_courtship';
    train = {'wild/movie1','wild/movie12','wild/movie13','wild/movie14', ...
    'wild/movie15','wild/movie16','wild/movie17','wild/movie18', ...
    'wild/movie19','wild/movie20','wild/movie21','hyper/movie4', ...
    'hyper/movie5','hyper/movie9','hyper/movie10','hyper/movie11'};
    test = {'wild/movie2','wild/movie3','wild/movie4','wild/movie5', ...
    'wild/movie6','wild/movie7','wild/movie8','wild/movie9', ...
    'wild/movie10','wild/movie11','hyper/movie1','hyper/movie2', ...
    'hyper/movie3','hyper/movie6','hyper/movie12'};
    FPS = 30;
end

export_to_structured_svm([], [], FPS, datadir, export_dir, [], train, test)

end