
Example training command:
  ./svm_struct_learn -c 1000 -e 5 -F ../fly_data/Params/WalkSVMFeatureParams.txt -B ../fly_data/Params/WalkSVMBehaviorParams.txt ../fly_data/TrainingData/WalkTrainList.txt ../fly_data/Model/WalkSVMStructLearnedModel.txt

Example texting command:
  ./svm_struct_test -F ../fly_data/Params/WalkSVMFeatureParams.txt -B ../fly_data/Params/WalkSVMBehaviorParams.txt -O ../fly_data/TrainingData/WalkTestList.txt ../fly_data/Model/WalkSVMStructLearnedModel.txt
  
