
Example training command:
  ./svm_struct_learn -c 1000 -e 5 -F ../fly_data/Params/WalkSVMFeatureParams.txt -B ../fly_data/Params/WalkSVMBehaviorParams.txt ../fly_data/TrainingData/WalkTrainList.txt ../fly_data/Model/WalkSVMStructLearnedModel.txt

Example texting command:
  ./svm_struct_test -F ../fly_data/Params/WalkSVMFeatureParams.txt -B ../fly_data/Params/WalkSVMBehaviorParams.txt -O ../fly_data/TrainingData/WalkTestList.txt ../fly_data/Model/WalkSVMStructLearnedModel.txt
  


-c 1000 -e 5 -F ../fly_data.sharpTurns.absdtheta_velmag.mess_here/Params/SharpturnSVMFeatureParams.absdtheta_velmag.txt -B ../fly_data.sharpTurns.absdtheta_velmag.mess_here/Params/SharpturnSVMBehaviorParams.txt ../fly_data.sharpTurns.absdtheta_velmag.mess_here/TrainingData/SharpturnTrainList.8.txt ../fly_data.sharpTurns.absdtheta_velmag.mess_here/Model/SharpturnSVMStructLearnedModel.8.txt

