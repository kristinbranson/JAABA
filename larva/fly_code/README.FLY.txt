
Example training command:
  ./svm_struct_learn -c 1000 -e 5 -F ../fly_data/Params/WalkSVMFeatureParams.txt -B ../fly_data/Params/WalkSVMBehaviorParams.txt ../fly_data/TrainingData/WalkTrainList.txt ../fly_data/Model/WalkSVMStructLearnedModel.txt

Example texting command:
  ./svm_struct_test -F ../fly_data/Params/WalkSVMFeatureParams.txt -B ../fly_data/Params/WalkSVMBehaviorParams.txt -O ../fly_data/TrainingData/WalkTestList.txt ../fly_data/Model/WalkSVMStructLearnedModel.txt
  


-c 1000 -e 5 -F ../fly_data.sharpTurns.absdtheta_velmag.mess_here/Params/SharpturnSVMFeatureParams.absdtheta_velmag.txt -B ../fly_data.sharpTurns.absdtheta_velmag.mess_here/Params/SharpturnSVMBehaviorParams.txt ../fly_data.sharpTurns.absdtheta_velmag.mess_here/TrainingData/SharpturnTrainList.8.txt ../fly_data.sharpTurns.absdtheta_velmag.mess_here/Model/SharpturnSVMStructLearnedModel.8.txt



With Debug Output:
Training Arguments:
-c 1000 -e 5 -F ../../../data/fly_data.sharpTurns.absdtheta_velmag/Params/SharpturnSVMFeatureParams.absdtheta_velmag.txt -B ../../../data/fly_data.sharpTurns.absdtheta_velmag/Params/SharpturnSVMBehaviorParams.txt  -Dwfpm ../../../data/fly_data.sharpTurns.absdtheta_velmag/debug ../../../data/fly_data.sharpTurns.absdtheta_velmag/TrainingData/SharpturnTrainList.8.txt ../../../data/fly_data.sharpTurns.absdtheta_velmag/Model/SharpturnSVMStructLearnedModel.8.txt

Testing Arguments: 
-F ../../../data/fly_data.sharpTurns.absdtheta_velmag/Params/SharpturnSVMFeatureParams.absdtheta_velmag.txt -B ../../../data/fly_data.sharpTurns.absdtheta_velmag/Params/SharpturnSVMBehaviorParams.txt  -Dwfpm ../../../data/fly_data.sharpTurns.absdtheta_velmag/debug_classify -O ../../../data/fly_data.sharpTurns.absdtheta_velmag/TrainingData/SharpturnTrainList.8.txt ../../../data/fly_data.sharpTurns.absdtheta_velmag/Model/SharpturnSVMStructLearnedModel.8.txt