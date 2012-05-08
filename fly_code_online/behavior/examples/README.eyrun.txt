=========================================================================
Training and testing
=========================================================================

Synthetic data
** train command:
bin/debug_static/svm_fly_behavior_sequence.out -B Data/synth_flies/Params/LungeWalkSVMBehaviorParams.txt -F Data/synth_flies/Params/LungeWalkSVMFeatureParams.txt -d Data/synth_flies/TrainingData/traindata.txt -o Data/synth_flies/Model/LungeWalkSVMStructLearntModel.txt

** test command:
bin/debug_static/svm_fly_behavior_sequence.out -B Data/synth_flies/Params/LungeWalkSVMBehaviorParams.txt -F Data/synth_flies/Params/LungeWalkSVMFeatureParams.txt -t Data/synth_flies/TrainingData/testdata.txt Data/synth_flies/TrainingData/testdata_prediction.txt -i Data/synth_flies/Model/LungeWalkSVMStructLearntModel.txt

-------------------------------------------------------------------------

Boy meets boy data
** train command:
bin/debug_static/svm_fly_behavior_sequence.out -B Data/midres_flies/Params/BoyMeetsBoySVMBehaviorParams_morebehs.txt -F Data/midres_flies/Params/BoyMeetsBoySVMFeatureParams.txt -d Data/midres_flies/TrainingData_morebehs/traindata_5.txt -o Data/midres_flies/Model/BoyMeetsBoySVMStructLearntModel.txt

** test command:
bin/debug_static/svm_fly_behavior_sequence.out -B Data/midres_flies/Params/BoyMeetsBoySVMBehaviorParams_morebehs.txt -F Data/midres_flies/Params/BoyMeetsBoySVMFeatureParams.txt -t Data/midres_flies/TrainingData_morebehs/traindata_5.txt Data/midres_flies/TrainingData_morebehs/testdata_prediction.txt -i Data/midres_flies/Model/BoyMeetsBoySVMStructLearntModel.txt

-------------------------------------------------------------------------

Boy meets boy subset
** train command:
bin/debug_static/svm_fly_behavior_sequence.out -B Data/midres_flies/Params/BoyMeetsBoySVMBehaviorParams_morebehs.txt -F Data/midres_flies/Params/BoyMeetsBoySVMFeatureParams.txt -d Data/midres_flies/mini_dataset_morebehs/traindata_5_mini.txt -o Data/midres_flies/Model/BoyMeetsBoySVMStructLearntModel.txt

** test command: 
bin/debug_static/svm_fly_behavior_sequence.out -B Data/midres_flies/Params/BoyMeetsBoySVMBehaviorParams_morebehs.txt -F Data/midres_flies/Params/BoyMeetsBoySVMFeatureParams.txt -t Data/midres_flies/mini_dataset_morebehs/traindata_5_mini.txt Data/midres_flies/TrainingData_morebehs/testdata_prediction.txt -i Data/midres_flies/Model/BoyMeetsBoySVMStructLearntModel.txt


=========================================================================
Comparison of ground truth to prediction
=========================================================================

**To compare ground truth to predicted output:
In Matlab, run compare_pred_gt.m
