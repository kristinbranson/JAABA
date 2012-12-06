=========================================================================
Training and testing
=========================================================================

Synthetic data
** train command:
bin/debug_static/svm_fly_behavior_sequence.out -B Data/synth_flies/Params/LungeWalkSVMBehaviorParams.txt -F Data/synth_flies/Params/LungeWalkSVMFeatureParams.txt -d Data/synth_flies/TrainingData/traindata.txt -o Data/synth_flies/Model/LungeWalkSVMStructLearntModel.txt

** test command:
bin/debug_static/svm_fly_behavior_sequence.out -B Data/synth_flies/Params/LungeWalkSVMBehaviorParams.txt -F Data/synth_flies/Params/LungeWalkSVMFeatureParams.txt -t Data/synth_flies/TrainingData/testdata.txt Data/synth_flies/TrainingData/testdata_prediction.txt -i Data/synth_flies/Model/LungeWalkSVMStructLearntModel.txt

-------------------------------------------------------------------------
-P 8099 -w 100 1.5
Boy meets boy data new
** train command:
bin/debug_static/svm_fly_behavior_sequence.out -B Data/midres_flies/Params/BoyMeetsBoySVMBehaviorParams.txt -F Data/midres_flies/Params/FeatureParamsEyrunTest.txt -d Data/midres_flies/TrainingDataEyrunTest/traindata_2.txt -o Data/midres_flies/Model/TrainingDataEyrunTest.txt -P 8099 -w 100 1.5

** test command:
bin/debug_static/svm_fly_behavior_sequence.out -B Data/midres_flies/Params/BoyMeetsBoySVMBehaviorParams.txt -F Data/midres_flies/Params/FeatureParamsEyrunTest.txt -t Data/midres_flies/TrainingDataEyrunTest/traindata_2.txt Data/midres_flies/TrainingDataEyrunTest/testdata_prediction.txt -i Data/midres_flies/Model/TrainingDataEyrunTest.txt -A modelid -a iterid
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

-------------------------------------------------------------------------

Mouse data
** train command:
bin/debug_static/svm_fly_behavior_sequence.out -B Data/mice/Params/MouseBehaviorParams10.txt -F Data/mice/Params/MouseFeatureParams.txt -d Data/mice/XaviData10behs/traindata.txt -o Data/mice/Model/10behs_100cost.txt -P 8099 -w 100 1.5

** test command:
bin/debug_static/svm_fly_behavior_sequence.out -B Data/mice/Params/MouseBehaviorParams10.txt -F Data/mice/Params/MouseFeatureParams.txt -t Data/mice/XaviData10behs/testdata.txt Data/mice/XaviData10behs/testdata_prediction.txt -i Data/mice/Model/10behs_100cost.txt.13 -A 1 -a 13

-------------------------------------------------------------------------

Eric's data
** train command:
bin/debug_static/svm_fly_behavior_sequence.out -B Data/eric_flies/Params/BehaviorParams_more.txt -F Data/eric_flies/Params/FeatureParams.txt -d Data/eric_flies/many_features_and_behs/traindata_special.txt -o Data/eric_flies/Model/many_features_and_behs.txt -P 8099 -w 100 1.5

** test command:
bin/debug_static/svm_fly_behavior_sequence.out -B Data/eric_flies/Params/BehaviorParams_more.txt -F Data/eric_flies/Params/FeatureParams.txt -t Data/eric_flies/many_features_and_behs/testdata_special.txt Data/eric_flies/many_features_and_behs/testdata_prediction.txt -i Data/eric_flies/Model/many_features_and_behs.txt.13 -A 1 -a 13

-------------------------------------------------------------------------

Honeybee data
** train command:
bin/debug_static/svm_fly_behavior_sequence.out -B Data/honeybee/Params/BehaviorParams.txt -F Data/honeybee/Params/FeatureParams.txt -d Data/honeybee/honeybee_data/train_minus_1.txt -o Data/honeybee/Model/honeybee_data.txt -P 8099 -w 100 1.5

** test command:
bin/debug_static/svm_fly_behavior_sequence.out -B Data/honeybee/Params/BehaviorParams.txt -F Data/honeybee/Params/FeatureParams.txt -t Data/honeybee/honeybee_data/train_1.txt Data/honeybee/honeybee_data/testdata_prediction.txt -i Data/honeybee/Model/honeybee_data.txt.13 -A 1 -a 13

=========================================================================
Comparison of ground truth to prediction
=========================================================================

**To compare ground truth to predicted output:
In Matlab, run compare_pred_gt.m
